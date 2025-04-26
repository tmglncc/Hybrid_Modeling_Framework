import numpy as np
import pysindy_local as ps
import matplotlib.pyplot as plt
import os
from scipy.integrate import odeint
from sklearn.gaussian_process.kernels import ExpSineSquared, Matern, RBF, WhiteKernel, RationalQuadratic, DotProduct, ConstantKernel
from ModelPlots import ModelPlots
from ModelCalibration import ModelCalibration
from ModelSelection import ModelSelection
from DataDenoising import DataDenoising

def is_new_model(model_set, model, n_vars, precision):
	for model_element in model_set:
		flag = True
		for i in range(n_vars):
			if model_element.equations(precision = precision)[i] != model.equations(precision = precision)[i]:
				flag = False
				break

		if flag:
			return False

	return True

experiment_id = 0

# Train and test data parameters
step = 1

# Method parameters
fd_order = 2
poly_degrees = range(2, 3)
fourier_nfreqs = range(1, 2)
optimizer_method = "STLSQ+SA"
precision = 3

plot_sse = True
plot_sse_correlation = False
plot_relative_error = True
plot_Ftest = True
plot_qoi = True
plot_musig = False
plot_simulation = False
plot_derivative = False
calibration_mode = "Bayes"

stlsq_alphas = [0.001, 0.01, 0.1, 1.0, 10.0]

# Read train data
data = np.genfromtxt("data.csv", dtype = float, delimiter = ',', names = True)
t = data["frameTime"][::step]
X = np.stack((data["meanNecrotic"][::step], data["meanQuiescent"][::step], data["meanProliferative2"][::step] + data["meanProliferative1"][::step]), axis = -1)

dd = DataDenoising(X, t, ["N", "Q", "P"])
dd.plot_sma([5, 10, 20])
dd.plot_ema([0.1, 0.2, 0.3], [False])
dd.plot_l2r([10.0, 100.0, 1000.0])
dd.plot_tvr([0.001], [0.25, 0.5, 1.0])
# dd.plot_gpr([RBF(length_scale_bounds = (1.0e3, 1.0e4)),
# 	ConstantKernel(constant_value_bounds = (1.0e-2, 1.0e-1))*RBF(length_scale_bounds = (1.0e3, 1.0e4)),
# 	ConstantKernel(constant_value_bounds = (1.0e-2, 1.0e-1))*RBF(length_scale_bounds = (1.0e3, 1.0e4)) + WhiteKernel(noise_level_bounds = (1.0e-10, 1.0e-9))],
# 	[10], [1.0e-10], ["RBF", "ConstantKernel*RBF", "ConstantKernel*RBF + WhiteKernel"], 0
# )
X, X_min, X_max = dd.gaussian_process_regression(
	RBF(length_scale_bounds = (1.0e3, 1.0e4)),
	10, 1.0e-10
)

X_plot = X
t_plot = t

X = np.delete(X, slice(100), 0)
t = np.delete(t, slice(100), None)

X0 = X[0, :]
t_steps = len(t)

# Read test data
data_test = np.genfromtxt("data.csv", dtype = float, delimiter = ',', names = True)
t_test = data_test["frameTime"][::step]
X_test = np.stack((data_test["meanNecrotic"][::step], data_test["meanQuiescent"][::step], data_test["meanProliferative2"][::step] + data_test["meanProliferative1"][::step]), axis = -1)

dd = DataDenoising(X_test, t_test, ["N", "Q", "P"])
X_test, X_test_min, X_test_max = X, X_min, X_max

t_test = np.delete(t_test, slice(100), None)

X0_test = X_test[0, :]

model_set = []
for poly_degree in poly_degrees:
	for fourier_nfreq in fourier_nfreqs:
		for stlsq_alpha in stlsq_alphas:
			experiment_id += 1
			print("Experimento " + str(experiment_id) 
				+ ": Grau = " + str(poly_degree) 
				+ ", Frequência = " + str(fourier_nfreq) 
				+ ", alpha = " + str(stlsq_alpha) + "\n"
			)

			# Define method properties
			# differentiation_method = ps.FiniteDifference(order = fd_order)
			differentiation_method = ps.SmoothedFiniteDifference()
			feature_library = ps.PolynomialLibrary(degree = poly_degree) # + ps.FourierLibrary(n_frequencies = fourier_nfreq)
			optimizer = ps.STLSQ(
				alpha = stlsq_alpha,
				fit_intercept = False,
				verbose = True,
				window = 3,
				epsilon = 65.0,
				time = t,
				sa_times = np.array([100.0, 300.0, 500.0, 700.0]),
				non_physical_features = [['1', 'N', 'N^2'],
					['1', 'N', 'N^2'],
					['1', 'N', 'N^2']
				]
			)

			# Compute sparse regression
			model = ps.SINDy(
				differentiation_method = differentiation_method,
				feature_library = feature_library,
				optimizer = optimizer,
				feature_names = ["N", "Q", "P"]
			)
			model.fit(X, t = t)
			model.print(precision = precision)
			print("\n")

			# Generate model plots
			mp = ModelPlots(model, optimizer_method, experiment_id)
			if plot_sse:
				mp.plot_sse()
			if plot_sse_correlation:
				mp.plot_sse_correlation()
			if plot_relative_error:
				mp.plot_relative_error()
			if plot_Ftest:
				mp.plot_Ftest()
			if plot_qoi:
				mp.plot_qoi()
			if plot_musig:
				mp.plot_musig()
			if plot_simulation:
				mp.plot_simulation(X, t, X0, precision = precision)
			if plot_derivative:
				mp.plot_derivative(X, t)

			# Add model to the set of models
			if not model_set or is_new_model(model_set, model, len(model.feature_names), precision):
				model_set.append(model)

del model_set[1]

# Compute number of terms
ms = ModelSelection(model_set, t_steps)
ms.compute_k()

for model_id, model in enumerate(model_set):
	print("Modelo " + str(model_id+1) + "\n")
	model.print(precision = precision)
	print("\n")

	dd = DataDenoising(X_test, t_test, model.feature_names)

	# Compute derivative
	X_dot_test = model.differentiate(X_test, t_test)
	dd.plot_derivative(X_dot_test, t_test, 0, X0_test)

	# Simulate with another initial condition
	try:
		if calibration_mode is None:
			simulation = model.simulate(X0_test, t = t_test)
		elif calibration_mode == "LM":
			mc = ModelCalibration(model, model_id, X_test, t, X0_test, 0)
			mc.levenberg_marquardt()
			model.print(precision = precision)
			print("\n")

			simulation = model.simulate(X0_test, t = t_test)
		elif calibration_mode == "Bayes":
			mc = ModelCalibration(model, model_id, X_test, t, X0_test, 0)
			mc.bayesian_calibration()
			mc.traceplot()
			mc.plot_posterior()
			mc.plot_pair()
			X0_test = mc.summary()
			print("\n")
			model.print(precision = precision)
			print("\n")

			simulation, simulation_min, simulation_max = mc.get_simulation()
	except:
		print("Modelo " + str(model_id+1) + " não pode ser simulado ou recalibrado" + "\n")
		continue

	# Generate figures
	orange = "#F46E2B"
	blue = "#06608F"
	green = "#00B849"

	plt.rcParams.update({'font.size': 20})
	fig, ax = plt.subplots(1, 1, dpi = 300)
	ax.plot(t_plot, X_plot[:,0], "ko", label = r"Data $\mathcal{N}$", alpha = 0.5, markersize = 3)
	ax.plot(t_plot, X_plot[:,1], "k^", label = r"Data $\mathcal{Q}$", alpha = 0.5, markersize = 3)
	ax.plot(t_plot, X_plot[:,2], "ks", label = r"Data $\mathcal{P}$", alpha = 0.5, markersize = 3)
	ax.plot(t_test, simulation[:,0], color = orange, label = r"Model $\mathcal{N}$", alpha = 1.0, linewidth = 1)
	ax.plot(t_test, simulation[:,1], color = blue, label = r"Model $\mathcal{Q}$", alpha = 1.0, linewidth = 1)
	ax.plot(t_test, simulation[:,2], color = green, label = r"Model $\mathcal{P}$", alpha = 1.0, linewidth = 1)
	if calibration_mode == "Bayes":
		ax.fill_between(t, simulation_min[:,0], simulation_max[:,0], color = orange, alpha = 0.4)
		ax.fill_between(t, simulation_min[:,1], simulation_max[:,1], color = blue, alpha = 0.4)
		ax.fill_between(t, simulation_min[:,2], simulation_max[:,2], color = green, alpha = 0.4)
	ax.set(xlabel = r"Time (h)", ylabel = r"Number of cells",
		xlim = [t_plot[0], t_plot[-1]+1], ylim = [0, 501],
		xticks = np.arange(t_plot[0], t_plot[-1]+1, step=100), yticks = np.arange(0, 501, step=100)
		# title = "N' = " + model.equations(precision = precision)[0] + "\n"
		# + "Q' = " + model.equations(precision = precision)[1] + "\n"
		# + "P' = " + model.equations(precision = precision)[2] + "\n"
		# + "Initial condition = " + str(X0_test)
	)
	ax.legend(prop={'size': 14})
	# fig.suptitle(optimizer_method + " - Model " + str(model_id+1), fontsize = 16, y = 0.99)
	# fig.show()
	plt.savefig(os.path.join("output", "model" + str(model_id+1) + "_ic0.png"), bbox_inches = 'tight')
	plt.close()

	# Compute SSE
	sse = ms.compute_SSE(X_test.reshape(simulation.shape), simulation)

	# Set mean SSE to the model
	ms.set_model_SSE(model_id, sse)

# Compute AIC and AICc
best_AIC_model = ms.compute_AIC()
best_AICc_model = ms.compute_AICc()
best_BIC_model = ms.compute_BIC()

# Get best model
print("Melhor modelo AIC = " + str(best_AIC_model+1) + "\n")
print("Melhor modelo AICc = " + str(best_AICc_model+1) + "\n")
print("Melhor modelo BIC = " + str(best_BIC_model+1) + "\n")

# Write results
ms.write_output()
ms.write_AICc_weights()
ms.write_pareto_curve(optimizer_method)
