from casadi import *
import numpy as np
import matplotlib.pyplot as plt
import csv

chem_protocol = "C"

t0 = 100.0
tf = 400.0
Nt = 300

time_ode_sol = np.linspace(t0, tf, Nt+1)

n1 = -6.393e-5
n2 = 9.130e-5
q1 = -7.012e-3
q2 = 3.615e-2
q3 = -2.634e-5
p1 = 1.210e-2
p2 = -2.229e-2
eta = 0.05634

x0 = [0.161, 16.431, 11.271]

N = MX.sym('N')
Q = MX.sym('Q')
P = MX.sym('P')
x = vertcat(N, Q, P)

xdot = vertcat(n1*x[0]*x[1] + n2*x[1]**2.0,
	q1*x[1] + q2*x[2] + q3*x[0]*x[1],
	p1*x[1] + p2*x[2]
)

ode = {}
ode['x']   = x
ode['ode'] = xdot

x_plot = np.array([x0])
for i in range(Nt):
	F = integrator('F', 'cvodes', ode, time_ode_sol[i], time_ode_sol[i+1])
	res = F(x0=x_plot[-1])
	x_plot = np.vstack((x_plot, res["xf"].full().flatten()))

# plt.figure()
# plt.plot(time_ode_sol, x_plot[:,0], 'b-', label=r'$N(t)$')
# plt.plot(time_ode_sol, x_plot[:,1], 'g-', label=r'$Q(t)$')
# plt.plot(time_ode_sol, x_plot[:,2], 'r-', label=r'$P(t)$')
# plt.legend(loc='best')
# plt.xlabel(r'$t$ (h)')
# plt.ylabel('State variables')

# plt.show()

t0 = tf
tf = t0 + 300.0
Nt = 300

time_ocp = np.linspace(t0, tf, Nt+1)

w_u_list = [1.0] # [1.0, 100.0]
w_T_list = [1.0] # [1.0, 100.0]
D_ub_list = [2] # [2, 16]
use_tumor_at_final_time_list = [True] # [False, True]
use_D_list = [True]

with open("experiments.csv", mode='w') as experiments_file:
	experiments_writer = csv.writer(experiments_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

	experiments_writer.writerow(["w_u_val", "w_T_val", "D_ub_val", "use_tumor_at_final_time", "use_D", "lambda_val", "fig_number", "final_obj"])

	exp_number = 1
	fig_number = 1
	for w_u_val in w_u_list:
		for w_T_val in w_T_list:
			for D_ub_val in D_ub_list:
				for use_tumor_at_final_time in use_tumor_at_final_time_list:
					for use_D in use_D_list:
						if use_D:
							lambda_list = [0.20733444] # [0.24649569]
						else:
							lambda_list = [0.0]

						for lambda_val in lambda_list:
							print("Running experiment " + str(exp_number))
							exp_number += 1

							u_max = 1.0
							w_u = w_u_val
							w_T = w_T_val
							lambda_ = lambda_val
							D_ub = D_ub_val

							Nc = MX.sym('Nc')
							Qc = MX.sym('Qc')
							Pc = MX.sym('Pc')
							D = MX.sym('D')
							Ju = MX.sym('Ju')
							x = vertcat(Nc, Qc, Pc, D, Ju)
							u = MX.sym('u')

							if use_D:
								xdot = vertcat(n1*x[0]*x[1] + n2*x[1]**2.0,
									q1*x[1] + q2*x[2] + q3*x[0]*x[1],
									p1*x[1] + p2*x[2] - lambda_*x[3]*x[2],
									u*u_max - eta*x[3],
									w_u*(u*u_max)**2.0
								)
							else:
								xdot = vertcat(n1*x[0]*x[1] + n2*x[1]**2.0,
									q1*x[1] + q2*x[2] + q3*x[0]*x[1],
									p1*x[1] + p2*x[2] - u*u_max*x[2],
									u*u_max,
									w_u*(u*u_max)**2.0
								)

							if use_tumor_at_final_time:
								L = w_u*(u*u_max)**2.0
							else:
								L = w_u*(u*u_max)**2.0 + w_T*(x[1] + x[2])**2.0

							# ode = {'x': x, 'p': u, 'ode': xdot, 'quad': L}
							# F = integrator('F', 'cvodes', ode, t0, (tf - t0)/Nt)

							M = 4
							DT = (tf - t0)/Nt/M

							f = Function('f', [x, u], [xdot, L])

							X0 = MX.sym('X0', 5)
							D = MX.sym('D')

							X = X0
							Q = 0
							for j in range(M):
							   k1, k1_q = f(X, D)
							   k2, k2_q = f(X + DT/2 * k1, D)
							   k3, k3_q = f(X + DT/2 * k2, D)
							   k4, k4_q = f(X + DT * k3, D)
							   X = X + DT/6*(k1 + 2*k2 + 2*k3 + k4)
							   Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q)
							F = Function('F', [X0, D], [X, Q], ['x0', 'p'], ['xf', 'qf'])

							Fk = F(x0=[x_plot[-1,0], x_plot[-1,1], x_plot[-1,2], 0, 0], p=0)
							print(Fk['xf'])
							print(Fk['qf'])

							u_start = [DM(0)]*Nt

							x0 = [x_plot[-1,0], x_plot[-1,1], x_plot[-1,2], 0, 0]
							xk = DM(x0)
							x_start = [xk]
							for k in range(Nt):
								xk = F(x0=xk, p=u_start[k])['xf']
								x_start += [xk]

							w = []
							w0 = []
							lbw = []
							ubw = []
							J = 0
							g = []
							lbg = []
							ubg = []

							Xk = MX.sym('X0', 5)
							w += [Xk]
							lbw += x0
							ubw += x0
							w0 += [x_start[0]]

							for k in range(Nt):
								Dk = MX.sym('D_' + str(k))
								w   += [Dk]
								lbw += [0]
								ubw += [1]
								w0  += [u_start[k]]

								Fk = F(x0=Xk, p=Dk)
								Xk_end = Fk['xf']
								J = J + Fk['qf']

								Xk = MX.sym('X_' + str(k+1), 5)
								w   += [Xk]
								lbw += [0, 0, 0, 0, 0]
								ubw += [inf,  inf, inf, D_ub, inf]
								w0  += [x_start[k+1]]

								g   += [Xk_end - Xk]
								lbg += [0, 0, 0, 0, 0]
								ubg += [0, 0, 0, 0, 0]
							if use_tumor_at_final_time:
								J = J + w_T*(Xk_end[1] + Xk_end[2])**2.0

							try:
								prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)}
								solver = nlpsol('solver', 'ipopt', prob)

								sol = solver(x0=vertcat(*w0), lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
								J_opt = sol['f'].full().flatten()
								w_opt = sol['x'].full().flatten()

								Nc_opt = w_opt[0::6]
								Qc_opt = w_opt[1::6]
								Pc_opt = w_opt[2::6]
								D_opt = w_opt[3::6]
								Ju_opt = w_opt[4::6]
								u_opt = w_opt[5::6]
							except:
								experiments_writer.writerow([w_u_val, w_T_val, D_ub_val, use_tumor_at_final_time, use_D, lambda_val, "Solution not found", "Solution not found"])
								continue

							time_solution = np.concatenate([time_ode_sol[:-1], time_ocp])
							u_solution = np.concatenate([np.zeros(time_ode_sol.shape[0])[:-1], u_opt])
							N_solution = np.concatenate([x_plot[:,0][:-1], Nc_opt])
							Q_solution = np.concatenate([x_plot[:,1][:-1], Qc_opt])
							P_solution = np.concatenate([x_plot[:,2][:-1], Pc_opt])

							if use_D:
								D_solution = np.concatenate([np.zeros(time_ode_sol.shape[0])[:-1], D_opt])

								plt.rcParams.update({'font.size': 20})
								fig, ax = plt.subplots(1, 1, dpi = 300)
								plt.plot(time_solution, D_solution, '-', color = "tab:purple")
								ax.set(xlim = [0, time_solution[-1]+1], ylim = [-0.1, 2.1],
									xticks = np.arange(0, time_solution[-1]+1, step=100), yticks = np.arange(0, 2.1, step=0.25)
								)
								plt.xlabel(r'Time (h)')
								plt.ylabel(r'Drug bioavailability')
								plt.savefig("Figure_" + str(fig_number) + "_control_bioavailability.pdf", bbox_inches = 'tight')
								plt.close()

								with open('data_control_bioavailability.txt', 'w') as f:
									f.write(f'time = {time_solution[np.where(D_solution != 0.0)].astype(int).tolist()}; \n')
									f.write(f'D = {D_solution[np.where(D_solution != 0.0)].tolist()}; \n')
									f.write(f'protocol = {[chem_protocol]*len(np.where(D_solution != 0.0)[0])}; \n')

							plt.rcParams.update({'font.size': 20})
							fig, ax = plt.subplots(1, 1, dpi = 300)
							plt.plot(time_solution, np.append(np.nan, u_solution), '-', color = "tab:pink")
							ax.set(xlim = [0, time_solution[-1]+1], ylim = [-0.05, 1.05],
								xticks = np.arange(0, time_solution[-1]+1, step=100), yticks = np.arange(0, 1.1, step=0.2)
							)
							plt.xlabel(r'Time (h)')
							plt.ylabel(r'Drug dose')
							plt.savefig("Figure_" + str(fig_number) + "_control.pdf", bbox_inches = 'tight')
							plt.close()

							with open('data_control.txt', 'w') as f:
								f.write(f'time = {time_solution[np.where(u_solution != 0.0)].astype(int).tolist()}; \n')
								f.write(f'u = {u_solution[np.where(u_solution != 0.0)].tolist()}; \n')
								f.write(f'protocol = {[chem_protocol]*len(np.where(u_solution != 0.0)[0])}; \n')

							orange = "#F46E2B"
							blue = "#06608F"
							green = "#00B849"

							plt.rcParams.update({'font.size': 20})
							fig, ax = plt.subplots(1, 1, dpi = 300)
							plt.plot(time_solution, N_solution, '-', color = orange, label=r'$\mathcal{N}$')
							plt.plot(time_solution, Q_solution, '-', color = blue, label=r'$\mathcal{Q}$')
							plt.plot(time_solution, P_solution, '-', color = green, label=r'$\mathcal{P}$')
							ax.set(xlim = [0, time_solution[-1]+1], ylim = [0, 501],
								xticks = np.arange(0, time_solution[-1]+1, step=100), yticks = np.arange(0, 501, step=100)
							)
							ax.legend(prop={'size': 14}, loc='upper left')
							plt.xlabel(r'Time (h)')
							plt.ylabel('Number of cells')
							plt.savefig("Figure_" + str(fig_number) + "_tumor.pdf", bbox_inches = 'tight')
							plt.close()

							experiments_writer.writerow([w_u_val, w_T_val, D_ub_val, use_tumor_at_final_time, use_D, lambda_val, fig_number, J_opt[-1]])
							fig_number += 1