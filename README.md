# Hybrid Modeling Framework

Despite advances in oncology, cancer persists as a major global health burden, with treatment efficacy often limited by incomplete understanding of its multiscale dynamics spanning tissues, cells, and molecular processes. To address this challenge, we develop a high-fidelity hybrid multiscale model that combines continuous descriptions at the tissue and molecular levels with a discrete individual-based approach at the cellular scale. To optimize chemotherapy treatment within this multiscale framework, we first derive a low-fidelity surrogate model using a data-driven approach integrating Sparse Identification of Nonlinear Dynamics (SINDy) and global sensitivity analysis (SA). This reduced model enables efficient application of optimal control techniques for treatment optimization. Our framework provides a promising strategy for designing effective therapeutic protocols while minimizing drug toxicity in multiscale cancer scenarios.

![Our computational pipeline analyzes the hybrid multiscale model through an optimal control point of view. Using tumor cell dynamics as input data, the SINDy-SA framework derives a reduced-order nonlinear dynamical system. This surrogate model enables the formulation and solution of an optimal control problem, yielding treatment strategies that can be transferred to the original high-fidelity system.](https://drive.google.com/uc?export=view&id=1HJhMmMCZGMdcFLdtgasNlACGFACpEUsy)

## Requirements and running

Our experiments have been performed using a **Hybrid Multiscale Model** (version 1.0), the **SINDy-SA** framework (version 1.0), and an **Optimal Control** Solver. The implementation of the hybrid multiscale model is predominantly in C++, while both SINDy-SA framework and optimal control solver are built in Python.

For requirements and instructions on how to run the hybrid multiscale model, please access [https://github.com/tmglncc/Hybrid_Multiscale_Model](https://github.com/tmglncc/Hybrid_Multiscale_Model).

For requirements and instructions on how to run the SINDy-SA framework, please access [https://github.com/tmglncc/SINDy-SA](https://github.com/tmglncc/SINDy-SA).

The optimal control solver uses **Python 3.8.2** and the following modules are required:

- CasADi 3.6.5 ([https://web.casadi.org/](https://web.casadi.org/))
- NumPy 1.20.3 ([https://numpy.org/](https://numpy.org/))
- Matplotlib 3.7.2 ([https://matplotlib.org/](https://matplotlib.org/))

To run the optimal control solver, please follow the instructions below:

1. Clone this repository directly from terminal:
	 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`$ git clone https://github.com/tmglncc/Hybrid_Modeling_Framework.git`
	
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**OR**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Download the _.zip_ file and decompress it.

2. Enter the _Optimal Control_ folder.

3. Enter the directory of the experiment you would like to run:
   - _chem-400-1-B_: experiment with application of cisplatin after _t = 400_ h;
   - _chem-400-1-C_: experiment with application of taxotere after _t = 400_ h; or
   - _no-chem_: experiment with the application of any drug after _t = 400_ h.
  
4. Run the Python script, for example:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`$ python3 direct_multiple_shooting.py`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Notice that the Python script filename may be different for different experiments.

5. Check out the results.

## Repository organization

This GitHub repository is organized according to the following directories and files:

- **Hybrid Multiscale Model** folder: contains the experiments that implement a hybrid multiscale model for tumor growth with chemotherapy.

	- _chem-400-1-B_ subfolder: experiment with application of cisplatin at t = 400 h;
	- _chem-400-1-C_ subfolder: experiment with application of taxotere at t = 400 h;
	- _no-chem_ subfolder: control experiment.

 - **Low-Fidelity Surrogate Model** folder: contains the Python scripts that implement the SINDy-SA framework to identify a low-fidelity surrogate model from simulated data of the high-fidelity hybrid multiscale model.

 - **Optimal Control** folder: contains the experiments that solve an optimal control problem based on the identified low-fidelity surrogate model.

	- _chem-400-1-B_ subfolder: experiment with application of cisplatin after _t = 400_ h;
	- _chem-400-1-C_ subfolder: experiment with application of taxotere after _t = 400_ h;
	- _no-chem_ subfolder: experiment with the application of any drug after _t = 400_ h.

 - **Hybrid Multiscale Model with Optimal Control** folder: contains the experiments that implement the hybrid multiscale model for tumor growth using the optimal chemotherapy administration protocol.

	- _chem-400-1-B_ subfolder: experiment with application of the optimized cisplatin;
	- _chem-400-1-C_ subfolder: experiment with application of the optimized taxotere.

## Cite as

Naozuka, G.T.; Lima, E.A.B.F.; Almeida, R.C. Hybrid Modeling Framework, 2025. Version 1.0. Available online: [https://github.com/tmglncc/Hybrid_Modeling_Framework](https://github.com/tmglncc/Hybrid_Modeling_Framework) (accessed on 13 May 2025), doi: [10.5281/zenodo.15399549](https://doi.org/10.5281/zenodo.15399549)
