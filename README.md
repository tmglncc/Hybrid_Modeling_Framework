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

2. Enter the directory of the experiment you would like to run:
   - _chem-400-1-B_: experiment with application of cisplatin after _t = 400_ h;
   - _chem-400-1-C_: experiment with application of taxotere after _t = 400_ h; or
   - _no-chem_: experiment with the application of any drug after _t = 400_ h.
  
3. Run the Python script. For example:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`$ python3 direct_multiple_shooting.py`

4. Check out the results.

## Repository organization

## Cite as
