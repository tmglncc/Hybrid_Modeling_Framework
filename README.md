# Hybrid Modeling Framework

Despite advances in oncology, cancer persists as a major global health burden, with treatment efficacy often limited by incomplete understanding of its multiscale dynamics spanning tissues, cells, and molecular processes. To address this challenge, we develop a high-fidelity hybrid multiscale model that combines continuous descriptions at the tissue and molecular levels with a discrete individual-based approach at the cellular scale. To optimize chemotherapy treatment within this multiscale framework, we first derive a low-fidelity surrogate model using a data-driven approach integrating Sparse Identification of Nonlinear Dynamics (SINDy) and global sensitivity analysis (SA). This reduced model enables efficient application of optimal control techniques for treatment optimization. Our framework provides a promising strategy for designing effective therapeutic protocols while minimizing drug toxicity in multiscale cancer scenarios.

![Our computational pipeline analyzes the hybrid multiscale model through an optimal control point of view. Using tumor cell dynamics as input data, the SINDy-SA framework derives a reduced-order nonlinear dynamical system. This surrogate model enables the formulation and solution of an optimal control problem, yielding treatment strategies that can be transferred to the original high-fidelity system.](https://drive.google.com/uc?export=view&id=1HJhMmMCZGMdcFLdtgasNlACGFACpEUsy)

## Requirements and running

## Repository organization

## Cite as
