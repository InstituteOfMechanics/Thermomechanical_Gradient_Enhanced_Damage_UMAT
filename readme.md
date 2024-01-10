# Thermomechanical Gradient Enhanced Damage UMAT
This repository contains a user material implementation for Abaqus (`UMAT` and `UMATHT`) for a thermo-mechanically coupled gradient-enhanced damage model. The framework can be adopted to several multi-field problems. 
A major benefit of the user material implementation (e.g. compared to a user element implementation) is the applicability of several integrated Abaqus features, such as contact algorithms, element formulations and solver structures. 

The implementation is provided by Lennart Sobisch (<lennart.sobisch@tu-dortmund.de>) and is documented in the publication Sobisch et al. (2024).
If you publish results based on this framework, please cite the relevant paper.

The repository contains the following folders:
- `numerical_examples`: A folder containing the source code for three numerical examples (one-element test, reduced section plate test, indentation test)
- `src`: A folder contatining the source code of the user material implementation, including a detailed user material.

The implementation framework is a direct continuation of Ostwald et al. (2019), where a user material implementation of two-field problems was proposed.

---
## References
- L. Sobisch, T. Kaiser, T. Furlan, A. Menzel, A user material approach for the solution of multi-field problems in Abaqus: Theoretical foundations, gradient-enhanced damage mechanics and thermo-mechanical coupling, Finite Elements in Analysis & Design ?? (2024) ??-??. doi:??.
- R. Ostwald, E. Kuhl, A. Menzel, On the implementation of finite deformation gradient-enhanced damage models, Computational Mechanics 64 (847-877) (2019). <https://doi.org/10.1007/s00466-019-01684-5>.
