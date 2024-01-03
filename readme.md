# Thermomechanical Gradient Enhanced Damage UMAT
This repository contains a user material implementation for Abaqus (UMAT and UMATHT) for a thermo-mechanically coupled gradient-enhanced damage model. The framework can be adopted to several multi-field problems. 
A major benefit of the user material implementation (e.g. compared to a user element implementation) is the applicability of several integrated Abaqus features, such as contact algorithms, element formulations and solver structures. 

The implementation is provided by Lennart Sobisch (<lennart.sobisch@tu-dortmund.de>) and is documented in the publication  [^1].
If you publish results based on any of these models, please cite the relevant paper.

The repository contains the following folders:

    src: A folder containing the source code necessary to build the model
    doc: A folder containing documentation for the models
    ref.bib: Latex bibliography file with reference(s) to the appropriate paper(s)

---
## References
[^1]: L. Sobisch, T. Kaiser, T. Furlan, A. Menzel, A user material approach for the solution of multi-field problems in Abaqus: Theoretical foundations, gradient-enhanced damage mechanics and thermo-mechanical coupling, Finite Elements in Analysis & Design ?? (2024) ??-??. doi:??
