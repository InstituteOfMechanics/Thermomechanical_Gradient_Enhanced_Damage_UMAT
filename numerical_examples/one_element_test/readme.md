# One-element test
An One-element test is provided. 

Geometry, loading conditions and material parameters can be adoptet in *createElementTest.py*.

The file *abaqus_v6.env* needs to be placed in the work directory in order to enable the fortran compiler option for cray pointers.

Create the input file from the provided *python* script with:
> abaqus cae noGui=createElementTest.py

Start the Abaqus job with:
> abaqus job=elementTest user=../../src/UMAT_DamThermMech_1_H.f -interactive
