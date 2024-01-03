# Reduced section plate test
A reduced section plate test is provided. 

Geometry, loading conditions and material parameters can be adoptet in *createRCPTest.py*.

The file *abaqus_v6.env* needs to be placed in the work directory in order to enable the fortran compiler option for cray pointers.

Create the input file from the provided *python* script with:
> abaqus cae noGui=createRCPTest.py

Start the Abaqus job with:
> abaqus job=tensileTest user=../../src/UMAT_DamThermMech_1_H.f -interactive
