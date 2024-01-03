from abaqus import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from mesh import *
from interaction import *
from load import *
import mesh
from job import *
import numpy as np


# User definition
E           = 210000.            # Youngs modulus               [Mpa] = [N/mm^2]
nu          = 0.3                # Poissons ration              [-]
dens        = 7.8e-3             # Density                      [t/mm^3]
capacity    = 2.0e0              # Specific heat capacity       [mJ/tK]
alpha       = 1e-5               # Thermal expansion coeff.     [1/K]
K0          = 50.0e-2            # Thermal conductivity         [mW/mmK]
thetaZ      = 293.               # Reference temperature        [K]
betaD       = 1e3                # Damage penalty param.        [1/Mpa]
kappaD      = 5.0e1              # Damage threshold param.      [Mpa]
etaD        = 5.0e-3             # Damage evolution param.      [1/Mpa]
cD          = 0.75               # Regularisation param.        [mm^2/Mpa]
cp          = 5e-1               # Damping param.               [-]
gamma       = 1.0                # Gradient damage switch       [-]

time        = 10                 # Simulation time              [s]
dtMax       = 1e0                # Max time increment           [s]

totElems    = 2**10              # Number of elements           [-]

xLG         = 80.0               # Total length of specimen     [mm]
xL1         = 20.0               # Length of inner part         [mm] 
yLG         = 20.0               # Total width of specimen      [mm]
yL1         = 3.75               # Width of inner part          [mm]
zL          = 5.0                # Thickness of specimen        [mm]

disp        = 10.                # Displacement load            [mm]

pos1        = 5.                 # Position of Export point A   [mm]
pos2        = 10.                # Position of Export point B   [mm]

# ----------------------------------------------------------------------------
#  Model
# ----------------------------------------------------------------------------
jobName     = 'tensileTest'
Mdb()
modelName   = 'tensileTest'
mdb.models.changeKey(fromName='Model-1', toName=modelName)
model       = mdb.models[modelName]

# ----------------------------------------------------------------------------
#  Parts
# ----------------------------------------------------------------------------
# Primary-Part
sketch = model.ConstrainedSketch(name='primarySketch', sheetSize=10.)
sketch.Line(point1=(0.0, 0.0), point2=(xLG/2., 0.0))
sketch.Line(point1=(xLG/2.0, 0.0), point2=(xLG/2., yLG/2.))
sketch.Line(point1=(xLG/2., yLG/2.), point2=(xLG/2. - xL1, yLG/2.))
sketch.Line(point1=(0.0, yLG/2. - yL1), point2=(0.0, 0.0))
sketch.Arc3Points(point1=(-xLG/2. + xL1, yLG/2.), 
                  point2=(xLG/2. - xL1, yLG/2.), 
                  point3=(0.0, yLG/2. - yL1))
sketch.autoTrimCurve(curve1=sketch.geometry[6], point1=(-xLG/2. + xL1, yLG/2.))

primaryPart = model.Part(name='primaryPart', dimensionality=THREE_D, type=DEFORMABLE_BODY)
primaryPart.BaseSolidExtrude(sketch=sketch, depth=zL/2.)

del model.sketches['primarySketch']

# Set with all primary cells
primaryPart.Set(cells=primaryPart.cells, name='allPrimaryCells')
# Set with y0 Face
x0Face = primaryPart.faces.findAt(((0.0,yLG/4,zL/4),))
primaryPart.Set(faces=x0Face, name='x0Face')
# Set with y0 Face
y0Face = primaryPart.faces.findAt(((xLG/4,0.0,zL/4),))
primaryPart.Set(faces=y0Face, name='y0Face')
# # Set with x1 Face
x1Face = primaryPart.faces.findAt(((xLG/2.,yLG/4,zL/4),))
primaryPart.Set(faces=x1Face, name='x1Face')
z0Face = primaryPart.faces.findAt(((xLG/4,yLG/4,0.0),))
primaryPart.Set(faces=z0Face, name='z0Face')

# Secondary-Part
sketch = model.ConstrainedSketch(name='secondarySketch', sheetSize=10.)
sketch.Line(point1=(0.0, 0.0), point2=(xLG/2., 0.0))
sketch.Line(point1=(xLG/2.0, 0.0), point2=(xLG/2., yLG/2.))
sketch.Line(point1=(xLG/2., yLG/2.), point2=(xLG/2. - xL1, yLG/2.))
sketch.Line(point1=(0.0, yLG/2. - yL1), point2=(0.0, 0.0))
sketch.Arc3Points(point1=(-xLG/2. + xL1, yLG/2.), 
                  point2=(xLG/2. - xL1, yLG/2.), 
                  point3=(0.0, yLG/2. - yL1))
sketch.autoTrimCurve(curve1=sketch.geometry[6], point1=(-xLG/2. + xL1, yLG/2.)) 
    
secondaryPart = model.Part(name='secondaryPart', dimensionality=THREE_D, type=DEFORMABLE_BODY)
secondaryPart.BaseSolidExtrude(sketch=sketch, depth=zL/2.)

del model.sketches['secondarySketch']

secondaryPart.Set(cells=secondaryPart.cells, name='allSecondaryCells')

# ----------------------------------------------------------------------------
#   Materials and Sections
# ----------------------------------------------------------------------------

primaryMat = model.Material(name='primaryMaterial')
primaryMat.Density(table=((dens, ), ))
primaryMat.Depvar(n=13)
primaryMat.UserMaterial(type=THERMOMECHANICAL, unsymm=ON,
                        mechanicalConstants=(0,E,nu,capacity,alpha,K0,thetaZ,dens,betaD,etaD,kappaD,gamma,cD,cp),
                        thermalConstants=(0,E,nu,capacity,alpha,K0,thetaZ,dens,betaD,etaD,kappaD,gamma,cD,cp))

secondaryMat = model.Material(name='secondaryMaterial')
secondaryMat.Density(table=((dens, ), ))
secondaryMat.Depvar(n=13)
secondaryMat.UserMaterial(type=THERMOMECHANICAL, unsymm=ON,
                          mechanicalConstants=(1,E,nu,capacity,alpha,K0,thetaZ,dens,betaD,etaD,kappaD,gamma,cD,cp), 
                          thermalConstants=(1,E,nu,capacity,alpha,K0,thetaZ,dens,betaD,etaD,kappaD,gamma,cD,cp))

model.HomogeneousSolidSection(name='primarySection', material='primaryMaterial', thickness=None)
model.HomogeneousSolidSection(name='secondarySection', material='secondaryMaterial', thickness=None)

primaryPart.SectionAssignment(region=primaryPart.sets['allPrimaryCells'], 
                              sectionName='primarySection', 
                              offset=0.0, 
                              offsetType=MIDDLE_SURFACE, 
                              offsetField='', 
                              thicknessAssignment=FROM_SECTION)

secondaryPart.SectionAssignment(region=secondaryPart.sets['allSecondaryCells'], 
                                sectionName='secondarySection', 
                                offset=0.0, 
                                offsetType=MIDDLE_SURFACE, 
                                offsetField='', 
                                thicknessAssignment=FROM_SECTION)

# ----------------------------------------------------------------------------
#  Assembly
# ----------------------------------------------------------------------------

assembly            = model.rootAssembly
primaryInstance     = assembly.Instance(name='primaryInstance', part=primaryPart, dependent=ON)
secondaryInstance   = assembly.Instance(name='secondaryInstance', part=secondaryPart, dependent=ON)
assembly.translate(instanceList=('secondaryInstance', ), vector=(xLG*2.0, 0.0, 0.0))

# ----------------------------------------------------------------------------
#  Step
# ----------------------------------------------------------------------------

step = model.CoupledTempDisplacementStep(name='Step-1', previous='Initial', 
                                         initialInc=dtMax, maxInc=dtMax, 
                                         deltmx=100000.0, maxNumInc=1000000,
                                         minInc=1e-15, nlgeom=ON)
step.setValues(timePeriod=time)
model.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'PE', 'PEEQ', 
                                                             'PEMAG', 'LE', 'U', 
                                                             'RF', 'NT', 'HFL', 
                                                             'RFL', 'SDV'))
model.fieldOutputRequests['F-Output-1'].setValues(numIntervals=time*10)

# # ----------------------------------------------------------------------------
# #  Mesh
# # ----------------------------------------------------------------------------
elemType1 = mesh.ElemType(elemCode=C3D8T, elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6T, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4T, elemLibrary=STANDARD)
    
primaryPart.setElementType(regions=primaryPart.sets['allPrimaryCells'], 
                           elemTypes=(elemType1,elemType2,elemType3))

primaryPart.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=xLG/2-xL1)
primaryPart.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=pos1)
primaryPart.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=pos2)

d = primaryPart.datums

primaryPart.PartitionCellByDatumPlane(datumPlane=d[8], cells=primaryPart.cells)
primaryPart.PartitionCellByDatumPlane(datumPlane=d[9], cells=primaryPart.cells)
primaryPart.PartitionCellByDatumPlane(datumPlane=d[10], cells=primaryPart.cells)

zNum = 4

workNum = int(np.sqrt(totElems))
lenA = int(workNum/zNum)
lenB = int(lenA/2)
lenC = int(lenB/2)
lenD = int(lenC/2)

h1 = primaryPart.edges.findAt(((0.,yLG/4,0.0),))
h3 = primaryPart.edges.findAt(((0.,yLG/4,zL/2.),))
h5 = primaryPart.edges.findAt(((xL1,yLG/4,0.0),))
h7 = primaryPart.edges.findAt(((xL1,yLG/4,zL/2.),))
h9 = primaryPart.edges.findAt(((pos1,yLG/4,0.0),))
h11 = primaryPart.edges.findAt(((pos1,yLG/4,zL/2.),))
h13 = primaryPart.edges.findAt(((pos2,yLG/4,0.0),))
h15 = primaryPart.edges.findAt(((pos2,yLG/4,zL/2.),))
h17 = primaryPart.edges.findAt(((xLG/2,yLG/4,0.0),))
h19 = primaryPart.edges.findAt(((xLG/2,yLG/4,zL/2.),))

horizontalLines = (h1,h3,h5,h7,h9,h11,h13,h15,h17,h19)
primaryPart.Set(edges=horizontalLines, name='horizontalLines')

v1 = primaryPart.edges.findAt(((0.,0.,zL/4),))
v2 = primaryPart.edges.findAt(((0.,yLG/2-yL1,zL/4),))
v3 = primaryPart.edges.findAt(((xLG/2.,0.,zL/4),))
v4 = primaryPart.edges.findAt(((xLG/2.,yLG/2.,zL/4),))
v5 = primaryPart.edges.findAt(((xLG/2.-xL1,0.,zL/4),))
v6 = primaryPart.edges.findAt(((xLG/2.-xL1,yLG/2.,zL/4),))
v7 = primaryPart.edges.findAt(((pos1,6.476881,zL/4),))
v8 = primaryPart.edges.findAt(((pos2,7.163213,zL/4),))
v9 = primaryPart.edges.findAt(((pos1,0.,zL/4),))
v10 = primaryPart.edges.findAt(((pos2,0.,zL/4),))

verticalLines = (v1,v2,v3,v4,v5,v6,v7,v8,v9,v10)
primaryPart.Set(edges=verticalLines, name='verticalLines')

d1 = primaryPart.edges.findAt(((pos1/2,0.,0.),))
d2 = primaryPart.edges.findAt(((pos1/2,0.,zL/2),))
d3 = primaryPart.edges.findAt(((pos1+pos1/2,0.,0.),))
d4 = primaryPart.edges.findAt(((pos1+pos1/2,0,zL/2),))
d5 = primaryPart.edges.findAt(((pos2+pos2/2,0,0.),))
d6 = primaryPart.edges.findAt(((pos2+pos2/2,0,zL/2),))
d7 = primaryPart.edges.findAt(((xLG/2-yL1/2,0,0.),))
d8 = primaryPart.edges.findAt(((xLG/2-yL1/2,0,zL/2),))
d9 = primaryPart.edges.findAt(((xLG/2-yL1/2,yLG/2,0.),))
d10 = primaryPart.edges.findAt(((xLG/2-yL1/2,yLG/2,zL/2),))

depthLines = (d7,d8,d9,d10)
primaryPart.Set(edges=depthLines, name='depthLines')

r1 = primaryPart.edges.findAt(((15.066912,8.345734,0.),))
r2 = primaryPart.edges.findAt(((15.066912,8.345734,zL/2),))

r5 = primaryPart.edges.findAt(((2.502573,6.306749,0.),))
r6 = primaryPart.edges.findAt(((2.502573,6.306749,zL/2),))
r7 = primaryPart.edges.findAt(((7.507847,6.762883,0.),))
r8 = primaryPart.edges.findAt(((7.507847,6.762883,zL/2),))

largeRoundLines = (r1,r2)
primaryPart.Set(edges=largeRoundLines, name='largeRoundLines')

smallRoundLines = (r5,r6,r7,r8)
primaryPart.Set(edges=smallRoundLines, name='smallRoundLines')

for i in horizontalLines:
    primaryPart.seedEdgeByNumber(edges=i, number=lenB, constraint=FIXED)
    
for i in verticalLines:
    primaryPart.seedEdgeByNumber(edges=i, number=zNum, constraint=FIXED)
    
for i in depthLines:
    primaryPart.seedEdgeByNumber(edges=i, number=lenA, constraint=FIXED)

for i in largeRoundLines:
    primaryPart.seedEdgeByNumber(edges=i, number=lenB, constraint=FIXED)

for i in smallRoundLines:
    primaryPart.seedEdgeByNumber(edges=i, number=lenC, constraint=FIXED)
    
primaryPart.generateMesh()


secondaryPart.setElementType(regions=secondaryPart.sets['allSecondaryCells'], 
                             elemTypes=(elemType1,elemType2,elemType3))

secondaryPart.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=xLG/2-xL1)
secondaryPart.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=pos1)
secondaryPart.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=pos2)

d = secondaryPart.datums

secondaryPart.PartitionCellByDatumPlane(datumPlane=d[4], cells=secondaryPart.cells)
secondaryPart.PartitionCellByDatumPlane(datumPlane=d[5], cells=secondaryPart.cells)
secondaryPart.PartitionCellByDatumPlane(datumPlane=d[6], cells=secondaryPart.cells)

h1 = secondaryPart.edges.findAt(((0.,yLG/4,0.0),))
h3 = secondaryPart.edges.findAt(((0.,yLG/4,zL/2.),))
h5 = secondaryPart.edges.findAt(((xL1,yLG/4,0.0),))
h7 = secondaryPart.edges.findAt(((xL1,yLG/4,zL/2.),))
h9 = secondaryPart.edges.findAt(((pos1,yLG/4,0.0),))
h11 = secondaryPart.edges.findAt(((pos1,yLG/4,zL/2.),))
h13 = secondaryPart.edges.findAt(((pos2,yLG/4,0.0),))
h15 = secondaryPart.edges.findAt(((pos2,yLG/4,zL/2.),))
h17 = secondaryPart.edges.findAt(((xLG/2,yLG/4,0.0),))
h19 = secondaryPart.edges.findAt(((xLG/2,yLG/4,zL/2.),))

horizontalLines = (h1,h3,h5,h7,h9,h11,h13,h15,h17,h19)
secondaryPart.Set(edges=horizontalLines, name='horizontalLines')

v1 = secondaryPart.edges.findAt(((0.,0.,zL/4),))
v2 = secondaryPart.edges.findAt(((0.,yLG/2-yL1,zL/4),))
v3 = secondaryPart.edges.findAt(((xLG/2.,0.,zL/4),))
v4 = secondaryPart.edges.findAt(((xLG/2.,yLG/2.,zL/4),))
v5 = secondaryPart.edges.findAt(((xLG/2.-xL1,0.,zL/4),))
v6 = secondaryPart.edges.findAt(((xLG/2.-xL1,yLG/2.,zL/4),))
v7 = secondaryPart.edges.findAt(((pos1,6.476881,zL/4),))
v8 = secondaryPart.edges.findAt(((pos2,7.163213,zL/4),))
v9 = secondaryPart.edges.findAt(((pos1,0.,zL/4),))
v10 = secondaryPart.edges.findAt(((pos2,0.,zL/4),))

verticalLines = (v1,v2,v3,v4,v5,v6,v7,v8,v9,v10)
secondaryPart.Set(edges=verticalLines, name='verticalLines')

d1 = secondaryPart.edges.findAt(((pos1/2,0.,0.),))
d2 = secondaryPart.edges.findAt(((pos1/2,0.,zL/2),))
d3 = secondaryPart.edges.findAt(((pos1+pos1/2,0.,0.),))
d4 = secondaryPart.edges.findAt(((pos1+pos1/2,0,zL/2),))
d5 = secondaryPart.edges.findAt(((pos2+pos2/2,0,0.),))
d6 = secondaryPart.edges.findAt(((pos2+pos2/2,0,zL/2),))
d7 = secondaryPart.edges.findAt(((xLG/2-yL1/2,0,0.),))
d8 = secondaryPart.edges.findAt(((xLG/2-yL1/2,0,zL/2),))
d9 = secondaryPart.edges.findAt(((xLG/2-yL1/2,yLG/2,0.),))
d10 = secondaryPart.edges.findAt(((xLG/2-yL1/2,yLG/2,zL/2),))

depthLines = (d7,d8,d9,d10)
secondaryPart.Set(edges=depthLines, name='depthLines')

r1 = secondaryPart.edges.findAt(((15.066912,8.345734,0.),))
r2 = secondaryPart.edges.findAt(((15.066912,8.345734,zL/2),))

r5 = secondaryPart.edges.findAt(((2.502573,6.306749,0.),))
r6 = secondaryPart.edges.findAt(((2.502573,6.306749,zL/2),))
r7 = secondaryPart.edges.findAt(((7.507847,6.762883,0.),))
r8 = secondaryPart.edges.findAt(((7.507847,6.762883,zL/2),))

largeRoundLines = (r1,r2)
secondaryPart.Set(edges=largeRoundLines, name='largeRoundLines')

smallRoundLines = (r5,r6,r7,r8)
secondaryPart.Set(edges=smallRoundLines, name='smallRoundLines')

for i in horizontalLines:
    secondaryPart.seedEdgeByNumber(edges=i, number=lenB, constraint=FIXED)
    
for i in verticalLines:
    secondaryPart.seedEdgeByNumber(edges=i, number=zNum, constraint=FIXED)
    
for i in depthLines:
    secondaryPart.seedEdgeByNumber(edges=i, number=lenA, constraint=FIXED)

for i in largeRoundLines:
    secondaryPart.seedEdgeByNumber(edges=i, number=lenB, constraint=FIXED)

for i in smallRoundLines:
    secondaryPart.seedEdgeByNumber(edges=i, number=lenC, constraint=FIXED)

secondaryPart.generateMesh()
assembly.regenerate()

# ----------------------------------------------------------------------------
# Linear constraints
# ----------------------------------------------------------------------------

for nodeIter in range(len(primaryPart.nodes)):    
    assembly.SetFromNodeLabels(name='P-'+str(nodeIter), 
                               nodeLabels=(('primaryInstance', (primaryPart.nodes[nodeIter].label,)),))
    
for nodeIter in range(len(secondaryPart.nodes)):
    assembly.SetFromNodeLabels(name='S-'+str(nodeIter), 
                               nodeLabels=(('secondaryInstance', (secondaryPart.nodes[nodeIter].label,)),))
    
for nodeIter in range(len(primaryPart.nodes)):
    model.MultipointConstraint(name='P-S-'+str(nodeIter), controlPoint=assembly.sets['P-'+str(nodeIter)], 
                               surface=assembly.sets['S-'+str(nodeIter)], 
                               mpcType=PIN_MPC, userMode=DOF_MODE_MPC, userType=0, csys=None)
    
 
# ----------------------------------------------------------------------------
# Boundary and initial conditions
# ----------------------------------------------------------------------------

model.DisplacementBC(name='xFix', createStepName='Step-1', 
                     region=assembly.sets['primaryInstance.x0Face'], 
                     u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
                     fieldName='', localCsys=None)

model.DisplacementBC(name='yFix', createStepName='Step-1', 
                     region=assembly.sets['primaryInstance.y0Face'], 
                     u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
                     fieldName='', localCsys=None)   
    
model.DisplacementBC(name='zFix', createStepName='Step-1', 
                     region=assembly.sets['primaryInstance.z0Face'], 
                     u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
                     fieldName='', localCsys=None)
    
model.DisplacementBC(name='Load', createStepName='Step-1', 
                     region=assembly.sets['primaryInstance.x1Face'], 
                     u1=disp, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
                     fieldName='', localCsys=None)

model.Temperature(name='Theta0', createStepName='Initial', 
                  region=assembly.sets['primaryInstance.allPrimaryCells'], 
                  distributionType=UNIFORM, 
                  crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, 
                  magnitudes=(293.0,))

model.Temperature(name='phi0', createStepName='Initial', 
                  region=assembly.sets['secondaryInstance.allSecondaryCells'], 
                  distributionType=UNIFORM, 
                  crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, 
                  magnitudes=(0.0,))

# ----------------------------------------------------------------------------
# Initialization of damage variable
# ----------------------------------------------------------------------------

centPoint = primaryInstance.vertices.findAt(((0.,0.,zL/2.),))
assembly.Set(nodes=(centPoint[0].getNodes(),), name='center')

pos1Point = primaryInstance.vertices.findAt(((pos1,0.,zL/2.),))
assembly.Set(nodes=(pos1Point[0].getNodes(),), name='pos1')

pos2Point = primaryInstance.vertices.findAt(((pos2,0.,zL/2.),))
assembly.Set(nodes=(pos2Point[0].getNodes(),), name='pos2')

centPointSec = secondaryInstance.vertices.findAt(((xLG*2,0.,zL/2),))
assembly.Set(nodes=(centPointSec[0].getNodes(),), name='centerSec')

pos1PointSec = secondaryInstance.vertices.findAt(((xLG*2+pos1,0.,zL/2.),))
assembly.Set(nodes=(pos1PointSec[0].getNodes(),), name='pos1Sec')

pos2PointSec = secondaryInstance.vertices.findAt(((xLG*2+pos2,0.,zL/2.),))
assembly.Set(nodes=(pos2PointSec[0].getNodes(),), name='pos2Sec')

numEl = len(primaryPart.elements)
numElQp = 8*numEl

primaryMat.userMaterial.setValues(mechanicalConstants=(0.0, E, nu, capacity, alpha, K0, thetaZ, dens, 
                                                       betaD, etaD, kappaD, gamma, cD, cp, numEl,numElQp), 
                                  thermalConstants=(0.0, E, nu, capacity, alpha, K0, thetaZ, dens, 
                                                    betaD, etaD, kappaD, gamma, cD, cp, numEl,numElQp))

secondaryMat.userMaterial.setValues(mechanicalConstants=(1.0, E, nu, capacity, alpha, K0, thetaZ, dens, 
                                                         betaD, etaD, kappaD, gamma, cD, cp, numEl,numElQp), 
                                  thermalConstants=(1.0, E, nu, capacity, alpha, K0, thetaZ, dens, 
                                                    betaD, etaD, kappaD, gamma, cD, cp, numEl,numElQp))
                                      

# # ----------------------------------------------------------------------------
# # Job
# # ----------------------------------------------------------------------------

mdb.Job(name=jobName, 
        model='tensileTest', 
        description='', 
        type=ANALYSIS, 
        atTime=None, 
        waitMinutes=0, 
        waitHours=0, 
        queue=None, 
        memory=90, 
        memoryUnits=PERCENTAGE, 
        getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, 
        nodalOutputPrecision=FULL,    
        echoPrint=OFF, 
        modelPrint=OFF, 
        contactPrint=OFF, 
        historyPrint=OFF, 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numGPUs=0)

mdb.jobs[jobName].writeInput(consistencyChecking=OFF)
