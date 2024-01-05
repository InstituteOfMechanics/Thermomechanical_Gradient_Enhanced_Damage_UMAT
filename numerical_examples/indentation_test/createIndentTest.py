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
etaD        = 5e-3               # Damage evolution param.      [1/Mpa]
cD          = 0.5e0              # Regularisation param.        [mm^2/Mpa]
cp          = 1e-2               # Damping param.               [-]
gamma       = 1.0                # Gradient damage switch       [-]

alpha_conv  = 1e2                # Surf. coeff. of heat cond.   [mW/Kmm^2]
K_sphere    = 50e-2              # Thermal conductivity Sphere  [mW/mmK]
Esphere     = 1e6                # Youngs modulus sphere        [Mpa] = [N/mm^2]

time        = 10                 # Simulation time              [s]
dtMax       = 1e0                # Max. time increment          [s]

xLG         = 30.                # Length of geometry           [mm]
yLG         = 30.                # Width of geometry            [mm]
zL          = 10.                # Height of geometry           [mm]
radius      = 8.                 # Radius of sphere             [mm]

disp        = 2.5                # Displacement load            [mm]

lcQuad      = 2.                 # char. Element length plate   [mm]
lcSphere    = 1.5                # char. Element length sphere  [mm]          


# ----------------------------------------------------------------------------
#  Model
# ----------------------------------------------------------------------------

jobName     = 'indentationTest'
Mdb()
modelName   = 'indentationTest'
mdb.models.changeKey(fromName='Model-1', toName=modelName)
model       = mdb.models[modelName]

# ----------------------------------------------------------------------------
#  Parts
# ----------------------------------------------------------------------------

# Sphere part
sketch = model.ConstrainedSketch(name='sphereSketch', sheetSize=10.)
sketch.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
sketch.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(radius, 0.0))
sketch.Line(point1=(0.0, radius), point2=(0.0, 0.0))
sketch.Line(point1=(0.0, 0.0), point2=(radius, 0.0))
sketch.autoTrimCurve(curve1=sketch.geometry[sketch.geometry.keys()[1]], point1=(0.0, -radius))
sketch.autoTrimCurve(curve1=sketch.geometry[sketch.geometry.keys()[-1]], point1=(-radius, 0.))

spherePart = model.Part(name='spherePart', dimensionality=THREE_D, type=DEFORMABLE_BODY)
spherePart.BaseSolidRevolve(sketch=sketch, angle=90, flipRevolveDirection=OFF)

del model.sketches['sphereSketch']

# Set with all sphere cells
spherePart.Set(cells=spherePart.cells, name='allSphereCells')
# Set with z0 Face
x1Point = spherePart.vertices.findAt(((0.,0.,radius),))
spherePart.Set(vertices=x1Point, name='x1Point')
# Set with x0 Face
x0Face = spherePart.faces.findAt(((0.0,radius/2,radius/2),))
spherePart.Set(faces=x0Face, name='x0Face')
# Set with y0 Face
y0Face = spherePart.faces.findAt(((radius/2,0.0,radius/2),))
spherePart.Set(faces=y0Face, name='y0Face')
# Set with x0 Face
z0Face = spherePart.faces.findAt(((radius/2,radius/2,0.0),))
spherePart.Set(faces=z0Face, name='z0Face')


# Primary-Part
sketch = model.ConstrainedSketch(name='primarySketch', sheetSize=10.)
sketch.Line(point1=(0.0, 0.0), point2=(xLG/2, 0.0))
sketch.Line(point1=(xLG/2, 0.0), point2=(xLG/2, yLG/2))
sketch.Line(point1=(xLG/2, yLG/2), point2=(0.0, yLG/2))
sketch.Line(point1=(0.0, yLG/2), point2=(0.0, 0.0))

primaryPart = model.Part(name='primaryPart', dimensionality=THREE_D, type=DEFORMABLE_BODY)
primaryPart.BaseSolidExtrude(sketch=sketch, depth=zL)

del model.sketches['primarySketch']

# Set with all primary cells
primaryPart.Set(cells=primaryPart.cells, name='allPrimaryCells')
# Set with x0 Face
x0Face = primaryPart.faces.findAt(((0.0,yLG/2,zL/2),))
primaryPart.Set(faces=x0Face, name='x0Face')
# Set with y0 Face
y0Face = primaryPart.faces.findAt(((xLG/2,0.0,zL/2),))
primaryPart.Set(faces=y0Face, name='y0Face')
# Set with z0 Face
z0Face = primaryPart.faces.findAt(((xLG/2,yLG/2,zL),))
primaryPart.Set(faces=z0Face, name='z0Face')


# Secondary-Part
sketch = model.ConstrainedSketch(name='secondarySketch', sheetSize=10.)
sketch.Line(point1=(0.0, 0.0), point2=(xLG/2, 0.0))
sketch.Line(point1=(xLG/2, 0.0), point2=(xLG/2, yLG/2))
sketch.Line(point1=(xLG/2, yLG/2), point2=(0.0, yLG/2))
sketch.Line(point1=(0.0, yLG/2), point2=(0.0, 0.0))
    
secondaryPart = model.Part(name='secondaryPart', dimensionality=THREE_D, type=DEFORMABLE_BODY)
secondaryPart.BaseSolidExtrude(sketch=sketch, depth=zL)

del model.sketches['secondarySketch']

secondaryPart.Set(cells=secondaryPart.cells, name='allSecondaryCells')


# ----------------------------------------------------------------------------
#   Materials and Sections
# ----------------------------------------------------------------------------

sphereMat = model.Material(name='sphereMaterial')
sphereMat.Density(table=((dens, ), ))
sphereMat.Elastic(table=((Esphere, 0.3), ))
sphereMat.Conductivity(table=((K_sphere, ), ))
sphereMat.SpecificHeat(table=((capacity/dens, ), ))

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
                          mechanicalConstants=(1,E,nu,capacity,alpha,K0,thetaZ,1.,betaD,etaD,kappaD,gamma,cD,cp), 
                          thermalConstants=(1,E,nu,capacity,alpha,K0,thetaZ,1.,betaD,etaD,kappaD,gamma,cD,cp))


model.HomogeneousSolidSection(name='sphereSection', material='sphereMaterial', thickness=None)
model.HomogeneousSolidSection(name='primarySection', material='primaryMaterial', thickness=None)
model.HomogeneousSolidSection(name='secondarySection', material='secondaryMaterial', thickness=None)


spherePart.SectionAssignment(region=spherePart.sets['allSphereCells'], 
                              sectionName='sphereSection', offset=0.0, 
                              offsetType=MIDDLE_SURFACE, offsetField='', 
                              thicknessAssignment=FROM_SECTION)

primaryPart.SectionAssignment(region=primaryPart.sets['allPrimaryCells'], 
                              sectionName='primarySection', offset=0.0, 
                              offsetType=MIDDLE_SURFACE, offsetField='', 
                              thicknessAssignment=FROM_SECTION)

secondaryPart.SectionAssignment(region=secondaryPart.sets['allSecondaryCells'], 
                                sectionName='secondarySection', offset=0.0, 
                                offsetType=MIDDLE_SURFACE, offsetField='', 
                                thicknessAssignment=FROM_SECTION)

# ----------------------------------------------------------------------------
#  Assembly
# ----------------------------------------------------------------------------

assembly = model.rootAssembly

sphereInstance      = assembly.Instance(name='sphereInstance', part=spherePart, dependent=ON)
primaryInstance     = assembly.Instance(name='primaryInstance', part=primaryPart, dependent=ON)
secondaryInstance   = assembly.Instance(name='secondaryInstance', part=secondaryPart, dependent=ON)

assembly.translate(instanceList=('secondaryInstance', ), vector=(xLG*2.0, 0.0, 0.0))
assembly.translate(instanceList=('sphereInstance', ), vector=(0.0, 0.0, -radius))

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

step.control.setValues(allowPropagation=OFF, resetDefaultValues=OFF, 
                       timeIncrementation=(4.0, 8.0, 9.0, 50., 40.0, 4.0, 12.0, 5.0, 6.0, 3.0, 50.0))

# # ----------------------------------------------------------------------------
# #  Mesh
# # ----------------------------------------------------------------------------

    
elemType1 = mesh.ElemType(elemCode=C3D8T, elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6T, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4T, elemLibrary=STANDARD)

spherePart.setElementType(regions=spherePart.sets['allSphereCells'], elemTypes=(elemType1,elemType2,elemType3))

spherePart.seedPart(size=lcSphere, deviationFactor=0.1, minSizeFactor=0.1)
spherePart.generateMesh()

primaryPart.setElementType(regions=primaryPart.sets['allPrimaryCells'], elemTypes=(elemType1,elemType2,elemType3))
primaryPart.seedPart(size=lcQuad, deviationFactor=0.1, minSizeFactor=0.1)
primaryPart.generateMesh()

secondaryPart.setElementType(regions=secondaryPart.sets['allSecondaryCells'], elemTypes=(elemType1,elemType2,elemType3))
secondaryPart.seedPart(size=lcQuad, deviationFactor=0.1, minSizeFactor=0.1)
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
    model.MultipointConstraint(name='P-S-'+str(nodeIter), 
                               controlPoint=assembly.sets['P-'+str(nodeIter)], 
                               surface=assembly.sets['S-'+str(nodeIter)], 
                               mpcType=PIN_MPC, userMode=DOF_MODE_MPC, 
                               userType=0, csys=None)
    
contact = model.ContactProperty('IntProp-1')

contact.TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC, 
                           slipRateDependency=OFF, pressureDependency=OFF, 
                           temperatureDependency=OFF, dependencies=0, 
                           table=((0.2, ), ), shearStressLimit=None, 
                           maximumElasticSlip=FRACTION, fraction=0.005, 
                           elasticSlipStiffness=None)

contact.NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON, 
                       constraintEnforcementMethod=DEFAULT)

contact.ThermalConductance(definition=TABULAR, clearanceDependency=ON, 
                           pressureDependency=OFF, temperatureDependencyC=OFF, 
                           massFlowRateDependencyC=OFF, dependenciesC=0, 
                           clearanceDepTable=((alpha_conv, 0.0), (0.0, 0.01)))

generalCont = model.ContactStd(name='Int-1', createStepName='Initial')
generalCont.includedPairs.setValuesInStep(stepName='Initial', useAllstar=ON)
generalCont.contactPropertyAssignments.appendInStep(stepName='Initial', 
                                                    assignments=((GLOBAL, SELF, 'IntProp-1'), ))
 
# ----------------------------------------------------------------------------
# Boundary and initial conditions
# ----------------------------------------------------------------------------

model.DisplacementBC(name='xPFix', createStepName='Initial', 
                     region=assembly.sets['primaryInstance.x0Face'], 
                     u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
                     fieldName='', localCsys=None)

model.DisplacementBC(name='yPFix', createStepName='Initial', 
                     region=assembly.sets['primaryInstance.y0Face'], 
                     u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
                     fieldName='', localCsys=None)

model.DisplacementBC(name='zPFix', createStepName='Initial', 
                     region=assembly.sets['primaryInstance.z0Face'], 
                     u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
                     fieldName='', localCsys=None)

model.DisplacementBC(name='xSFix', createStepName='Initial', 
                     region=assembly.sets['sphereInstance.x0Face'], 
                     u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
                     fieldName='', localCsys=None)

model.DisplacementBC(name='ySFix', createStepName='Initial', 
                     region=assembly.sets['sphereInstance.y0Face'], 
                     u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
                     fieldName='', localCsys=None)

model.DisplacementBC(name='zLoad', createStepName='Step-1', 
                     region=assembly.sets['sphereInstance.z0Face'], 
                     u1=UNSET, u2=UNSET, u3=disp, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
                     fieldName='', localCsys=None)

model.Temperature(name='Theta0_sphere', createStepName='Initial', 
                  region=assembly.sets['sphereInstance.allSphereCells'], 
                  distributionType=UNIFORM, 
                  crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, 
                  magnitudes=(393.0,))

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

numEl = len(primaryPart.elements)
numElQp = 8*numEl

primaryMat.userMaterial.setValues(mechanicalConstants=(0.0, E, nu, capacity, alpha, K0, thetaZ, dens, 
                                                       betaD, etaD, kappaD, gamma, cD, cp, numEl,numElQp), 
                                  thermalConstants=(0.0, E, nu, capacity, alpha, K0, thetaZ, dens, 
                                                    betaD, etaD, kappaD, gamma, cD, cp, numEl,numElQp))

secondaryMat.userMaterial.setValues(mechanicalConstants=(1.0, E, nu, capacity, alpha, K0, thetaZ, 1., 
                                                         betaD, etaD, kappaD, gamma, cD, cp, numEl,numElQp), 
                                  thermalConstants=(1.0, E, nu, capacity, alpha, K0, thetaZ, 1., 
                                                    betaD, etaD, kappaD, gamma, cD, cp, numEl,numElQp))
                                      

# # ----------------------------------------------------------------------------
# # Job
# # ----------------------------------------------------------------------------

mdb.Job(name=jobName, 
        model='indentationTest', 
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
