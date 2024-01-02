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
E  	    = 210000. 	# Youngs modulus 		       [Mpa] = [N/mm^2]
nu	    = 0.3 		# Poissons ratio 		       [-]
dens 	= 7.8e-3  	# Density 			           [t/mm^3]
c 	    = 2.0 		# Specific heat capacity       [mJ/tK]
alpha   = 1e-5		# Thermal expansion coeff.	   [1/K]
K0 	    = 50.0e-2	# Thermal conductivity		   [mW/mmK]
thetaZ  = 293.0		# Reference Temperature        [K]
betaD 	= 1e3		# Damage penalty parameter	   [1/Mpa]
kappaD  = 5.0e1		# Damage threshold parameter   [Mpa]		
etaD 	= 1.0e-3	    # Damage evolution parameter	   [1/Mpa]
cD      = 1e-1		# Regularisation parameter 	   [mm^2/Mpa]
cP 	    = 1e-2		# Damping parameter		       [-]

gamma 	= 1 		# Gradient damage switch

time 	= 10.0		# Simulation time		       [s]
dtMax 	= 0.1		# Max. Time Step		           [s]
lc 	    = 1.0		# Charact. element length	   [mm]

disp 	= 1.0		# Displacement		           [mm]



# ----------------------------------------------------------------------------
#  Model
# ----------------------------------------------------------------------------
jobName = 'elementTest'
Mdb()
modelName = 'elementTest'
mdb.models.changeKey(fromName='Model-1', toName=modelName)
model = mdb.models[modelName]


# ----------------------------------------------------------------------------
#  Part
# ----------------------------------------------------------------------------
# Primary-Part

sketch = model.ConstrainedSketch(name='primarySketch', sheetSize=10.)
sketch.rectangle(point1=(0.0, 0.0), point2=(lc, lc))
primaryPart = model.Part(name='primaryPart', 
                         dimensionality=THREE_D, 
                         type=DEFORMABLE_BODY)
primaryPart.BaseSolidExtrude(sketch=sketch, depth=lc)
del model.sketches['primarySketch']
# Set with all primary cells
primaryPart.Set(cells=primaryPart.cells, name='allPrimaryCells') 	
# Set with y0 Face
y0Face = primaryPart.faces.findAt(((lc/2,0,lc/2),))			
primaryPart.Set(faces=y0Face, name='y0Face')
# Set with x0 Face
x0Face = primaryPart.faces.findAt(((0,lc/2,lc/2),))
primaryPart.Set(faces=x0Face, name='x0Face')
# Set with z0 Face
z0Face = primaryPart.faces.findAt(((lc/2,lc/2,0),))
primaryPart.Set(faces=z0Face, name='z0Face')
# Set with y1 Face
y1Face = primaryPart.faces.findAt(((lc/2,lc,lc/2),))
primaryPart.Set(faces=y1Face, name='y1Face')



# Secondary-Part

sketch = model.ConstrainedSketch(name='secondarySketch', sheetSize=10.)
sketch.rectangle(point1=(0.0, 0.0), point2=(lc, lc))
secondaryPart = model.Part(name='secondaryPart', 
                           dimensionality=THREE_D, 
                           type=DEFORMABLE_BODY)
secondaryPart.BaseSolidExtrude(sketch=sketch, depth=lc)
del model.sketches['secondarySketch']
# Set with all secondary cells
secondaryPart.Set(cells=secondaryPart.cells, name='allSecondaryCells')

# ----------------------------------------------------------------------------
#  Materials and Sections
# ----------------------------------------------------------------------------
primaryMat = model.Material(name='primaryMaterial')
primaryMat.Density(table=((dens, ), ))
primaryMat.Depvar(n=13)
primaryMat.UserMaterial(type=THERMOMECHANICAL, 
                        mechanicalConstants=(0,E,nu,c,alpha,K0,thetaZ,dens,betaD,etaD,kappaD,gamma,cD,cP), 
                        thermalConstants=(0,E,nu,c,alpha,K0,thetaZ,dens,betaD,etaD,kappaD,gamma,cD,cP))

secondaryMat = model.Material(name='secondaryMaterial')
secondaryMat.Density(table=((dens, ), ))
secondaryMat.Depvar(n=13)
secondaryMat.UserMaterial(type=THERMOMECHANICAL, 
                          mechanicalConstants=(1,E,nu,c,alpha,K0,thetaZ,dens,betaD,etaD,kappaD,gamma,cD,cP), 
                          thermalConstants=(1,E,nu,c,alpha,K0,thetaZ,dens,betaD,etaD,kappaD,gamma,cD,cP))

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

assembly = model.rootAssembly

assembly.Instance(name='primaryInstance', part=primaryPart, dependent=ON)
assembly.Instance(name='secondaryInstance', part=secondaryPart, dependent=ON)
assembly.translate(instanceList=('secondaryInstance', ), vector=(2.0, 0.0, 0.0))

# ----------------------------------------------------------------------------
#  Step
# ----------------------------------------------------------------------------

step = model.CoupledTempDisplacementStep(name='Step-1', 
                                         previous='Initial', 
                                         initialInc=dtMax, 
                                         maxInc=dtMax, 
                                         deltmx=1e5, 
                                         maxNumInc=int(1e4),
                                         nlgeom=ON)

step.setValues(timePeriod=time)
# step.setValues(matrixStorage=SOLVER_DEFAULT, solutionTechnique=SEPARATED) 
model.fieldOutputRequests['F-Output-1'].setValues(numIntervals=int(10*time))
model.fieldOutputRequests['F-Output-1'].setValues(variables=('S','PE','PEEQ','PEMAG', 
                                                             'LE','U','RF','CF','CSTRESS', 
                                                             'CDISP','NT','HFL','RFL','SDV'))

# ----------------------------------------------------------------------------
#  Mesh
# ----------------------------------------------------------------------------

    
elemType1 = mesh.ElemType(elemCode=C3D8T,elemLibrary=STANDARD,secondOrderAccuracy=OFF,distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6T,elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4T,elemLibrary=STANDARD)
    
primaryPart.setElementType(regions=primaryPart.sets['allPrimaryCells'], 
                           elemTypes=(elemType1, elemType2, elemType3))
secondaryPart.setElementType(regions=secondaryPart.sets['allSecondaryCells'], 
                             elemTypes=(elemType1, elemType2, elemType3))

primaryPart.seedPart(size=lc,deviationFactor=0.1,minSizeFactor=0.1)
secondaryPart.seedPart(size=lc,deviationFactor=0.1,minSizeFactor=0.1)

primaryPart.generateMesh()
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
        mpcType=PIN_MPC, 
        userMode=DOF_MODE_MPC, 
        userType=0, 
        csys=None)
    
    
# ----------------------------------------------------------------------------
# Boundary and initial conditions
# ----------------------------------------------------------------------------

model.DisplacementBC(name='yFix', 
                     createStepName='Step-1', 
                     region=assembly.sets['primaryInstance.y0Face'], 
                     u1=UNSET, u2=0.0, u3=UNSET, 
                     ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, 
                     fixed=OFF, 
                     distributionType=UNIFORM, 
                     fieldName='', 
                     localCsys=None)
    
model.DisplacementBC(name='xFix', 
                     createStepName='Step-1', 
                     region=assembly.sets['primaryInstance.x0Face'], 
                     u1=0.0, u2=UNSET, u3=UNSET, 
                     ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, 
                     fixed=OFF, 
                     distributionType=UNIFORM, 
                     fieldName='', 
                     localCsys=None)
    
    
model.DisplacementBC(name='zFix', 
                     createStepName='Step-1', 
                     region=assembly.sets['primaryInstance.z0Face'], 
                     u1=UNSET, u2=UNSET, u3=0.0, 
                     ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                     amplitude=UNSET, 
                     fixed=OFF, 
                     distributionType=UNIFORM, 
                     fieldName='', 
                     localCsys=None)
    


model.DisplacementBC(name='Load', 
                      createStepName='Step-1', 
                      region=assembly.sets['primaryInstance.y1Face'], 
                      u1=UNSET, u2=disp, u3=UNSET, 
                      ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                      amplitude=UNSET, 
                      fixed=OFF, 
                      distributionType=UNIFORM, 
                      fieldName='', 
                      localCsys=None)


model.Temperature(name='Theta0', 
                  createStepName='Initial', 
                  region=assembly.sets['primaryInstance.allPrimaryCells'], 
                  distributionType=UNIFORM, 
                  crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, 
                  magnitudes=(293.0,))

model.Temperature(name='phi0', 
                  createStepName='Initial', 
                  region=assembly.sets['secondaryInstance.allSecondaryCells'], 
                  distributionType=UNIFORM, 
                  crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, 
                  magnitudes=(0.0,))
# ----------------------------------------------------------------------------
# Initialization of damage variable
# ----------------------------------------------------------------------------

numEl = len(secondaryPart.elements)
numElQp = 8*numEl

primaryMat.userMaterial.setValues(mechanicalConstants=(0.0,E,nu,c,alpha,K0,thetaZ,dens, 
                                                       betaD,etaD,kappaD,gamma,cD,cP,numEl,numElQp), 
                                  thermalConstants=(0.0,E,nu,c,alpha,K0,thetaZ,dens, 
                                                    betaD,etaD,kappaD,gamma,cD,cP,numEl,numElQp))

secondaryMat.userMaterial.setValues(mechanicalConstants=(1.0,E,nu,c,alpha,K0,thetaZ,dens, 
                                                         betaD,etaD,kappaD,gamma,cD,cP,numEl,numElQp), 
                                    thermalConstants=(1.0,E,nu,c,alpha,K0,thetaZ,dens, 
                                                      betaD,etaD,kappaD,gamma,cD,cP,numEl,numElQp))

# ----------------------------------------------------------------------------
# Job
# ----------------------------------------------------------------------------

mdb.Job(name=jobName, 
        model='elementTest', 
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
        nodalOutputPrecision=SINGLE,    
        echoPrint=OFF, 
        modelPrint=OFF, 
        contactPrint=OFF, 
        historyPrint=OFF, 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numGPUs=0)

mdb.jobs[jobName].writeInput(consistencyChecking=OFF)
