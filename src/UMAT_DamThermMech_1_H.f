!     ABAQUS Subroutine for gradient-enhanced damage model with
!     Thermo-mechanical coupling

        SUBROUTINE UFIELD(FIELD,KFIELD,NSECPT,KSTEP,KINC,TIME,NODE,
     &                    COORDS,TEMP,DTEMP,NFIELD)

            INCLUDE 'ABA_PARAM.INC'

            DIMENSION FIELD(NSECPT,NFIELD), TIME(2), COORDS(3),
     &                TEMP(NSECPT), DTEMP(NSECPT)



            write(*,*) 'UFIELD'


            RETURN
        END SUBROUTINE
      
      
!-----------------------------------------------------------------------
!
!                   U M A T - S U B R O U T I N E
!
!-----------------------------------------------------------------------

        SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     &                  DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,
     &                  DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,
     &                  PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,
     &                  DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
     
            IMPLICIT REAL*8(a-h,o-z)
              
            INCLUDE 'SMAAspUserArrays.hdr'
        
            CHARACTER*80 CMNAME
            DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     &                DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),
     &                DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     &                PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),
     &                DFGRD1(3,3)
            
            
            ! Character definitions
            CHARACTER(LEN=3)                :: tangent
            CHARACTER(LEN=256)              :: OUTDIR, path
              
            ! Internal variables
            DOUBLE PRECISION                :: E, nu, mu, lam, c, alpha
            DOUBLE PRECISION                :: K0, thetaZ, dens, betaD
            DOUBLE PRECISION                :: etaD, dD, gamma, cD, cP, detJ         
            INTEGER                         :: elemForm, totElems
            INTEGER                         :: totIntPoints 
            INTEGER                         :: LENOUTDIR,NOEL
            DOUBLE PRECISION, DIMENSION(13) :: matpars

            ! Get the Work directory for the file transfer
            CALL GETOUTDIR(OUTDIR, LENOUTDIR)
            
            ! Choose tangent type
            !tangent = 'ana'
            tangent = 'num'
     
            !---------------------------------------------------------
            ! Read material paramters
              
            ! The switch for the considered material model is stored
            ! in the material parameters.
            ! PROPS(1) = 0 will trigger the "master-element" with the
            !              temperature field.
            ! PROPS(1) = 1 will trigger the "slave-element" with the
            !              damage field.
              
            ! element Switch
            elemForm = PROPS(1)
            ! youngs modulus
            E        = PROPS(2)
            ! poissions ratio
            nu       = PROPS(3)
            ! specific heat capacity
            c        = PROPS(4)
            ! thermal expansion coefficient
            alpha    = PROPS(5)
            ! thermal conductivity
            K0       = PROPS(6)
            ! reference temperature
            thetaZ   = PROPS(7)
            ! density
            dens     = PROPS(8)
            ! penalty parameter
            betaD    = PROPS(9)
            ! damage saturation parameter
            etaD     = PROPS(10)
            ! damage initiation parameter
            dD       = PROPS(11)
            ! regularisation switch
            gamma    = PROPS(12)
            ! regularisation parameter
            cD       = PROPS(13)
            ! numerical capacity for damage
            cP       = PROPS(14)
            
            ! convert 'E' and 'nu' to Lame parameters (mue and lambda)
            mu       = E/(2.0d0*(1.0d0+nu))
            lam      = E*nu/((1.0d0+nu)*(1.0d0-2.0d0*nu))
            
            ! create the reduced array of material parameters
            matpars(1)  = mu
            matpars(2)  = lam
            matpars(3)  = c
            matpars(4)  = alpha
            matpars(5)  = K0
            matpars(6)  = thetaZ
            matpars(7)  = dens
            matpars(8)  = betaD
            matpars(9)  = etaD
            matpars(10) = dD
            matpars(11) = gamma 
            matpars(12) = cD
            matpars(13) = cP
            
            !---------------------------------------------------------
            ! Read element information
            totElems = PROPS(15)
            totIntPoints = PROPS(16)
            
            !---------------------------------------------------------
            ! Run the considered material routine

            IF (elemForm == 0) THEN
                ! Run the mechanical UMAT Interface for the temperature field
                CALL UMAT_TEMP(DFGRD1,DFGRD0,TEMP+DTEMP,DTIME,STRESS,
     &                         RPL,DDSDDE,DDSDDT,DRPLDE,DRPLDT,tangent,
     &                         matpars,NOEL,NPT,totElems,totIntPoints,
     &                         STATEV)
            ELSE IF (elemForm == 1) THEN
                ! Run the mechanical UMAT Interface for the damage field
                CALL UMAT_DAMAGE(DFGRD1,DFGRD0,Temp+DTEMP,STATEV,STRESS,RPL,
     &                           DDSDDE,DDSDDT,DRPLDE,DRPLDT,PNEWDT,
     &                           matpars,NOEL,NPT,totElems,totIntPoints,
     &                           DTIME,tangent)
                
            ELSE
                WRITE(*,*) 'User Error:'
                WRITE(*,*) 'The flag for the master-slave decision is wrong'
            END IF
       
	    STATEV(5) = DFGRD1(1,1)
	    STATEV(6) = DFGRD1(1,2)
	    STATEV(7) = DFGRD1(1,3)
	    STATEV(8) = DFGRD1(2,1)
	    STATEV(9) = DFGRD1(2,2)
	    STATEV(10) = DFGRD1(2,3)
	    STATEV(11) = DFGRD1(3,1)
	    STATEV(12) = DFGRD1(3,2)
	    STATEV(13) = DFGRD1(3,3)
	    
        END SUBROUTINE


!-----------------------------------------------------------------------
!
!                 U M A T H T - S U B R O U T I N E
!
!-----------------------------------------------------------------------

        SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,STATEV,TEMP,DTEMP,
     &                    DTEMDX,TIME,DTIME,PREDEF,DPRED,CMNAME,NTGRD,
     &                    NSTATV,PROPS,NPROPS,COORDS,PNEWDT,NOEL,NPT,
     &                    LAYER,KSPT,KSTEP,KINC)
     
            IMPLICIT REAL*8(a-h,o-z)
        
            CHARACTER*80 CMNAME
            DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     &                DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     &                TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)
        
	        DOUBLE PRECISION, DIMENSION(3) :: y, Q
            DOUBLE PRECISION, DIMENSION(3,3) :: F
            DOUBLE PRECISION :: c, cond0, cD, cP, DU, dens, densT, detJ, DUDT, U
            INTEGER          :: elemForm,I
            
            !---------------------------------------------------------
            ! Read current element type
            elemForm  = PROPS(1)
            
            !---------------------------------------------------------
            ! Read material parameters for thermal problem
            
            ! specific heat capacity
            c         = PROPS(4)
            ! thermal conductivity
            cond0     = PROPS(6)
            dens      = PROPS(8)
            ! regularisation parameter
            cD        = PROPS(13)
            ! numerical heat capacity
            cP        = PROPS(14)

            F(1,1) = STATEV(5)
            F(1,2) = STATEV(6)
            F(1,3) = STATEV(7)
            F(2,1) = STATEV(8)
            F(2,2) = STATEV(9)
            F(2,3) = STATEV(10)
            F(3,1) = STATEV(11)
            F(3,2) = STATEV(12)
            F(3,3) = STATEV(13)
            CALL CALCDET33(F,detJ)
            densT = dens/detJ

            !---------------------------------------------------------
            ! Run the thermal material routine

            IF (elemForm == 0) THEN
                
                ! derivative of the internal energy w.r.t temperature
                DUDT = c/dens
                DUDG = 0.0d0
                
                ! Update Increment of thermal energy
                DU  = DUDT*DTEMP
                U   = U + DU
        
                ! Initialize derivative of flux w.r.t temperature
                DFDT = 0.0d0

		! Push operation of heat flux vector
		Q = MATMUL(DTEMDX, F)/detJ
                DO I=1, NTGRD
                    ! heat flux according to Fourier's law:
                    FLUX(I) = -cond0*Q(I)
                    ! derivative of heat flux w.r.t gradient of temp
                    DFDG(I,I) = -cond0
                END DO
        
            ELSE IF (elemForm == 1) THEN
                
                ! Calculate the damage quantaties
                ! Derivative of internal energy w.r.t phi
                DUDT = cP/dens
                DUDG = 0.0d0
                
                ! Increment of internal energy
                DU = DUDT*DTEMP
                ! Update of internal energy
                U = U+DU
                
                ! Initialize derivative of flux w.r.t phi
                DFDT = 0.0d0
                
		! Push operation of micromorphic flux
		y = MATMUL(DTEMDX, F)/detJ              
                DO I=1, NTGRD
                    ! damage flux
                    FLUX(I) = - cD*y(I)
                    ! derivative of damage flux w.r.t gradient phi
                    DFDG(I,I) = -cD
                END DO

          ELSE
                WRITE(*,*)  'User Error:'
                WRITE(*,*)  'The flag for the master-slave decision is wrong'
          END IF

        END SUBROUTINE
        
        
!------------------------------------------------------------------------
!
!         U M A T   I N T E R F A C E   R O U T I N E 
!
!------------------------------------------------------------------------
        
        SUBROUTINE UMAT_TEMP(F,Fn,theta,dt,STRESS,RPL,DDSDDE,DDSDDT,
     &                       DRPLDE,DRPLDT,tangent,matpars,elemNum,
     &                       intPNum,totElems,totIntP,STATEV)
            
            INCLUDE 'SMAAspUserArrays.hdr'
            !INCLUDE 'SMAAspUserSubroutines.hdr'
            
            ! Abaqus return variables
            DOUBLE PRECISION, DIMENSION(6)       :: STRESS,DDSDDT,DRPLDE
            DOUBLE PRECISION, DIMENSION(6,6)     :: DDSDDE
            DOUBLE PRECISION                     :: RPL,DRPLDT
            
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: F,Fn
            DOUBLE PRECISION, DIMENSION(13)      :: matpars
            DOUBLE PRECISION, DIMENSION(13)       :: STATEV
            DOUBLE PRECISION                     :: theta,dt
            INTEGER                              :: elemNum,intPNum
            INTEGER                              :: totElems,totIntP
            
            ! Internal variables
            DOUBLE PRECISION, DIMENSION(3,3,3,3) :: DDSDDE4T,DPDF
            DOUBLE PRECISION, DIMENSION(3,3,3,3) :: DsigmaDg
            DOUBLE PRECISION, DIMENSION(3,3)     :: FinvT,Finv,FT,Ident
            DOUBLE PRECISION, DIMENSION(3,3)     :: P,DRPLDF,DRPLDg
            DOUBLE PRECISION, DIMENSION(3,3)     :: DPDtheta,sigma
            DOUBLE PRECISION, DIMENSION(3,3)     :: DsigDtheta,DRPLDF_N
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV0
            DOUBLE PRECISION                     :: phi,detJ,d,dn,ddot,Y
            DOUBLE PRECISION                     :: PNEWDT,rmech,ddot0,DRPLDF_NORM
            INTEGER                              :: curIntP
            CHARACTER(LEN=3)                     :: tangent,LENOUTDIR
            CHARACTER(LEN=256)                   :: OUTDIR  
            
            ! Global arrays
            DOUBLE PRECISION, DIMENSION(totIntP) :: phiArray
            POINTER(ptrPhi,phiArray)
            DOUBLE PRECISION, DIMENSION(totIntP) :: dArray
            POINTER(ptrd,dArray)
            DOUBLE PRECISION, DIMENSION(totIntP) :: ddotArray
            POINTER(ptrddot,ddotArray)
            DOUBLE PRECISION, DIMENSION(totIntP) :: thetaArray
            POINTER(ptrTheta,thetaArray)
              
            ! Get the Work directory for the data transfer
            CALL GETOUTDIR(OUTDIR, LENOUTDIR)
            
            !-----------------------------------------------------------
            ! Read material parameters
            mu      = matpars(1)    
            lam     = matpars(2)
            c       = matpars(3)      
            alpha   = matpars(4) 
            K0      = matpars(5)    
            thetaZ  = matpars(6) 
            dens    = matpars(7)
            betaD   = matpars(8)
            etaD    = matpars(9)
            dD      = matpars(10)
            
            !---------------------------------------------------------
            ! Read and write the data transfer

            ! mapping of local to global int Point number
            curIntP = 8*(elemNum-1) + intPNum
            
            !IF (ptrPhi .EQ. 0 ) ptrPhi = SMALocalFloatArrayCreate(1,totIntP,0.0d0)
            IF (ptrPhi .EQ. 0 ) ptrPhi = SMAFloatArrayCreate(1,totIntP,0.0d0)
            !ptrPhi = SMALocalFloatArrayAccess(1)
            ptrPhi = SMAFloatArrayAccess(1)
            phi = phiArray(curIntP)
            
            !IF (ptrd .EQ. 0 ) ptrd = SMALocalFloatArrayCreate(2,totIntP,0.0d0)
            IF (ptrd .EQ. 0 ) ptrd = SMAFloatArrayCreate(2,totIntP,0.0d0)
            !ptrd = SMALocalFloatArrayAccess(2)
            ptrd = SMAFloatArrayAccess(2)
            d = dArray(curIntP)
              
            !IF (ptrddot .EQ. 0 ) ptrddot = SMALocalFloatArrayCreate(3,totIntP,0.0d0)
            IF (ptrddot .EQ. 0 ) ptrddot = SMAFloatArrayCreate(3,totIntP,0.0d0)
            !ptrddot = SMALocalFloatArrayAccess(3)
            ptrddot = SMAFloatArrayAccess(3)
            ddot = ddotArray(curIntP)
            dn = d-ddot*dt
                
            !IF (ptrTheta .EQ. 0 ) ptrTheta = SMALocalFloatArrayCreate(4,totIntP,thetaZ)
            IF (ptrTheta .EQ. 0 ) ptrTheta = SMAFloatArrayCreate(4,totIntP,thetaZ)
            !ptrTheta = SMALocalFloatArrayAccess(4)
            ptrTheta = SMAFloatArrayAccess(4)
            thetaArray(curIntP) = theta
            
            
            !---------------------------------------------------------
            ! Set up some required tensors
              
            ! Identity matrix
            CALL IDENTITY(Ident) 
			
            ! Jacobi-Determinant
            CALL CALCDET33(F,detJ)	
			
            ! Transpose of Deform. Gradient
            FT = TRANSPOSE(F)	
			
            ! Inverse of Deform. Gradient
            CALL CALCINV33(F,Finv)		
			
            ! Transpose of Inverse of Deformation Gradient
            FinvT = TRANSPOSE(Finv)
            
            STATEV0 = STATEV
            ddot0 = ddot
            ! Evaluate constitutive model
            CALL KConstModel(F,Fn,Finv,phi,theta,matpars,STATEV,PNEWDT,ddot,dt,Y,rmech,P)
            
            RPL = rmech/detJ
            !---------------------------------------------------------
            ! Calculate stresses
			
            ! Push-Forward
            sigma = MATMUL(F,P)/detJ
    
            ! transform 2nd order tensor to appropriate Abaqus Voigt format
            CALL TENSOR33TOVOIGT(sigma,STRESS)
             
            ! --------------------------------------------------------
            ! Calculate DDSDDE
            ! Derivative of stresses w.r.t. strains
            IF (tangent=='num') THEN
                CALL DPDF_NUM(F,Fn,phi,matpars,STATEV0,theta,ddot0,dt,DPDF)
                !CALL DPDF_ANA(F,Finv,theta,matpars,DPDF,STATEV0,phi)
            ELSE IF (tangent=='ana') THEN
                CALL DPDF_ANA(F,Finv,theta,matpars,DPDF,STATEV0,phi)
            ELSE
                WRITE(*,*) 'User Error:'
                WRITE(*,*) 'Choose a correct tangent type'
            END IF
            
            ! Push Forward in spatial metric
            CALL PUSHFORWARDDPDF(DPDF,F,sigma,DsigmaDg)
            
            ! Jaumann-Correction for Abaqus
            CALL JaumannCorrection(sigma,DsigmaDg,DDSDDE4T)
                                                            
            ! Transfor 4th order tensor in Voigt notation
            CALL T3333TOVOIGT(DDSDDE4T, DDSDDE)

            ! --------------------------------------------------------
            ! Calculate DDSDDT
              
            ! Derivative of stresse w.r.t temperature
            IF (tangent=='num') THEN
                CALL DPDT_NUM(F,Fn,phi,matpars,STATEV,STATEV0,theta,ddot0,dt,DPDTheta)
            ELSE IF (tangent=='ana') THEN
                CALL DPDT_ANA(F,Finv,matpars,DPDTheta,STATEV0,phi,theta)
            ELSE
                WRITE(*,*) 'User Error:'
                WRITE(*,*) 'Choose a correct tangent type'
            END IF
            
            
            ! push forward in spatial metric
            DsigDTheta = MATMUL(F, DPDTheta)/detJ
            
            ! convert to Abaqus Voigt format
            CALL TENSOR33TOVOIGT(DsigDTheta,DDSDDT)
            
            ! --------------------------------------------------------
            ! Calculate DRPLDT
              
            ! Derivative of the volumetric heat source w.r.t temperature
            IF (tangent=='num') THEN
                CALL DRDT_NUM(F,Fn,phi,matpars,STATEV,STATEV0,theta,ddot0,dt,DRPLDT)
            ELSE IF (tangent=='ana') THEN
                CALL DRDT_ANA(F,Fn,theta,matpars,dt,DRPLDT,d,dn,phi)
            ELSE
                WRITE(*,*) 'User Error:'
                WRITE(*,*) 'Choose a correct tangent type'
            END IF
            
            DRPLDT = DRPLDT/detJ
            
            ! --------------------------------------------------------
            ! Calculate DRPLDE
            
            ! Derivative of the volumetric heat source w.r.t strains
            IF (tangent=='num') THEN
                CALL DRDF_NUM(F,Fn,phi,matpars,STATEV,STATEV0,theta,ddot0,dt,DRPLDF_N)
            ELSE IF (tangent=='ana') THEN
                CALL DRDF_ANA(F,Finv,Fn,theta,matpars,dt,DRPLDF,STATEV,dn,phi)
            ELSE
                WRITE(*,*) 'User Error:'
                WRITE(*,*) 'Choose a correct tangent type'
            END IF

            
            ! push forward DRPLDC to obtain DRPLDg (g = spatial metrik)
            DRPLDg = MATMUL(F,DRPLDF)/detJ
            
            ! convert spatial DRPLDg to Abaqus Voigt format
            CALL TENSOR33TOVOIGT(DRPLDg, DRPLDE) 
            
        END SUBROUTINE

        SUBROUTINE UMAT_DAMAGE(F,Fn,phi,STATEV,STRESS,RPL,DDSDDE,DDSDDT,
     &                         DRPLDE,DRPLDT,PNEWDT,matpars,elemNum,
     &                         intPNum,totElems,totIntP,dt,tangent)
                            
              INCLUDE 'SMAAspUserArrays.hdr'
              !INCLUDE 'SMAAspUserSubroutines.hdr'
            
              ! Abaqus return variables
              DOUBLE PRECISION, DIMENSION(6)   :: STRESS,DDSDDT,DRPLDE
              DOUBLE PRECISION, DIMENSION(13)  :: STATEV
              DOUBLE PRECISION, DIMENSION(6,6) :: DDSDDE
              DOUBLE PRECISION                 :: RPL,DRPLDT,PNEWDT
            
              ! Input variables
              DOUBLE PRECISION, DIMENSION(3,3) :: F,Fn
              DOUBLE PRECISION, DIMENSION(13)  :: matpars
              DOUBLE PRECISION                 :: phi,dt
              INTEGER                          :: elemNum,intPNum
              INTEGER                          :: totElems,totIntP
            
              ! Internal variables 
              DOUBLE PRECISION, DIMENSION(3,3) :: DRPLDF,DRPLDg,Finv,P
              DOUBLE PRECISION, DIMENSION(3,3) :: DPDPhi,DsigDPhi
              DOUBLE PRECISION, DIMENSION(13)  :: STATEV0
              DOUBLE PRECISION                 :: theta,ddot,dn,Y,rmech
              DOUBLE PRECISION                 :: ddot0,detJ
              INTEGER                          :: locElemNum
              INTEGER                          :: curIntP,LENOUTDIR
              CHARACTER(LEN=3)                 :: tangent
            
              ! Data transfer variables
              CHARACTER(LEN=256)                   :: OUTDIR
              DOUBLE PRECISION, DIMENSION(totIntP) :: phiArray
              POINTER(ptrPhi,phiArray)
              DOUBLE PRECISION, DIMENSION(totIntP) :: dArray
              POINTER(ptrd,dArray)
              DOUBLE PRECISION, DIMENSION(totIntP) :: ddotArray
              POINTER(ptrddot,ddotArray)
              DOUBLE PRECISION, DIMENSION(totIntP) :: thetaArray
              POINTER(ptrTheta,thetaArray)
            
              ! Get the Output directory
              CALL GETOUTDIR(OUTDIR, LENOUTDIR)
            
              !---------------------------------------------------------
              ! Read and write data transfer

              ! mapping of local int Point number to global int Point number
              locElemNum = elemNum - totElems
              curIntP = 8*(locElemNum-1) + intPNum
            
              !IF (ptrPhi .EQ. 0 ) ptrPhi = SMALocalFloatArrayCreate(1,totIntP,0.0d0)
              IF (ptrPhi .EQ. 0 ) ptrPhi = SMAFloatArrayCreate(1,totIntP,0.0d0)
              !ptrPhi = SMALocalFloatArrayAccess(1)
              ptrPhi = SMAFloatArrayAccess(1)
              phiArray(curIntP) = phi
        
              thetaZ = matpars(6)
              !IF (ptrTheta .EQ. 0 ) ptrTheta = SMALocalFloatArrayCreate(4,totIntP,thetaZ)
              IF (ptrTheta .EQ. 0 ) ptrTheta = SMAFloatArrayCreate(4,totIntP,thetaZ)
              !ptrTheta = SMALocalFloatArrayAccess(4)
              ptrTheta = SMAFloatArrayAccess(4)
              theta = thetaArray(curIntP)            
            
              !---------------------------------------------------------
              ! evaluate constitutive response
            
              ! inverse of F
              CALL CALCINV33(F, Finv)
              ! Jacobi-Determinant
              CALL CALCDET33(F,detJ)	
            
            
              STATEV0 = STATEV
              ddot0 = ddot
              CALL KConstModel(F,Fn,Finv,phi,theta,matpars,STATEV,PNEWDT,ddot,dt,Y,rmech,P)
              
              RPL = Y

              ! Even though, there is a stress result for the local const
              ! model, the stresses will be set to zero for the slave
              ! eleme. The stresses will be calculated for the Master elem
              STRESS = 0.0d0

              CALL DPDPhi_NUM(F,Fn,phi,theta,matpars,STATEV,STATEV0,ddot,dt,DPDPhi)
              
              ! push forward in spatial metric
              DsigDPhi = MATMUL(F, DPDPhi)/detJ
            
              ! convert to Abaqus Voigt format
              CALL TENSOR33TOVOIGT(DsigDPhi,DDSDDT)
              
              ! Also the mechanical tangent modulus will be set to zero for
              ! the slave elem
              DDSDDE = 0.0d0
            
              !---------------------------------------------------------
              ! Calculate DRPLDT
              IF (tangent=='num') THEN
                CALL DYDPhi_NUM(F,Fn,phi,matpars,STATEV,STATEV0,theta,ddot0,dt,DRPLDT)
              ELSE IF (tangent=='ana') THEN
                CALL DYDPhi_ANA(F,matpars,phi,theta,STATEV,DRPLDT)
              ELSE
                WRITE(*,*) 'User Error:'
                WRITE(*,*) 'Choose a correct tangent type'
              END IF
            
              DRPLDT = DRPLDT
                  
              !---------------------------------------------------------
              ! Calculate DRPLDF
              IF (tangent=='num') THEN
                !CALL DYDF_NUM(F,Fn,phi,matpars,STATEV,STATEV0,theta,ddot0,dt,DRPLDF)
                CALL DYDF_ANA(matpars,STATEV0,F,Finv,phi,theta,DRPLDF)
              ELSE IF (tangent=='ana') THEN
                CALL DYDF_ANA(matpars,STATEV0,F,Finv,phi,theta,DRPLDF)
              ELSE
                WRITE(*,*) 'User Error:'
                WRITE(*,*) 'Choose a correct tangent type'
              END IF
            
            
              ! Push forward DRPLDF to obtain DRPLDg (g = spatial metrik)
              DRPLDg = MATMUL(F,DRPLDF)
            
              ! convert spatial DRPLDg to Abaqus Voigt format
              CALL TENSOR33TOVOIGT(DRPLDg, DRPLDE)
            
              !---------------------------------------------------------
              ! Read and write local damage variable
            
              IF (ptrd .EQ. 0 ) ptrd = SMAFloatArrayCreate(2,totIntP,0.0d0)
              ptrd = SMAFloatArrayAccess(2)
              dArray(curIntP) = STATEV(1)

              IF (ptrddot .EQ. 0 ) ptrddot = SMAFloatArrayCreate(3,totIntP,0.0d0)
              ptrddot = SMAFloatArrayAccess(3)
              ddotArray(curIntP) = ddot
              
          END SUBROUTINE UMAT_DAMAGE
!------------------------------------------------------------------------
!
!  L O C A L  T E M P E R A T U R E  C O N S T I T U T I V E   M O D E L
!
!------------------------------------------------------------------------

        ! Constitutive Routine to calculate 1st Piola stresses
        SUBROUTINE CALCP(F,Finv,matpars,theta,P,fd)
            
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3) :: F,Finv
            DOUBLE PRECISION, DIMENSION(13)  :: matpars
            DOUBLE PRECISION                 :: theta,fd
        
            ! Output variables
            DOUBLE PRECISION, DIMENSION(3,3) :: P
            
            ! Internal variables
            DOUBLE PRECISION, DIMENSION(3,3) :: FinvT,Pela
            DOUBLE PRECISION                 :: mu,lam,c,alpha,K0,thetaZ 
            DOUBLE PRECISION                 :: dens,betaD,etaD,dD,K
            DOUBLE PRECISION                 :: detJ
                       
            !-----------------------------------------------------------
            ! read material parameters
            mu      = matpars(1)    
            lam     = matpars(2)
            c       = matpars(3)      
            alpha   = matpars(4) 
            K0      = matpars(5)    
            thetaZ  = matpars(6) 
            dens    = matpars(7)
            betaD   = matpars(8)
            etaD    = matpars(9)
            dD      = matpars(10)
            
            !-----------------------------------------------------------
            ! Set up some require quantities
            
            ! Jacobi-Determinant
            CALL CALCDET33(F, detJ)	
			
            ! Transpose of Inverse of Deformation Gradient
            FinvT = TRANSPOSE(Finv)

            ! bulk modulus
            K = lam + 2.0d0/3.0d0*mu
                        
            !-----------------------------------------------------------
            ! Calculate the Piola stresses
            
            ! Calculate the elastic Piola stresses
            Pela = (lam*LOG(detJ)-mu)*FinvT + mu*F 
     &  - 3.0d0/detJ*alpha*K*(theta-thetaZ)*(1.0d0-LOG(detJ))*FinvT
            
            ! Calculate the inelastic Piola stresses
            P = fd*Pela
            
        END SUBROUTINE

        ! Constitutive Routine to calculate volumetric heat source
        SUBROUTINE CALCRPL(F,Fn,theta,matpars,dt,RPL,d,dn,phi)
            
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3) :: F,Fn
            DOUBLE PRECISION, DIMENSION(13)  :: matpars
            DOUBLE PRECISION                 :: dt,theta,d,dn,phi
            
            ! Output variables
            DOUBLE PRECISION                 :: RPL
            
            ! Internal variables
            DOUBLE PRECISION                 :: detJ,Jn,Jdot,rEla,fd,IH
            DOUBLE PRECISION                 :: mu,lam,c,alpha,K0,dens
            DOUBLE PRECISION                 :: etaD,betaD,thetaZ,K,dD
            DOUBLE PRECISION                 :: gamma,GJ,DfdDd,PsiEla
            DOUBLE PRECISION                 :: DPsiElaDtheta,FDDF,ddot
            DOUBLE PRECISION                 :: DDfdDDd,PsiMech
            DOUBLE PRECISION                 :: DPsiMechDtheta
            
            !---------------------------------------------------------
            ! read material parameters
            
            mu     = matpars(1)
            lam    = matpars(2)
            c      = matpars(3)      
            alpha  = matpars(4) 
            K0     = matpars(5)    
            thetaZ = matpars(6) 
            dens   = matpars(7)
            betaD  = matpars(8)
            etaD   = matpars(9)
            dD     = matpars(10)
            gamma  = matpars(11)

            !---------------------------------------------------------
            ! Set up some required quantaties
        
            CALL CALCDET33(F,detJ)
        
            ! Jacobian of previous time step
            CALL CALCDET33(Fn,Jn)
            
            ! Rate of Jacobian
            Jdot = (detJ-Jn)/dt
              
            ! Rate of damage
            ddot = (d-dn)/dt
            
            ! bulk modulus
            K = lam + 2.0d0/3.0d0*mu
              
            ! Double contraction of F with F
            CALL T2DDT233(F,F,FDDF)
            
            ! Calculate the damage function
            CALL DAMAGEFUNC(etaD,dD,d,fd,DfdDd,DDfdDDd)
            
            !---------------------------------------------------------
            ! Gough Joule contribution:
              
            ! elastic heat source
            rEla = -3.0d0*theta*alpha*K*detJ**(-2.0d0)*(1.0d0-LOG(detJ))*Jdot
            
            ! Gough-Joule
            GJ = fd*rEla
              
            !-----------------------------------------------------------
            ! Inelastic heat contribution:

            ! Calculate the elastic free energy
            PsiEla = mu/2.0d0*(FDDF-3.0d0) - mu*LOG(detJ) 
     &               + lam/2.0d0*LOG(detJ)**(2.0d0)
     &               - 3.0d0*alpha*K*(theta-thetaZ)*LOG(detJ)/detJ
     &               + c*(theta-thetaZ-theta*LOG(theta/thetaZ))
            
            PsiMech = mu/2.0d0*(FDDF-3.0d0) - mu*LOG(detJ) 
     &                + lam/2.0d0*LOG(detJ)**(2.0d0)
     &                - 3.0d0*alpha*K*(theta-thetaZ)*LOG(detJ)/detJ
     
            ! Derivative of Psi_ela w.r.t. theta
            DPsiElaDtheta = -c*LOG(theta/thetaZ)-3.0d0*alpha*K*LOG(detJ)/detJ
            DPsiMechDtheta = -3.0d0*alpha*K*LOG(detJ)/detJ
            
            ! derivative of damage function
            DfdDd = -etaD*fd
              
            ! Inelastic heat
            IH = (theta*DfdDd*DPsiMechDtheta-DfdDd*PsiMech)*ddot
              
              
            !---------------------------------------------------------
            ! Calculate volumetric heat source
              
            RPL = (GJ + IH)*dens
              
        END SUBROUTINE

!-----------------------------------------------------------------------
!
!         L O C A L  D A M A G E  C O N S T I T U T I V E   M O D E L
!
!-----------------------------------------------------------------------
          
        ! Definition of damage related local constitutive model
        SUBROUTINE KConstModel(F,Fn,Finv,phi,theta,matpars,STATEV,PNEWDT,ddot,dt,Y,rmech,P)
          
            ! Finv und FinvT und Ident raus 
            
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3) :: F,Finv,Fn
            DOUBLE PRECISION, DIMENSION(13)  :: matpars
            DOUBLE PRECISION                 :: phi,theta,dt
            
            ! Output variables
            DOUBLE PRECISION, DIMENSION(3,3) :: P
            DOUBLE PRECISION, DIMENSION(13)  :: STATEV
            DOUBLE PRECISION                 :: PNEWDT,ddot,Y,rmech
            
            ! Internal variables
            DOUBLE PRECISION                 :: mu,lam,c,alpha,K0,dens
            DOUBLE PRECISION                 :: thetaZ,betaD,etaD,dD,K
            DOUBLE PRECISION                 :: gamma,fd,detJ,dn,dk,FDDF
            DOUBLE PRECISION                 :: q,DfdDd,PsiEla,res,dLam
            DOUBLE PRECISION                 :: resMin,DDfdDDd,dResMin
            DOUBLE PRECISION                 :: DRDLambda,nd,PsiMech
            INTEGER                          :: newtonIter
            
            !---------------------------------------------------------
            ! Read material parameters
                        
            mu     = matpars(1)
            lam    = matpars(2)
            c      = matpars(3)      
            alpha  = matpars(4) 
            K0     = matpars(5)    
            thetaZ = matpars(6) 
            dens   = matpars(7)
            betaD  = matpars(8)
            etaD   = matpars(9)
            dD     = matpars(10)
            gamma  = matpars(11)
          
            !---------------------------------------------------------
            ! Set up some required qunataties
                                    
            ! determinant of J
            CALL CALCDET33(F, detJ)
            
            ! Double contraction of F with F
            CALL T2DDT233(F,F,FDDF)
                 
            ! bulk modulus
            K = lam + 2.0d0/3.0d0*mu
            
            !---------------------------------------------------------
            ! Calculate the damage driving force
            
            ! read scalar damage variable kappa from sdv
            dn = STATEV(1);
        
            ! calculate damage function
            CALL DAMAGEFUNC(etaD,dD,dn,fd,DfdDd,DDfdDDd)
              
            ! Calculate the elastic free energy
            PsiEla = mu/2.0d0*(FDDF - 3.0d0) - mu*LOG(detJ) 
     &  + lam/2.0d0*LOG(detJ)**(2.0d0)
     &  - 3.0d0*alpha*K*(theta-thetaZ)*LOG(detJ)/detJ
     &  + c*(theta-thetaZ-theta*LOG(theta/thetaZ))
            
            PsiMech = mu/2.0d0*(FDDF - 3.0d0) - mu*LOG(detJ) 
     &  + lam/2.0d0*LOG(detJ)**(2.0d0)
     &  - 3.0d0*alpha*K*(theta-thetaZ)*LOG(detJ)/detJ
     
            ! calculate damage driving force
            q = PsiMech + betaD*gamma*(phi-gamma*dn)/(etaD*fd)
            !-----------------------------------------------------------  
            ! Predictor corrector scheme
            
            ! check damage condition
            dk = dn
            dLam = 0.0d0
            
            STATEV(3) = 0.0d0
            STATEV(4) = q
        
            
            IF (q - dn > 0.0d0) THEN

                ! inelastic state
                resMin = 1.0d0
                
                ! Newton-Raphso-scheme
                DO newtonIter = 1,20
                    
                    ! calculate current damage function
                    CALL DAMAGEFUNC(etaD,dD,dk,fd,DfdDd,DDfdDDd)
                    
                    ! calculate damage driving force
                    q = PsiMech + betaD*gamma*(phi-gamma*dk)/(etaD*fd)
                    
                    ! update residual function value
                    res = q - dk
                    
                    ! compute tangent of residual, i.e. DPhidDkappa
                    DRDLambda = betaD*(-(phi-dk)*DfdDd-1.0d0)/(etaD*fd)-1.0d0
                    
                    ! calculate Newton update
                    dLam = dLam - res/DRDLambda
                    dk = dn + dLam
                    
                    ! Break-ou<t criterion
                    IF (abs(res) < 10.0d0**(-8.0d0)) EXIT
                    
                    ! try to improve local algorithm's stability by tracing the lowest residual
                    ! function value and the associated value for kappa => kapparkmin
                    IF (abs(res) < abs(resMin)) THEN
                        resMin = res
                        dResMin = dk
                    END IF
                    
                END DO
                
                !-------------------------------------------------------
                ! If local Newton scheme did not converge, request a lower
                ! time increment from ABAQUS
                IF (newtonIter > 19) THEN
                    ! Residuum lower than 1e-5 is accepted without throwing a warning
					! in this case, dont change the proposed time increment (PNEWDT = 1)
                    PNEWDT = 1.0d0
                    IF (abs(resMin) > 10.0d0**(-5.0d0)) THEN
                        WRITE(*,*) 'Local problem! Lowest res/d:',
     &                              resMin, dResMin
                        ! we didnt converge, so use the specific value of d
                        ! that was associated to the lowest residual to proceed:
                        dk = dResMin
                        ! lower the proposed time increment
                        PNEWDT = 0.5d0
                        WRITE(*,*) 'Requesting finer timestep: PNEWDT=', PNEWDT
                    END IF
                END IF
                ! Update state variable for the inealstic case
                STATEV(1) = dk
            END IF
            
            ! after local convergence: get final values of the damage function fd
            CALL DAMAGEFUNC(etaD,dD,dk,fd,DfdDd,DDfdDDd)
            
            ! save the damage function 'f_d' as SDV for Postprocessing
            STATEV(2) = 1.0d0 - fd
            ddot = (dk-dn)/dt
            
            !-----------------------------------------------------------
            ! Calculate damage source term
            Y = -betaD*(phi - gamma*dk)
            !STATEV(3) = Y
            !-----------------------------------------------------------
            ! Calculate the Piola stresses
            CALL CALCP(F,Finv,matpars,theta,P,fd)

            !---------------------------------------------------------
            ! Calculate volumetric heat generation
            CALL CALCRPL(F,Fn,theta,matpars,dt,rmech,STATEV(1),dn,phi)
          
	    STATEV(3) = rmech
        END SUBROUTINE

        ! Definition of damage function
        SUBROUTINE DAMAGEFUNC(etaD,dD,d,fd,DfdDd,DDfdDDd)
            
            ! Input variables
            DOUBLE PRECISION :: etaD,d,dD
        
            ! Output variables
            DOUBLE PRECISION :: fd,DfdDd,DDfdDDd,H
        
            ! evaluate damage function and derivatives
            fd      = exp(-etaD*max(0.0d0,d-dD))
    
            IF (d-dD > 0.0d0) THEN
                H = 1.0d0
            ELSE
                H = 0.0d0
            END IF
            
            DfdDd   = -etaD*fd*H
            DDfdDDd = etaD**(2.0d0)*fd*H
            
            
          END SUBROUTINE

        SUBROUTINE DAMGEFUNC2(etaD,dD,d,fd,DfdDd,DDfdDDd)
        
            DOUBLE PRECISION :: etaD,d,dD
        
            DOUBLE PRECISION :: fd,DfdDd,DDfdDDd,dMax
        
            dMax = 0.2d0
            fd = 1.0d0 - dMax*(1.0d0 - exp(-etaD*d))
        
            DfdDd = -dMax*etaD*exp(-etaD*d)
            DDfdDDd = dMax*etaD**(2.0d0)*exp(-etaD*d)
        
        END SUBROUTINE
!------------------------------------------------------------------------
!
!             D A M A G E   R E L A T E D   T A N G E N T S
!
!------------------------------------------------------------------------

        ! Tangent of the volumetric source w.r.t. the global damage var
        SUBROUTINE DYDPhi_ANA(F,matpars,phi,theta,STATEV,DRPLDPHI)
          
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3) :: F
            DOUBLE PRECISION, DIMENSION(13)  :: matpars
            DOUBLE PRECISION                 :: phi,theta
            DOUBLE PRECISION, DIMENSION(13)  :: STATEV
        
            ! Output variables
            DOUBLE PRECISION                 :: DRPLDPHI
        
            ! Internal variables
            DOUBLE PRECISION                 :: mu,lam,c,alpha,K0,thetaZ
            DOUBLE PRECISION                 :: dens,betaD,etaD,dD,gamma
            DOUBLE PRECISION                 :: detJ,FDDF,K,d,fd,DDfdDDd
            DOUBLE PRECISION                 :: PsiEla,DYDphi,DYDd,dddphi
            DOUBLE PRECISION                 :: DfdDd
            !-------------------------------------------------------------
            ! read material parameters
            mu     = matpars(1)
            lam    = matpars(2)
            c      = matpars(3)      
            alpha  = matpars(4) 
            K0     = matpars(5)    
            thetaZ = matpars(6) 
            dens   = matpars(7)
            betaD  = matpars(8)
            etaD   = matpars(9)
            dD     = matpars(10)
            gamma  = matpars(11)
            
            !---------------------------------------------------------
            ! Set up some required qunataties
              
            ! determinant of J
            CALL CALCDET33(F, detJ)
            
            ! Double contraction of F with F
            CALL T2DDT233(F,F,FDDF)
                 
            ! bulk modulus
            K = lam + 2.0d0/3.0d0*mu
            
            ! Calculate the elastic free energy
            PsiEla = mu/2.0d0*(FDDF - 3.0d0) - mu*LOG(detJ) 
     &             + lam/2.0d0*LOG(detJ)**(2.0d0)
     &             - 3.0d0*alpha*K*(theta-thetaZ)*LOG(detJ)/detJ
     &             + c*(theta-thetaZ-theta*LOG(theta/thetaZ))
     
            !-------------------------------------------------------------
            ! Calculate derivative of RPL w.r.t. phi
          
            ! read current local damage variable
            d = STATEV(1)
          
            ! calculate damage function
            CALL DAMAGEFUNC(etaD,dD,d,fd,DfdDd,DDfdDDd)
            
            DYDphi  = -betaD
            DYDd    = betaD*gamma
            
            IF (d == 0.0d0) THEN
                dddphi = 0.0d0
            ELSE
                dddphi  = -betaD*gamma/(DDfdDDd*PsiEla+betaD*gamma**(2.0d0))!-nd*DfdDd*dD*(1.0d0-fd)**(nd-1.0d0))
            END IF
            
            ! Calculate tangent
            DRPLDPHI = DYDphi + DYDd*dddphi
            
        END SUBROUTINE

        ! Tangent of the volumetric source w.r.t. the global damage var
        SUBROUTINE DYDPhi_NUM(F,Fn,phi,matpars,STATEV,STATEVN,theta,ddot,dt,dYdPhi)
          
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: F,Fn
            DOUBLE PRECISION, DIMENSION(13)      :: matpars
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV,STATEVN
            DOUBLE PRECISION                     :: phi,theta,ddot,dt
            
            ! Output variables
            DOUBLE PRECISION                     :: dYdPhi
            
            ! Internal variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: P,Finv
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV0
            DOUBLE PRECISION                     :: eps,PNEWDT,rmech
            DOUBLE PRECISION                     :: dPhi1,dY1,Y
            DOUBLE PRECISION                     :: ddot0
            INTEGER                              :: i,j
                         
              
            eps = 10.0d0**(-8.0d0)
            
            dPhi1 = phi + eps
            
            ddot0 = ddot
            STATEV0 = STATEVN
            CALL CALCINV33(F,Finv)
            CALL KConstModel(F,Fn,Finv,phi,theta,matpars,STATEV0,PNEWDT,ddot0,dt,Y,rmech,P)
            
            STATEV0 = STATEVN
            ddot0 = ddot
            CALL KConstModel(F,Fn,Finv,dPhi1,theta,matpars,STATEV0,PNEWDT,ddot0,dt,dY1,rmech,P)
                    
            dYdPhi = (dY1-Y)/(eps)
                
        END SUBROUTINE
          
        ! Tangent of the volumetric source w.r.t F
        SUBROUTINE DYDF_ANA(matpars,STATEV,F,Finv,phi,theta,DRPLDF)
          
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3) :: F,Finv
            DOUBLE PRECISION, DIMENSION(13)  :: matpars
            DOUBLE PRECISION, DIMENSION(3)   :: STATEV  
            DOUBLE PRECISION                 :: phi,theta
              
            ! Output variables
            DOUBLE PRECISION, DIMENSION(3,3) :: DRPLDF
          
            ! Internal variables
            DOUBLE PRECISION, DIMENSION(3,3) :: FinvT,Pela,dddF
            DOUBLE PRECISION                 :: mu,lam,c,alpha,K0,thetaZ
            DOUBLE PRECISION                 :: dens,betaD,etaD,dD,gamma
            DOUBLE PRECISION                 :: K,fd,detJ,d,FDDF,DfdDd 
            DOUBLE PRECISION                 :: DDfdDDd,PsiEla,DYDd
              
            !-------------------------------------------------------------
            ! read material parameters
          
            mu     = matpars(1)
            lam    = matpars(2)
            c      = matpars(3)      
            alpha  = matpars(4) 
            K0     = matpars(5)    
            thetaZ = matpars(6) 
            dens   = matpars(7)
            betaD  = matpars(8)
            etaD   = matpars(9)
            dD     = matpars(10)
            gamma  = matpars(11)
    
            !-------------------------------------------------------------
            ! Set up some required quantaties
            
            ! transpose of inverse
            FinvT = transpose(Finv)
          
            ! determinant of F
            CALL CALCDET33(F, detJ)
     
            ! bulk modulus
            K = lam + 2.0d0/3.0d0*mu
              
            ! Double contraction of F with F
            CALL T2DDT233(F,F,FDDF)
            
            ! Calculate the elastic free energy
            PsiEla = mu/2.0d0*(FDDF - 3.0d0) - mu*LOG(detJ) 
     &               + lam/2.0d0*LOG(detJ)**(2.0d0)
     &               - 3.0d0*alpha*K*(theta-thetaZ)*LOG(detJ)/detJ
     &               + c*(theta-thetaZ-theta*LOG(theta/thetaZ))
     
            ! Calculate the elastic Piola stresses
            Pela = (lam*LOG(detJ)-mu)*FinvT + mu*F 
     &  - 3.0d0/detJ*alpha*K*(theta-thetaZ)*(1.0d0-LOG(detJ))*FinvT
     
            !-------------------------------------------------------------
            ! calculate derivative of DRPL w.r.t. F
          
            ! read current local damage variable
            d = STATEV(1)
          
            ! Calculate damage function
            CALL DAMAGEFUNC(etaD,dD,d,fd,DfdDd,DDfdDDd)
            DYDd    = betaD*gamma
            
            ! Calculate the derivative of d w.r.t. F
            dddF = -etaD*fd*Pela/(betaD*(phi-d)*etaD-betaD-etaD*fd)
            
            ! Calculate the tangent
            DRPLDF = DYDd*dddF
            
        END SUBROUTINE
          
        ! Tangent of the volumetric source w.r.t. the global damage var
        SUBROUTINE DYDF_NUM(F,Fn,phi,matpars,STATEV,STATEVN,theta,ddot,dt,dYdF)
          
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: F,Fn
            DOUBLE PRECISION, DIMENSION(13)      :: matpars
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV,STATEVN
            DOUBLE PRECISION                     :: phi,theta,ddot,dt
            
            ! Output variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: dYdF
            
            ! Internal variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: deltaF1
            DOUBLE PRECISION, DIMENSION(3,3)     :: deltaFinv1,Finv
            DOUBLE PRECISION, DIMENSION(3,3)     :: P
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV0
            DOUBLE PRECISION                     :: eps,PNEWDT,rmech
            DOUBLE PRECISION                     :: ddot0,dY1,Y
            INTEGER                              :: i,j
                         
              
            eps = 10.0d0**(-8.0d0)
            
            CALL CALCINV33(F,Finv)
            ddot0 = ddot
            STATEV0 = STATEVN
            CALL KConstModel(F,Fn,Finv,phi,theta,matpars,STATEV0,PNEWDT,ddot0,dt,Y,rmech,P)
            
            dYdF = 0.0d0
            Do i = 1,3
                DO j = 1,3
                    
                    deltaF1 = F
                    deltaF1(i,j) = deltaF1(i,j) + eps
                
                    CALL CALCINV33(deltaF1,deltaFinv1)
                    
                    STATEV0 = STATEVN
                    ddot0 = ddot
                    
                    CALL KConstModel(deltaF1,Fn,deltaFinv1,phi,theta,matpars,STATEV0,PNEWDT,ddot0,dt,dY1,rmech,P)
                    
                    dYdF(i,j) = dYdF(i,j) + (dY1-Y)/(eps)
                
                END DO
              END DO
            
        END SUBROUTINE

        ! Tangent of the volumetric source w.r.t. the global damage var
        SUBROUTINE DPDPhi_NUM(F,Fn,phi,theta,matpars,STATEV,STATEVN,ddot,dt,dPdPhi)
          
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: F,Fn
            DOUBLE PRECISION, DIMENSION(13)      :: matpars
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV,STATEVN
            DOUBLE PRECISION                     :: phi,dt,theta,ddot
            
            ! Output variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: dPdPhi
            
            ! Internal variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: dP1,P,Finv
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV0
            DOUBLE PRECISION                     :: eps,PNEWDT,Y,ddot0
            DOUBLE PRECISION                     :: deltaPhi1,rmech                   
              
            eps = 10.0d0**(-8.0d0)
            
            CALL CALCINV33(F,Finv)
                            
            dPdPhi = 0.0d0
            
            deltaPhi1 = phi + eps
            
            STATEV0 = STATEVN
            ddot0 = ddot
            CALL KConstModel(F,Fn,Finv,phi,theta,matpars,STATEV0,PNEWDT,ddot0,dt,Y,rmech,P)
            
            STATEV0 = STATEVN
            ddot0 = ddot
            CALL KConstModel(F,Fn,Finv,deltaPhi1,theta,matpars,STATEV0,PNEWDT,ddot0,dt,Y,rmech,dP1)
                    
            dPdPhi = (dP1-P)/eps
                        
        END SUBROUTINE
        
!------------------------------------------------------------------------
!
!        T E M P E R A T U R E   R E L A T E D   T A N G E N T S
!
!------------------------------------------------------------------------

        ! Tangent of the volumetric source w.r.t. the global damage var
        SUBROUTINE DPDF_NUM(F,Fn,phi,matpars,STATEV,theta,ddot,dt,dPdF)
          
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: F,Fn
            DOUBLE PRECISION, DIMENSION(13)      :: matpars
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV
            DOUBLE PRECISION                     :: phi,theta,ddot,dt
            
            ! Output variables
            DOUBLE PRECISION, DIMENSION(3,3,3,3) :: dPdF
            
            ! Internal variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: deltaF1
            DOUBLE PRECISION, DIMENSION(3,3)     :: deltaFinv1,Finv
            DOUBLE PRECISION, DIMENSION(3,3)     :: dP1,P
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV0
            DOUBLE PRECISION                     :: eps,PNEWDT,Y,rmech
            DOUBLE PRECISION                     :: ddot0
            INTEGER                              :: i,j
                         
              
            eps = 10.0d0**(-8.0d0)
            ddot0 = ddot
            STATEV0 = STATEV
            CALL CALCINV33(F,Finv)
            CALL KConstModel(F,Fn,Finv,phi,theta,matpars,STATEV0,PNEWDT,ddot0,dt,Y,rmech,P)
            
            dPdF = 0.0d0
            Do i = 1,3
                DO j = 1,3
                    
                    deltaF1 = F
                    deltaF1(i,j) = deltaF1(i,j) + eps
                
                    CALL CALCINV33(deltaF1,deltaFinv1)
                    
                    STATEV0 = STATEV
                    ddot0 = ddot
                    
                    CALL KConstModel(deltaF1,Fn,deltaFinv1,phi,theta,matpars,STATEV0,PNEWDT,ddot0,dt,Y,rmech,dP1)
                    
                    dPdF(:,:,i,j) = dPdF(:,:,i,j) + (dP1-P)/(eps)
                
                END DO
              END DO
            
        END SUBROUTINE
        
        ! Calculate the analytical tangent
        SUBROUTINE DPDF_ANA(F,Finv,theta,matpars,DPDF,STATEV,phi)

            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: F,Finv
            DOUBLE PRECISION, DIMENSION(13)      :: matpars
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV
            DOUBLE PRECISION                     :: theta,phi
            
            ! Output variables
            DOUBLE PRECISION, DIMENSION(3,3,3,3) :: DPDF

            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3,3,3) :: FITDYADFIT,FITUDYADFI
            DOUBLE PRECISION, DIMENSION(3,3,3,3) :: IODYADI,PDYADdddF
            DOUBLE PRECISION, DIMENSION(3,3,3,3) :: DPelaDF,DPDFPartial
            DOUBLE PRECISION, DIMENSION(3,3)     :: FT,FinvT,Pela,Ident
            DOUBLE PRECISION, DIMENSION(3,3)     :: dddF
            DOUBLE PRECISION                     :: mu,lam,c,alpha,K0,dD
            DOUBLE PRECISION                     :: thetaZ,dens,betaD,K
            DOUBLE PRECISION                     :: etaD,gamma,detJ,fd
            DOUBLE PRECISION                     :: FDDF,PsiEla,DfdDd
            DOUBLE PRECISION                     :: DDfdDDd,d

            !---------------------------------------------------------
            ! Read material parameters
            
            mu      = matpars(1)    
            lam     = matpars(2)
            c       = matpars(3)      
            alpha   = matpars(4) 
            K0      = matpars(5)    
            thetaZ  = matpars(6) 
            dens    = matpars(7)
            betaD   = matpars(8)
            etaD    = matpars(9)
            dD      = matpars(10)
            gamma   = matpars(11)
            
            !---------------------------------------------------------
            ! Set up some required quantaties
            d = STATEV(1)
            
            ! Identity matrix
            CALL IDENTITY(Ident)

            ! determinant of F
            CALL CALCDET33(F, detJ)
            
            ! transpose of F	
            FT = TRANSPOSE(F)	
            
            ! transpose of inverse of F
            FinvT = TRANSPOSE(Finv)
            
            ! bulk modulus
            K = lam + 2.0d0/3.0d0*mu
              
            ! Double contraction of F with F
            CALL T2DDT233(F,F,FDDF)

            ! Calculate damage function
            CALL DAMAGEFUNC(etaD,dD,d,fd,DfdDd,DDfdDDd)
            
            !---------------------------------------------------------
            ! Calculate some damage related quantaties
            
            ! Calculate the elastic Piola stresses
            Pela = (lam*LOG(detJ)-mu)*FinvT + mu*F 
     &  - 3.0d0/detJ*alpha*K*(theta-thetaZ)*(1.0d0-LOG(detJ))*FinvT
     
            ! Calculate the elastic free energy
            PsiEla = mu/2.0d0*(FDDF - 3.0d0) - mu*LOG(detJ) 
     &               + lam/2.0d0*LOG(detJ)**(2.0d0)
     &               - 3.0d0*alpha*K*(theta-thetaZ)*LOG(detJ)/detJ
     &               + c*(theta-thetaZ-theta*LOG(theta/thetaZ))
     
                        
            ! Calculate the derivative of d w.r.t. F
            dddF = -DfdDd*fd*etaD*Pela/(betaD*(phi-d)*etaD-betaD-etaD*fd)

            !---------------------------------------------------------
            ! Calculate some required tensor products
            
            ! Calculate standard dyadic product of FinvT w. FinvT
            CALL T2T2DYAD33(FinvT,FinvT,FITDYADFIT)
            
            ! Calculate the lower dyadic Product of FinvTrans w. Finv
            CALL T2T2UDYAD33(FinvT,Finv,FITUDYADFI)
            
            ! Calculate the upper dyadic product of Ident w. Ident
            CALL T2T2ODYAD33(Ident,Ident,IODYADI)
            
            ! Calculate the standard dyadic product of P0 w. P0
            CALL T2T2DYAD33(Pela,dddF,PDYADdddF)
            
            !---------------------------------------------------------
            ! Calculate DPDF
            
            ! Calculate the elastic part of DPDF
            DPelaDF = (lam+3.0d0/detJ*alpha*K*(theta-thetaZ)*(2.0d0-LOG(detJ)))*FITDYADFIT 
     &	+ (mu-lam*LOG(detJ)+3.0d0/detJ*alpha*K*(theta-thetaZ)*(1.0d0-LOG(detJ)))*FITUDYADFI 
     &  + mu*IODYADI
            
            DPDFPartial = fd*DPelaDF
             
            ! Calculate the inelastic DPDF
            DPDF = DPDFPartial + DfdDd*PDYADdddF
    

        END SUBROUTINE
        
        ! Tangent of the volumetric source w.r.t. the global damage var
        SUBROUTINE DPDT_NUM(F,Fn,phi,matpars,STATEV,STATEVN,theta,ddot,dt,dPdTheta)
          
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: F,Fn
            DOUBLE PRECISION, DIMENSION(13)      :: matpars
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV,STATEVN
            DOUBLE PRECISION                     :: phi,theta,ddot,dt
            
            ! Output variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: dPdTheta
            
            ! Internal variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: dP1,P,Finv
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV0
            DOUBLE PRECISION                     :: eps,PNEWDT,Y,rmech
            DOUBLE PRECISION                     :: ddot0
            DOUBLE PRECISION                     :: deltaTheta1                        
            
            eps = 10.0d0**(-8.0d0)
            
            CALL CALCINV33(F,Finv)
            
            STATEV0 = STATEVN
            ddot0 = ddot
            
            CALL KConstModel(F,Fn,Finv,phi,theta,matpars,STATEV0,PNEWDT,ddot0,dt,Y,rmech,P)
            
            dPdTheta = 0.0d0
            
            deltaTheta1 = theta + eps
            
            STATEV0 = STATEVN
            ddot0 = ddot
            
            CALL KConstModel(F,Fn,Finv,phi,deltaTheta1,matpars,STATEV0,PNEWDT,ddot0,dt,Y,rmech,dP1)

            dPdTheta = (dP1-P)/eps
            
            
        END SUBROUTINE
        
        ! Calculate the analytical tangent DPDT
        SUBROUTINE DPDT_ANA(F,Finv,matpars,DPDTheta,STATEV,phi,theta)
              
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3) :: F,Finv
            DOUBLE PRECISION, DIMENSION(13)  :: matpars
            DOUBLE PRECISION, DIMENSION(13)  :: STATEV
            DOUBLE PRECISION                 :: phi,theta
            
            ! Output variables
            DOUBLE PRECISION, DIMENSION(3,3) :: DPDTheta
            
            ! Internal variables
            DOUBLE PRECISION, DIMENSION(3,3) :: FinvT,Pela,DPelaDTheta
            DOUBLE PRECISION, DIMENSION(3,3) :: DPDThetaPartial,DPDd
            DOUBLE PRECISION                 :: mu,lam,c,alpha,K0,thetaZ
            DOUBLE PRECISION                 :: dens,betaD,etaD,dD,gamma
            DOUBLE PRECISION                 :: detJ,K,FDDF,DfdDd,fd
            DOUBLE PRECISION                 :: DDfdDDd,PsiEla,dddtheta
            DOUBLE PRECISION                 :: DPsiElaDtheta,d

            !---------------------------------------------------------
            ! read material parameters
        
            mu      = matpars(1)    
            lam     = matpars(2)
            c       = matpars(3)      
            alpha   = matpars(4) 
            K0      = matpars(5)    
            thetaZ  = matpars(6) 
            dens    = matpars(7)
            betaD   = matpars(8)
            etaD    = matpars(9)
            dD      = matpars(10)
            gamma   = matpars(11)

            !---------------------------------------------------------
            ! Set up some required quanaties
            d = STATEV(1)
            
            ! determinant of F
            CALL CALCDET33(F, detJ)
            
            ! transpose of inverse of F		
            FinvT = TRANSPOSE(Finv)
        
            ! bulk modulus
            K = lam + 2.0d0/3.0d0*mu
            
            ! Double contraction of F with F
            CALL T2DDT233(F,F,FDDF)
            
            ! Calculate damage function
            CALL DAMAGEFUNC(etaD,dD,d,fd,DfdDd,DDfdDDd)
            
            !---------------------------------------------------------
            ! Calculate some damage related qunataties
                       
            ! Calculate the elastic free energy
            PsiEla = mu/2.0d0*(FDDF - 3.0d0) - mu*LOG(detJ) 
     &               + lam/2.0d0*LOG(detJ)**(2.0d0)
     &               - 3.0d0*alpha*K*(theta-thetaZ)*LOG(detJ)/detJ
     &               + c*(theta-thetaZ-theta*LOG(theta/thetaZ))
     
     
            ! Calculate the derivative of free energy w.r.t. temp
            DPsiElaDtheta = c*(-LOG(theta/thetaZ))-3.0d0*alpha*K*LOG(detJ)/detJ
            
            IF (d == 0.0d0) THEN
                dddtheta = 0.0d0
            ELSE
                dddtheta = -DfdDd*DPsiElaDtheta/(DDfdDDd*PsiEla+betaD*gamma**(2.0d0))!-nd*DfdDd*dD*(1-fd)**(nd-1.0d0))
            END IF
            
            ! Calculate the elastic Piola stresses
            Pela = (lam*LOG(detJ)-mu)*FinvT + mu*F 
     &  - 3.0d0/detJ*alpha*K*(theta-thetaZ)*(1.0d0-LOG(detJ))*FinvT
            
            ! Calculate the elastic part of the tangent
            DPelaDTheta = -3.0d0/detJ*alpha*K*(1.0d0-LOG(detJ))*FinvT
            
            !---------------------------------------------------------
            ! Calculate DPDT
            DPDThetaPartial = fd*DPelaDTheta
            DPDd            = DfdDd*Pela
            
            ! Assemble the inelastic tangent
            DPDTheta = DPDThetaPartial + DPDd*dddtheta
            
        END SUBROUTINE

        ! Tangent of the volumetric source w.r.t. the global damage var
        SUBROUTINE DRDT_NUM(F,Fn,phi,matpars,STATEV,STATEVN,theta,ddot,dt,dRdTheta)
          
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: F,Fn
            DOUBLE PRECISION, DIMENSION(13)      :: matpars
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV,STATEVN
            DOUBLE PRECISION                     :: phi,theta,ddot,dt
            
            ! Output variables
            DOUBLE PRECISION                     :: dRdTheta
            
            ! Internal variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: Finv,P
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV0
            DOUBLE PRECISION                     :: eps,PNEWDT,Y
            DOUBLE PRECISION                     :: ddot0,drmech1,rmech
            DOUBLE PRECISION                     :: deltaTheta1                        
              
            eps = 10.0d0**(-8.0d0)
            
            CALL CALCINV33(F,Finv)
            
            STATEV0 = STATEVN
            ddot0 = ddot
            
            CALL KConstModel(F,Fn,Finv,phi,theta,matpars,STATEV0,PNEWDT,ddot0,dt,Y,rmech,P)
            
            dRdTheta = 0.0d0
            
            deltaTheta1 = theta + eps
            
            STATEV0 = STATEVN
            ddot0 = ddot
            
            CALL KConstModel(F,Fn,Finv,phi,deltaTheta1,matpars,STATEV0,PNEWDT,ddot0,dt,Y,drmech1,P)
            
            dRdTheta = (drmech1-rmech)/eps
            
            
        END SUBROUTINE
        
        ! Calculate the analytical tangent DrDtheta
        SUBROUTINE DRDT_ANA(F,Fn,theta,matpars,dt,dRdTheta,d,dn,phi)
          
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3) :: F,Fn
            DOUBLE PRECISION, DIMENSION(13)  :: matpars
            DOUBLE PRECISION                 :: theta,dt,d,dn,phi
              
            ! Output variables
            DOUBLE PRECISION                 :: dRdTheta
              
            ! Internal variables
            DOUBLE PRECISION                 :: mu,lam,c,alpha,K0,thetaZ
            DOUBLE PRECISION                 :: dens,betaD,etaD,dD,gamma
            DOUBLE PRECISION                 :: K,detJ,Jn,Jdot,ddot,FDDF
            DOUBLE PRECISION                 :: fd,PsiEla,DfdDd,DDfdDDd
            DOUBLE PRECISION                 :: DPsiElaDTheta,rElaGJ
            DOUBLE PRECISION                 :: rIH,D2PsiD2theta,DrGJDd
            DOUBLE PRECISION                 :: dddtheta,DrGJDtheta
            DOUBLE PRECISION                 :: DrIHDTheta,DrIHdd,nd

            !---------------------------------------------------------
            ! read material parameters
              
            mu      = matpars(1)    
            lam     = matpars(2)
            c       = matpars(3)      
            alpha   = matpars(4) 
            K0      = matpars(5)    
            thetaZ  = matpars(6) 
            dens    = matpars(7)
            betaD   = matpars(8)
            etaD    = matpars(9)
            dD      = matpars(10)
            gamma   = matpars(11)
            
            !---------------------------------------------------------
            ! Calculate some required quantaties
              
            ! determinant of F
            CALL CALCDET33(F, detJ)
              
            ! determinant of Fn	              
            CALL CALCDET33(Fn, Jn)
              
            ! rate of jacobian
            Jdot = (detJ - Jn)/dt
              
            ! rate of damage variable
            ddot = (d - dn)/dt
              
            ! bulk modulus
            K = lam + 2.0d0/3.0d0*mu
              
            ! Double contraction of F with F
            CALL T2DDT233(F,F,FDDF)
            
            ! Calculate damage function
            CALL DAMAGEFUNC(etaD,dD,d,fd,DfdDd,DDfdDDd)
            
            !-----------------------------------------------------------
            ! Calculate some damage related qunataties
            
            ! Calculate the elastic free energy
            PsiEla = mu/2.0d0*(FDDF - 3.0d0) - mu*LOG(detJ) 
     &               + lam/2.0d0*LOG(detJ)**(2.0d0)
     &               - 3.0d0*alpha*K*(theta-thetaZ)*LOG(detJ)/detJ
     &               + c*(theta-thetaZ-theta*LOG(theta/thetaZ))
     
              
            ! Calculate the derivative of free energy w.r.t. temp
            DPsiElaDTheta = c*(-LOG(theta/thetaZ))-3.0d0*alpha*K*LOG(detJ)/detJ
              
            ! c is wrong but works better
            D2PsiD2Theta = -c/theta
                 
            ! elastic heat source
            rElaGJ = -3.0d0*theta*alpha*K/(detJ**(2.0d0))*(1.0d0-LOG(detJ))*Jdot
        
            ! Inelastic heat
            rIH = (theta*DfdDd*DPsiElaDtheta 
     &                - DfdDd*PsiEla)*ddot
                    !+ betaD*gamma*(phi-gamma*d)
            
            nd = 2.0d0/3.0d0
            ! implicit derivative of d w.r.t theta
            IF (d == 0.0d0) THEN
                dddtheta = 0.0d0
            ELSE
                dddtheta = -DfdDd*DPsiElaDtheta/(DDfdDDd*PsiEla+betaD*gamma**(2.0d0))!-nd*DfdDd*dD*(1-fd)**(nd-1.0d0))
            END IF
            !---------------------------------------------------------
            !Calculate DRPDTheta     
     
            DrGJDtheta = fd*(-3.0d0/(detJ**(2.0d0))*alpha*K*(1.0d0-LOG(detJ))*Jdot)
            DrGJDd     = DfdDd*rElaGJ
            
            DrIHDtheta = (theta*DfdDd*D2PsiD2Theta - DfdDd*DPsiElaDtheta)*ddot
            DrIHDd = (theta*DDfdDDd*DPsiElaDtheta 
     &                - DDfdDDd*PsiEla)*ddot
     &                + rIH/dt
            

            dRdTheta = (DrGJDtheta + DrIHDtheta + (DrGJDd + DrIHDd)*dddtheta)*dens
        

        END SUBROUTINE

        ! Tangent of the volumetric source w.r.t. the global damage var
        SUBROUTINE DRDF_NUM(F,Fn,phi,matpars,STATEV,STATEVN,theta,ddot,dt,dRdF)
          
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: F,Fn
            DOUBLE PRECISION, DIMENSION(13)      :: matpars
            DOUBLE PRECISION, DIMENSION(13)      :: STATEV,STATEVN
            DOUBLE PRECISION                     :: phi,theta,ddot,dt
            
            ! Output variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: dRdF
            
            ! Internal variables
            DOUBLE PRECISION, DIMENSION(3,3)     :: deltaF1,P
            DOUBLE PRECISION, DIMENSION(3,3)     :: deltaFinv1,Finv
            DOUBLE PRECISION, DIMENSION(13)       :: STATEV0
            DOUBLE PRECISION                     :: eps,PNEWDT,Y
            DOUBLE PRECISION                     :: ddot0,drmech1,rmech
            INTEGER                              :: i,j
                         
              
            eps = 10.0d0**(-8.0d0)
            
            CALL CALCINV33(F,Finv)
            STATEV0 = STATEVN
            ddot0 = ddot
            CALL KConstModel(F,Fn,Finv,phi,theta,matpars,STATEV0,PNEWDT,ddot0,dt,Y,rmech,P)
            
            dRdF = 0.0d0
            Do i = 1,3
                DO j = 1,3
                    
                    deltaF1 = F
                    deltaF1(i,j) = deltaF1(i,j) + eps
                    
                    CALL CALCINV33(deltaF1,deltaFinv1)
                    
                    STATEV0 = STATEVN
                    ddot0 = ddot
                    
                    CALL KConstModel(deltaF1,Fn,deltaFinv1,phi,theta,matpars,STATEV0,PNEWDT,ddot0,dt,Y,drmech1,P)
                    
                    dRdF(i,j) = dRdF(i,j) + (drmech1-rmech)/(eps)
                
                END DO
              END DO
            
        END SUBROUTINE
        
        ! Calculate the analytical tangent DRDF
        SUBROUTINE DRDF_ANA(F,Finv,Fn,theta,matpars,dt,DRPLDF,STATEV,dn,phi)
        
            ! Input variables
            DOUBLE PRECISION, DIMENSION(3,3) :: F,Fn,Finv
            DOUBLE PRECISION, DIMENSION(13)  :: matpars
            DOUBLE PRECISION, DIMENSION(13)  :: STATEV
            DOUBLE PRECISION                 :: theta,dt,dn,phi
              
            ! Output variables
            DOUBLE PRECISION, DIMENSION(3,3) :: DRPLDF
              
            ! Internal variables
            DOUBLE PRECISION, DIMENSION(3,3) :: FinvT,Pela,DPelaDTheta
            DOUBLE PRECISION, DIMENSION(3,3) :: dddF,DrGJElaDF,DrGJDF
            DOUBLE PRECISION, DIMENSION(3,3) :: DrIHDF
            DOUBLE PRECISION                 :: mu,lam,c,alpha,K0,thetaZ
            DOUBLE PRECISION                 :: dens,betaD,etaD,dD,gamma
            DOUBLE PRECISION                 :: detJ,Jn,Jdot,ddot,K,FDDF
            DOUBLE PRECISION                 :: fd,PsiEla,DPsiElaDtheta
            DOUBLE PRECISION                 :: rElaGJ,DfdDd,DDfdDDd,rIH
            DOUBLE PRECISION                 :: DrGJDd,DrIHDd,d
            
            !---------------------------------------------------------
            ! read material parameters
              
            mu      = matpars(1)    
            lam     = matpars(2)
            c       = matpars(3)      
            alpha   = matpars(4) 
            K0      = matpars(5)    
            thetaZ  = matpars(6) 
            dens    = matpars(7)
            betaD   = matpars(8)
            etaD    = matpars(9)
            dD      = matpars(10)
            gamma   = matpars(11)
            
            !---------------------------------------------------------
            ! Set up some required quantaties
              
            ! determinant of F
            CALL CALCDET33(F, detJ)
              
            ! determinant of Fn	              
            CALL CALCDET33(Fn, Jn)
              
            ! rate of Jacobian
            Jdot = (detJ - Jn)/dt
              
            ! rate of d
            ddot = (d-dn)/dt
              
            ! calculate transpose of inverse		
            FinvT = TRANSPOSE(Finv)
              
            ! bulk modulus
            K = lam + 2.0d0/3.0d0*mu
              
            ! Double contraction of F with F
            CALL T2DDT233(F,F,FDDF)
            
            ! Calculate damage function
            CALL DAMAGEFUNC(etaD,dD,d,fd,DfdDd,DDfdDDd)

            !-----------------------------------------------------------
            ! Calculate some damage related qunataties
            
            ! Calculate the elastic free energy
            PsiEla = mu/2.0d0*(FDDF - 3.0d0) - mu*LOG(detJ) 
     &               + lam/2.0d0*LOG(detJ)**(2.0d0)
     &               - 3.0d0*alpha*K*(theta-thetaZ)*LOG(detJ)/detJ
     &               + c*(theta-thetaZ-theta*LOG(theta/thetaZ))
     
            ! Calculate the elastic Piola stresses
            Pela = (lam*LOG(detJ)-mu)*FinvT + mu*F 
     &  - 3.0d0/detJ*alpha*K*(theta-thetaZ)*(1.0d0-LOG(detJ))*FinvT
     
            ! Calculate the elastic part of the tangent
            DPelaDTheta = -3.0d0/detJ*alpha*K*(1.0d0-LOG(detJ))*FinvT
     
            ! Calculate the derivative of free energy w.r.t. temp
            DPsiElaDtheta = c*(-LOG(theta/thetaZ))-3.0d0*alpha*K*LOG(detJ)/detJ
     
            ! elastic heat source
            rElaGJ = -3.0d0*theta*alpha*K/(detJ**(2.0d0))*(1.0d0-LOG(detJ))*Jdot
              
            ! Calculate the derivative of d w.r.t. F
            dddF = -etaD*fd*Pela/(betaD*(phi-d)*etaD-betaD-etaD*fd)
            
            ! Inelastic heat
            rIH = (theta*DfdDd*DPsiElaDtheta 
     &                - DfdDd*PsiEla)*ddot
                     !+ betaD*gamma*(phi-gamma*d)
            
            DrGJElaDF = 3.0d0/detJ*alpha*K*theta*FinvT
     &  *((3.0d0-2.0d0*LOG(detJ))/detJ*Jdot-(1.0d0-LOG(detJ)/dt))
            !---------------------------------------------------------
            ! Calculate DRPLDF
            
            DrGJDF = fd*DrGJElaDF
            DrGJDd = DfdDd*rElaGJ
              
            DrIHDF = (theta*DfdDd*DPelaDTheta - DfdDd*Pela)*ddot
            DrIHDd = (theta*DDfdDDd*DPsiElaDtheta 
     &                - DDfdDDd*PsiEla)*ddot
     &                + rIH/dt
                      !- betaD*gamma**(2.0d0)
            
            DRPLDF = (DrGJDF + DrIHDF + (DrGJDd + DrIHDd)*dddF)*dens
            
        END SUBROUTINE

!-----------------------------------------------------------------------
!
!                   H E L P E R   R O U T I N E S
!
!-----------------------------------------------------------------------
        ! Compute the identity matrix
        SUBROUTINE IDENTITY(Ident)
              
              DOUBLE PRECISION, DIMENSION(3,3) :: Ident
              INTEGER :: i
	
              Ident = 0.0d0
              DO i = 1,3
                Ident(i,i) = 1.0d0
              END DO
	
        END SUBROUTINE
        
        ! Compute determinant of a 3x3 matrix
        SUBROUTINE CALCDET33(A, detA)
            
              DOUBLE PRECISION, DIMENSION(3,3) :: A
              DOUBLE PRECISION :: detA
            
              detA = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3))
     &          + A(1,2)*(A(3,1)*A(2,3) - A(2,1)*a(3,3))
     &          + A(1,3)*(A(2,1)*A(3,2) - A(3,1)*a(2,2))
        
        END SUBROUTINE

        ! Compute inverese of a 3x3 matrix
        SUBROUTINE CALCINV33(A, Ainv)
            
              DOUBLE PRECISION, dimension(3,3) :: A, Ainv
              DOUBLE PRECISION :: lwork(3)
              INTEGER :: info
              INTEGER :: n=3
              INTEGER :: ipiv(3)
              EXTERNAL DGETRI, DGETRF
		
              Ainv = A
		
              CALL DGETRF(n,n,Ainv,n,ipiv,info)
		
              If (info==0) Then
                CALL DGETRI(n,Ainv,n,ipiv,lwork,n,info)
              Else
                PRINT *, 'The matrix is singular'
              End If
        
        END SUBROUTINE
        
        SUBROUTINE FROBNORM33(A, Anorm)
        
              DOUBLE PRECISION, DIMENSION(3,3) :: A
              DOUBLE PRECISION                 :: Anorm
            
              INTEGER                          :: i,j
            
              Anorm = 0.0d0
              DO i = 1,3
                DO j = 1,3
                    Anorm = Anorm + A(i,j)**(2.0d0)
                END DO
              END DO
            
              Anorm = SQRT(Anorm)
        
        END SUBROUTINE
        
        SUBROUTINE FROBNORM66(A, Anorm)
        
              DOUBLE PRECISION, DIMENSION(6,6) :: A
              DOUBLE PRECISION                 :: Anorm
            
              INTEGER                          :: i,j
            
              Anorm = 0.0d0
              DO i = 1,6
                DO j = 1,6
                    Anorm = Anorm + A(i,j)**(2.0d0)
                END DO
              END DO
            
              Anorm = SQRT(Anorm)
        
        END SUBROUTINE
                
        ! Compute Voigt Notation of 3x3 tensor
        SUBROUTINE TENSOR33TOVOIGT(T, TVoi)
              
              DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: T
              DOUBLE PRECISION, DIMENSION(6),   INTENT(OUT) :: TVoi
              
              TVoi(1) = T(1,1)
              TVoi(2) = T(2,2)
              TVoi(3) = T(3,3)
              TVoi(4) = T(1,2)
              TVoi(5) = T(1,3)
              TVoi(6) = T(2,3)
        
        END SUBROUTINE
        
        ! Calculate the double contraction of two second order tensors
        SUBROUTINE T2DDT233(A,B,c)
            
              DOUBLE PRECISION, DIMENSION(3,3) :: A, B
              DOUBLE PRECISION :: c
              INTEGER :: i, j
            
              c = 0.0d0
            
              DO i = 1,3
                DO j = 1,3
                    c = c + A(i,j)*B(i,j)
                END DO
              END DO
            
        END SUBROUTINE
        
        ! Calculate the non standard upper bar dyadic product
        SUBROUTINE T2T2ODYAD33(A, B, C)
        
              DOUBLE PRECISION, DIMENSION(3,3) :: A, B
              DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C
              INTEGER :: i, j, k, l
              C = 0.0d0

              DO i = 1,3
                DO j = 1,3
                    DO k = 1,3
                        DO l = 1,3
                            C(i,j,k,l) = A(i,k)*B(j,l)
                        END DO
                    END DO
                END DO
              END DO
		
        END SUBROUTINE
        
        ! Calculate the standard dyadic product
        SUBROUTINE T2T2DYAD33(A, B, C)
            
              DOUBLE PRECISION, DIMENSION(3,3) :: A, B
              DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C
              INTEGER :: i, j, k, l
            
            C = 0.0d0

              DO i = 1,3
                DO j = 1,3
                    DO k = 1,3
                        DO l = 1,3
                            C(i,j,k,l) = A(i,j)*B(k,l)
                        END DO
                    END DO
                END DO
              END DO
		
        END SUBROUTINE

        ! Calculate the non standard lower bar dyadic product
        SUBROUTINE T2T2UDYAD33(A, B, C)
              
              DOUBLE PRECISION, DIMENSION(3,3) :: A, B
              DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C
              INTEGER :: i, j, k, l
            
              C = 0.0d0

              DO i = 1,3
                DO j = 1,3
                    DO k = 1,3
                        DO l = 1,3
                            C(i,j,k,l) = A(i,l)*B(j,k)
                        END DO
                    END DO
                END DO
              END DO
            
        END SUBROUTINE
        
        ! Calculate the standard dyadic product of two vectors
        SUBROUTINE VVDYAD3(a, b, C)
            
            DOUBLE PRECISION, DIMENSION(3) :: a, b
            DOUBLE PRECISION, DIMENSION(3,3) :: C
            INTEGER :: i,j
            
            C = 0.0d0
            
            DO i = 1,3
            DO j = 1,3
                C(i,j) = a(i)*b(j)
            END DO
            END DO
        END SUBROUTINE

        ! Perform the Push-Forward of the tangent DPDF to sigma related
        SUBROUTINE PUSHFORWARDDPDF(DPDF,F,sigma,DSIGDG)
            
              !Input
              DOUBLE PRECISION, DIMENSION(3,3,3,3)    :: DPDF
              DOUBLE PRECISION, DIMENSION(3,3)        :: F, sigma
            
              ! Output
              DOUBLE PRECISION, DIMENSION(3,3,3,3)    :: DSIGDG
            
              ! Internal variable
              DOUBLE PRECISION, DIMENSION(3,3)        :: ident, Ftrans
              DOUBLE PRECISION                        :: detJ
              INTEGER                                 :: i,j,k,l,m,n,o
            
            
              CALL IDENTITY(ident)
             
              CALL CALCDET33(F, detJ)
            
              Ftrans = TRANSPOSE(F)
            
              DSIGDG = 0.0d0
            
              DO i = 1,3
                DO j = 1,3
                    DO k = 1,3
                        DO l = 1,3
                            DO m = 1,3
                                DO n = 1,3
                                    DO o = 1,3
                                        DSIGDG(i,j,m,o) = DSIGDG(i,j,m,o) 
     &                 + ident(i,k)*F(j,l)*DPDF(k,l,m,n)*Ftrans(n,o)/detJ
                                    END DO
                                END DO
                            END DO
                        END DO
                    END DO
                END DO
              END DO 
            
              DO i = 1,3
                DO j = 1,3
                    DO m = 1,3
                        DO o = 1,3
                            DSIGDG(i,j,m,o) = DSIGDG(i,j,m,o) 
     &                                - ident(i,m)*sigma(j,o)
                        END DO
                    END DO
                END DO
              END DO
        
        END SUBROUTINE
                
        ! Perform the Pull-Back of the tangent DPDF to DSDC
        SUBROUTINE PULLBACKDPDF(DPDF,F,P, DSDC)
        
              ! Input
              DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: F, P
              DOUBLE PRECISION, DIMENSION(3,3,3,3) :: DPDF
            
              ! OUTPUT
              DOUBLE PRECISION, DIMENSION(3,3,3,3) :: DSDC
            
              ! Internally used variebles
              DOUBLE PRECISION, DIMENSION(3,3) :: C, Ident, Finv, FinvTrans
              DOUBLE PRECISION, DIMENSION(3,3) :: Ftrans, Cinv, S
              INTEGER :: i, j, k, l, q, r, m
            
              ! Set up some required tensors
              ! Identity matrix
              Ident = 0.0d0
              CALL IDENTITY(Ident) 
			
              ! Transpose of F
              Ftrans = TRANSPOSE(F)
            
              ! Inverse of F
              Finv = 0.0d0
              CALL CALCINV33(F, Finv)		
			
              ! Transpose of Inverse of F
              FinvTrans = TRANSPOSE(Finv)
            
              ! Left Cauchy green tensor
              C = MATMUL(Ftrans, F)
            
              ! Inverse of C
              Cinv = 0.0d0
              Call CALCINV33(C, Cinv) 
            
              ! Pull Back of P
              S = MATMUL(Finv, P)
            
              DSDC = 0.0d0
              ! Index notation of NFEM Script (4.36)
              DO i = 1,3
                DO j = 1,3
                    DO k = 1,3
                        DO l = 1,3
                            DO q = 1,3
                                DO r = 1,3
                                    DO m = 1,3
                                        DSDC(q,j,r,l) = DSDC(q,j,r,l) 
     &  + Finv(q,i)*DPDF(i,j,k,m)*FinvTrans(k,r)*Ident(m,l)
                                    END DO
                                END DO
                            END DO
                        END DO
                    END DO
                END DO
              END DO
            
              DO q = 1,3
                DO j = 1,3
                    DO r = 1,3
                        DO l = 1,3
                            DSDC(q,j,r,l) = DSDC(q,j,r,l)-Cinv(q,r)*S(j,l)
                        END DO
                    END DO
                END DO
              END DO
            
              DSDC = DSDC/2.0d0
        
        END SUBROUTINE
        
        ! Perform the Push forward of the DSDC tensor
        SUBROUTINE PUSHFORWARDDSDC(F,DSDC,DsigmaDg)
            
              DOUBLE PRECISION, DIMENSION(3,3,3,3) :: DSDC, DsigmaDg
              DOUBLE PRECISION, DIMENSION(3,3) :: F
              INTEGER :: i, j, k, l, m, n, o, p
            
              ! carry out the push-forward of DSDC (yielding DsigmaDg):
             DsigmaDg= 0.0d0
             DO i = 1,3
                DO j = 1,3
                    DO k = 1,3
                        DO l = 1,3
                            DO m = 1,3
                                DO n = 1,3
                                    DO o = 1,3
                                        DO p = 1,3
                            ! nutze "direkte indexnotation",
                            ! cf. Simo+Hughes, S.257, eq.(7.1.86)
                            DsigmaDg(i,j,k,l) = DsigmaDg(i,j,k,l)
     &                       + F(i,m)*F(j,n)*F(k,o)*F(l,p)*DSDC(m,n,o,p)
                                        END DO
                                    END DO
                                END DO
                            END DO
                        END DO
                    END DO
                END DO
              END DO
          
          END SUBROUTINE
          
        ! Perform the Jaumann Correction
        SUBROUTINE JaumannCorrection(sigma,DsigmaDg,DDSDDEtens)
            
              DOUBLE PRECISION, DIMENSION(3,3,3,3) :: DsigmaDg, DDSDDEtens
              DOUBLE PRECISION, DIMENSION(3,3) :: sigma, eye
              INTEGER :: i, j, k, l
            
              ! add Jaumann rate correction to DsigmaDg for the final
              ! Abaqus tangent DDSDDE (as 4th order tensor):
              ! required here: Cauchy stress (3.42)
            
              eye = 0.0d0
              forall(i = 1:3) eye(i,i) = 1.0d0
            
              DO i = 1,3
                DO j = 1,3
                    DO k = 1,3
                        DO l = 1,3
                            DDSDDEtens(i,j,k,l) = DsigmaDg(i,j,k,l)
     &                        + eye(i,k)*sigma(j,l)
     &                        + sigma(i,k)*eye(j,l)
                        END DO
                    END DO
                END DO
              END DO
            
        END SUBROUTINE
        
        
        ! Perform Voigt transformation
        SUBROUTINE T3333TOVOIGT(ABATAN, DDSDDE)
        
              ! transforms a 4th order abaqus tangent tensor (push-forward of
              ! dC/dS with added 'corotational correction tensor') to the
              ! final 6x6 Voigt-type DDSDDE return format required by Abaqus
              
              DOUBLE PRECISION, DIMENSION(3,3,3,3), INTENT(IN)  :: ABATAN
              DOUBLE PRECISION, DIMENSION(6,6),     INTENT(OUT) :: DDSDDE
              
              ! transfer to Abaqus Voigt form
              ! 11 kl
              DDSDDE(1,1) = ABATAN(1,1,1,1)
              DDSDDE(1,2) = ABATAN(1,1,2,2)
              DDSDDE(1,3) = ABATAN(1,1,3,3)
              DDSDDE(1,4) = ABATAN(1,1,1,2)
              DDSDDE(1,5) = ABATAN(1,1,1,3)
              DDSDDE(1,6) = ABATAN(1,1,2,3)
              ! 22 kl
              DDSDDE(2,1) = ABATAN(2,2,1,1)
              DDSDDE(2,2) = ABATAN(2,2,2,2)
              DDSDDE(2,3) = ABATAN(2,2,3,3)
              DDSDDE(2,4) = ABATAN(2,2,1,2)
              DDSDDE(2,5) = ABATAN(2,2,1,3)
              DDSDDE(2,6) = ABATAN(2,2,2,3)
              ! 33 kl
              DDSDDE(3,1) = ABATAN(3,3,1,1)
              DDSDDE(3,2) = ABATAN(3,3,2,2)
              DDSDDE(3,3) = ABATAN(3,3,3,3)
              DDSDDE(3,4) = ABATAN(3,3,1,2)
              DDSDDE(3,5) = ABATAN(3,3,1,3)
              DDSDDE(3,6) = ABATAN(3,3,2,3)
              ! 12 kl
              DDSDDE(4,1) = ABATAN(1,2,1,1)
              DDSDDE(4,2) = ABATAN(1,2,2,2)
              DDSDDE(4,3) = ABATAN(1,2,3,3)
              DDSDDE(4,4) = ABATAN(1,2,1,2)
              DDSDDE(4,5) = ABATAN(1,2,1,3)
              DDSDDE(4,6) = ABATAN(1,2,2,3)
              ! 13 kl
              DDSDDE(5,1) = ABATAN(1,3,1,1)
              DDSDDE(5,2) = ABATAN(1,3,2,2)
              DDSDDE(5,3) = ABATAN(1,3,3,3)
              DDSDDE(5,4) = ABATAN(1,3,1,2)
              DDSDDE(5,5) = ABATAN(1,3,1,3)
              DDSDDE(5,6) = ABATAN(1,3,2,3)
              ! 23 kl
              DDSDDE(6,1) = ABATAN(2,3,1,1)
              DDSDDE(6,2) = ABATAN(2,3,2,2)
              DDSDDE(6,3) = ABATAN(2,3,3,3)
              DDSDDE(6,4) = ABATAN(2,3,1,2)
              DDSDDE(6,5) = ABATAN(2,3,1,3)
              DDSDDE(6,6) = ABATAN(2,3,2,3)
        
        END SUBROUTINE

        ! Print 3x3 matrix
        SUBROUTINE Kprint33(str, mat)
        
              CHARACTER (LEN=6) :: str
              DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: mat
            
              WRITE(*,*) str, ':'
              WRITE(*,*) mat(1,1), mat(1,2), mat(1,3)
              WRITE(*,*) mat(2,1), mat(2,2), mat(2,3)
              WRITE(*,*) mat(3,1), mat(3,2), mat(3,3)
        
        END SUBROUTINE Kprint33
        
        ! Print 6x6 matrix
        SUBROUTINE Kprint66(str, mat)
        
              CHARACTER (LEN=6) :: str
              DOUBLE PRECISION, DIMENSION(6,6), INTENT(IN)  :: mat
            
              WRITE(*,*) '--'
              WRITE(*,*) str, ':', mat(1,1), mat(1,2), mat(1,3), mat(1,4), mat(1,5), mat(1,6)
              WRITE(*,*) str, ':', mat(2,1), mat(2,2), mat(2,3), mat(2,4), mat(2,5), mat(2,6)
              WRITE(*,*) str, ':', mat(3,1), mat(3,2), mat(3,3), mat(3,4), mat(3,5), mat(3,6)
              WRITE(*,*) str, ':', mat(4,1), mat(4,2), mat(4,3), mat(4,4), mat(4,5), mat(4,6)
              WRITE(*,*) str, ':', mat(5,1), mat(5,2), mat(5,3), mat(5,4), mat(5,5), mat(5,6)
              WRITE(*,*) str, ':', mat(6,1), mat(6,2), mat(6,3), mat(6,4), mat(6,5), mat(6,6)
        
        END SUBROUTINE Kprint66


