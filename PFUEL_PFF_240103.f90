!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ABAQUS user defined phase field element subroutine
!
! Stephan Roth, TU Bergakademie Freiberg, 30.07.2020
!
! 30.07.2020: Multi-phase multi-component
! 15.06.2021: concentrations independent on phases
! 01.11.2023: without NCP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                   mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                   ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                   pnewdt,jprops,njprop,period)

      USE ABQINTERFACE
      USE ABAModul
      
      !!!!!!!!!!
      ! Debugging
      !!!!!!!!!!
      USE BulkEnergyModule
      USE CrackSurfaceEnergyModule
	  !!!!!!!!!!
	  ! 
	  USE KineticEnergyModule ! only for displaying the energy, no impact on the UMAT

      IMPLICIT NONE
      

      ! UEL variables
      INTEGER(kind=AbqIK), INTENT(IN) :: nprops, njprop, mcrd, nnode, ndofel, &
                                         npredf, mlvarx, nrhs, nsvars, jtype, &
                                         kstep, kinc, jelem, ndload, mdload, &
                                         lflags(*), jdltyp(mdload,*), &
                                         jprops(njprop)
      REAL(kind=AbqRK), INTENT(IN) :: props(nprops), &
                                      coords(mcrd,nnode), u(ndofel), period, &
                                      v(ndofel), a(ndofel), time(2), dtime, &
                                      predef(2,npredf,nnode), du(mlvarx,*), &
                                      adlmag(mdload,*), ddlmag(mdload,*), params(*)
      REAL(kind=AbqRK), INTENT(INOUT) :: energy(8), pnewdt, svars(nsvars)
      REAL(kind=AbqRK), INTENT(OUT) :: amatrx(ndofel,ndofel), rhs(mlvarx,*)

      ! UMAT variables
      INTEGER(kind=AbqIK) :: ntens, layer, kspt, npt
      REAL(kind=AbqRK) :: predef_umat(npredf), dpred(npredf), rpl, drpldt, celent, F0(3,3), F1(3,3), GPcoords(3), drot(3,3)
      REAL(kind=AbqRK), ALLOCATABLE :: ddsddt(:), drplde(:)
      CHARACTER*45 :: cmname

      ! further variables
      REAL(kind=AbqRK) :: fint(ndofel), energy_gp(8), JacobiDet, prop_thickness, prop_debug
      REAL(kind=AbqRK), ALLOCATABLE :: Matrix_B(:,:), stran(:), dstran(:), Ct(:,:), stress(:)
      INTEGER(kind=AbqIK) :: i1, i2, i3
      INTEGER(kind=AbqIK) :: nNodalDOF, nPp
      LOGICAL :: isNlgeom, positivDetJ, k_RHS, k_K, k_M, positivDetJTemp, is_print_inc, is_print_elem, debug_Flag
      
      ! Dynamic STEPS
      LOGICAL :: isDynamic
      REAL(kind=AbqRK) :: kinED_gp
      REAL(kind=AbqRK) :: prop_dynamic, prop_density
      REAL(kind=AbqRK), ALLOCATABLE :: MassM(:,:), Matrix_N(:,:), vel_gp(:) ! DYNAMIC STEPS
      
      ! pure debug Variables
      INTEGER(kind=AbqIK), PARAMETER :: n_debug_elem = 1
      INTEGER(kind=AbqIK), PARAMETER :: n_debug_inc = 10
      INTEGER(kind=AbqIK), PARAMETER :: n_debug_ip = 1
      
      INTEGER(kind=AbqIK) :: debug_elements(n_debug_elem)
      INTEGER(kind=AbqIK) :: debug_increments(n_debug_inc)
      INTEGER(kind=AbqIK) :: debug_ipoints(n_debug_ip)
      
      DATA debug_elements /1/
      DATA debug_increments /1,2,3,4,5,500,501,502,503,504/
      DATA debug_ipoints /1/



!~ 		!!! Martinez-Paneda Vergleich print !!!
		! Zeilen für u11, u12, u21, u22, u31, u32, u41, u42
		! Das sind die alten Zeilen: 1, 2, 4, 5, 7, 8, 10, 11
		INTEGER, DIMENSION(8) :: row_u = (/1, 2, 4, 5, 7, 8, 10, 11/)
		INTEGER, DIMENSION(8) :: col_u = (/1, 2, 4, 5, 7, 8, 10, 11/)
		! Zeilen für phi1, phi2, phi3, phi4
		! Das sind die alten Zeilen: 3, 6, 9, 12
		INTEGER, DIMENSION(4) :: row_phi = (/3, 6, 9, 12/)
		INTEGER, DIMENSION(4) :: col_phi = (/3, 6, 9, 12/)
      
      ! initialisation
      positivDetJ=.TRUE.; isNlgeom=.FALSE.; k_RHS=.TRUE.; k_K=.TRUE.; k_M=.FALSE.; positivDetJTemp=.TRUE.
      
      debug_Flag=.FALSE.
	  ! Flags for Debug
	  is_print_inc = .FALSE.
	  is_print_elem = .FALSE.
	  
      ! umat dummies
      layer = 1; kspt = 1
      celent = one; F0 = zero; F1 = zero
      
      

      ! 'Wer bin ich, was kann ich ?'

      ! evaluation of lflags
      ! geometric nonlinear
      IF (lflags(2) .EQ. 1) isNlgeom = .TRUE.
      
      ! props:
      ! 
      ! 1:   0.
      ! 2:   nCpFloat
      ! 3:   0.
      ! 4:   Flag numerical Tangent (1 - yes, 0 - analytical)
      ! 5:   E-Modulus
      ! 6:   poissons ratio
      ! 7:   gamma star (star convex split, de Lorenzis)
      ! 8:   l0 (internal length)
      ! 9:   Gc (fracture toughness)
      ! 10:  Flag Solver (1 - monolithic, 2 - BFGS, 3 - staggered)
      ! 11:  Flag Coupling terms (1 - use coupling terms, 2 - no coupling terms)
      ! 12:  thickness (plain strain)
      ! 13:  Flag print (1 - print, 2 - no print)
      ! 14:  Flag check stress
      ! 15:  Flag check tangent
      ! 16:  viscous parameter Phase-Field (eta, viscous regularization, Miehe)
      ! 17:  GcI
      ! 18:  GcII
      ! 19:  Flag dynamic step (1 - dynamic step, 0 - no dynamic terms)
      ! 20:  density (for dynamic steps)
      !
      
      
      ! jprops
	  !
	  ! jprop1 = reduced Integration
	  ! integration: 1 - reduced, 0 - full (default), 2 - full (2nd integration scheme)
	  ! jprop2 = axi-/spheri-symmetric
	  ! jprop3 = kLinearity?
	  ! jprop4 = number properties (nprop)
	  ! jprop5 = dimension
	  ! jprop6 = 0
	  ! NDI, NSHR is determined in UELlib
	  ! 
	  IF (time(2) .EQ. zero .AND. jelem .EQ. one) THEN
		WRITE(6,*) '========================================'
		WRITE(6,*) ' Check props array PF:'
		WRITE(6,*) '========================================'
		WRITE(6,*) ' 1  dummy                 = ', props(1)
		WRITE(6,*) ' 2  nCpFloat             = ', props(2)
		WRITE(6,*) ' 3  dummy                = ', props(3)
		WRITE(6,*) ' 4  Flag numerical Tangent= ', props(4)
		WRITE(6,*) ' 5  E-Modulus            = ', props(5)
		WRITE(6,*) ' 6  Poisson ratio        = ', props(6)
		WRITE(6,*) ' 7  gamma star           = ', props(7)
		WRITE(6,*) ' 8  l0                   = ', props(8)
		WRITE(6,*) ' 9  Gc                   = ', props(9)
		WRITE(6,*) '10  Flag Solver          = ', props(10)
		WRITE(6,*) '11  Flag Coupling terms  = ', props(11)
		WRITE(6,*) '12  thickness            = ', props(12)
		WRITE(6,*) '13  Flag print           = ', props(13)
		WRITE(6,*) '14  Flag check stress    = ', props(14)
		WRITE(6,*) '15  Flag check tangent   = ', props(15)
		WRITE(6,*) '16  eta (viscous PF)     = ', props(16)
		WRITE(6,*) '17  GcI                  = ', props(17)
		WRITE(6,*) '18  GcII                 = ', props(18)
		WRITE(6,*) '19  Flag dynamic step    = ', props(19)
		WRITE(6,*) '20  density              = ', props(20)
		write(6,*) '========================================'
	  END IF
	  
	  
		
      SELECT CASE(lflags(3))
      CASE(1)
        k_RHS = .TRUE.
        k_K   = .TRUE.
        k_M   = .FALSE.
      CASE(2)
        k_RHS = .FALSE.
        k_K   = .TRUE.
        k_M   = .FALSE.
      CASE(3)
        k_RHS = .FALSE.
        k_K   = .FALSE.
        k_M   = .FALSE.
      CASE(4)
        k_RHS = .FALSE.
        k_K   = .FALSE.
        k_M   = .TRUE.
      CASE(5)
        k_RHS = .TRUE.
        k_K   = .FALSE.
        k_M   = .FALSE.
      CASE(6)
        k_RHS = .FALSE.
        k_K   = .FALSE.
        k_M   = .TRUE.
      CASE DEFAULT
        k_RHS = .FALSE.
        k_K   = .FALSE.
        k_M   = .FALSE.
        write(7,*) 'no output required'
      END SELECT
      
      ! DYNAMIC Parameters
      ! Initialize
      prop_density = 0.
      
      
      prop_dynamic = props(19)
      ! Dynamic Process?
      IF (prop_dynamic .GT. zero) THEN
		isDynamic = .TRUE.
		prop_density = props(20)
	  END IF
      
      
      
      ! Important to Make sure the Calculation of the Matrix_N is also conducted for Steps with k_M
      ! (see in Materialroutine: CALL SORT_SHAPEFUNC(...))
      IF (k_M) THEN
		IF (.NOT. isDynamic .AND. prop_density .GT. zero) THEN
		   WRITE(6,*) ''
		   WRITE(6,*) '*****************************************************'
		   WRITE(6,*) ' ACHTUNG: Statischer Step mit definierter Dichte!'
		   WRITE(6,*) ' Es wurde kein dynamischer Step gewaehlt.'
		   WRITE(6,*) ''
		   WRITE(6,*) '   isDynamic   = ', isDynamic
		   WRITE(6,*) '   prop_density= ', prop_density
		   WRITE(6,*) ''
		   WRITE(6,*) ' Die Masse wird daher ignoriert.'
		   WRITE(6,*) '*****************************************************'
		   WRITE(6,*) ''
		END IF
	  END IF
      
      ! number of nodal DOF
      nNodalDOF = ndofel/nnode
      nPp       = 1 ! number of order parameters, here: damage variable 

      ! number of tensor coordinates (generalised): strain, damage variable
      ntens=NDI+NSHR+nPp*(1+D)

      ! allocation of UMAT-matrices
      ALLOCATE(ddsddt(ntens), drplde(ntens), Matrix_B(ntens,ndofel), Ct(ntens,ntens), stran(ntens), dstran(ntens), stress(ntens))
      ddsddt=zero; drplde=zero; Matrix_B = zero; stran = zero; dstran = zero; Ct = zero; stress = zero
      
      ! ALLOCATION DYNAMIC STEP VARIABLES
      ALLOCATE(MassM(ndofel, ndofel))
      
	  ALLOCATE(Matrix_N(nNodalDOF, ndofel))
	  ALLOCATE(vel_gp(D))
	  vel_gp = zero
	  MassM = zero
	  Matrix_N = zero
      
      ! query at time 0: number of internal state variables per integration point and material parameters: answer from UMAT
      IF ( time(2) .EQ. zero ) THEN
        IF (nsvars .LT. NGP*numSDV) THEN
          WRITE(7,*) "A bigger number of solution-dependent state variables is required.", NGP*numSDV, nsvars
          WRITE(7,*) "Please note that this is a total number of SDVs, not per integration point."
          CALL XEXIT()
        END IF
        IF (nprops .LT. numMatPar) THEN
          WRITE(7,*) "A bigger number of real parameters is required.", numMatPar, nprops
          CALL XEXIT()
        END IF
        ! check of material parameters: answer from UMAT
        CALL CheckMaterialParameters(props)
      END IF

      ! thickness of 2D-Elements
      prop_thickness = one
      IF (D .EQ. 2) prop_thickness = props(thicknessIndex)

      ! Ende 'Wer bin ich, was kann ich ?'



      ! UMAT call only when neccessary
!      IF (k_RHS .OR. k_K) THEN
        ! loop over integration points
        energy = zero
        amatrx = zero

        fint   = zero
        !
        ! ANFANG PRINT
        !
		prop_debug = props(13)
		IF (prop_debug .GT. zero) THEN
			debug_Flag = .TRUE.
		END IF
		  
		  IF (debug_Flag) THEN
			  ! See which printouts
			  DO i1=1,n_debug_elem
			   IF (jelem.EQ.debug_elements(i1)) THEN
				is_print_elem = .true.
				EXIT
			   ENDIF
			  END DO
			  
			  DO i1=1,n_debug_inc
			   IF (kinc.EQ.debug_increments(i1)) THEN
				is_print_inc = .true.
				EXIT
			   ENDIF
			  END DO
		  END IF
		  
        DO npt=1,NGP
          ! compute JACOBI determinante and B-matrix
          CALL BMatrixJacobian(coords(1:D,1:nnode),u,D,nnode,ndofel,NDI,NSHR,ntens,njprop,ShapeFunc(GPPos(npt,1:D)), &
                               ShapeFuncDeriv(GPPos(npt,1:D)),jprops,prop_thickness,isNlgeom,Matrix_B,JacobiDet,drot,GPcoords,positivDetJ)

          IF (positivDetJ) THEN
          ! element OK
            ! compute local (transformed) separation
            stran  = matmul(Matrix_B, u)
            dstran = matmul(Matrix_B, du(1:ndofel,1))
            ! here extract stress from svars and compute F0,F1...
            !=======================================================
            stress = zero
            !stress = svars(2:1+ndi+nshr)
            energy_gp = zero
            vel_gp = zero
            !=======================================================
            ! evaluation constitutive law
			!
            !
            CALL umat(stress,svars(numSDV*(npt-1)+1:numSDV*npt),Ct,energy_gp(2), &
                      energy_gp(4),energy_gp(3),rpl,ddsddt,drplde,drpldt,stran, &
                      dstran,time,dtime,predef_umat(1),dpred(1),predef_umat,dpred,cmname,NDI,NSHR,ntens, &
                      numSDV,props,nprops,GPcoords,drot,pnewdt,celent,F0,F1,jelem,npt,layer,kspt,kstep,kinc)
			!
			! kinetic Energy, Velocity at GP
			!
			IF (isDynamic) THEN
				CALL SORT_SHAPEFUNC(Matrix_N, ShapeFunc(GPPos(npt,1:D)), nnode, nPp, D) ! Sicherstellen dass für k_M auch isDynamic sicher aktiviert ist!
				vel_gp = MATMUL(Matrix_N,v)
				kinED_gp = kinED(prop_density, vel_gp, D)
			END IF
			!
			!
			IF (k_M) THEN
		    ! mass matrix
		    ! Matrix_N is the Shapefunction-Matrix (explained in SORT_SHAPEFUNC) with evaluated N(npt) at GP npt
			  MassM = MassM + prop_density * GPWeight(npt)*JacobiDet * matmul(transpose(Matrix_N),Matrix_N)
!~ 			  MassM_debug = MassM_debug + prop_density * GPWeight(npt)*JacobiDet
			END IF
			! 
			! Energy Densities
			!
            energy_gp(1) = svars(numSDV*(npt-1)+17) ! (SDV 17) -> crack bulk Energy
            energy_gp(2) = svars(numSDV*(npt-1)+16) ! (SDV 16) -> stored bulk Energy
            energy_gp(3) = svars(numSDV*(npt-1)+17) ! (SDV 17) -> crack bulk Energy
            energy_gp(4) = svars(numSDV*(npt-1)+18) ! (SDV 18) -> viscous dissipation Energy
            energy_gp(5) = kinED_gp ! -> kinetic Energy
            energy_gp(6) = svars(numSDV*(npt-1)+19) ! (SDV 19) -> total Energy Bulk
            energy_gp(7) = 0. ! 

            ! integration over element
            ! energies
            !
            ! 1 = ALLKE -> reine Bruchenergie aus Bulk Bereich (s. PFUEL_PFF)
            ! 2 = ALLSE -> PFFCZ kombinierte linear elast. stored Energy (gesamt elast. gesp. Energie)
            ! 3 = ALLCD -> PFFCZ kombinierte gesp. Bruchenergie (gesamt Bruchenergie)
            ! 4 = ALLPD -> Kontaktenergie CZ bei negativer Normalseparation + viskose Regularisierung aus PFF (keine physikalische Interpretation. Nur um zu checken ob Energieterme auftreten)
            ! 5 = ALLVD -> kinetische Energie PF
            ! 6 = ALLAE -> totale Energie CZ + totale Energie PFF (gesamte Energie im System)
            ! 7 = ALLEE -> reine Bruchenergie aus CZ Bereich (s. CZUEL.f90)

            ! integration over element
            ! energies
            energy(:) = energy(:) + GPWeight(npt)*JacobiDet * energy_gp(:)
            ! stiffness matrix
            amatrx = amatrx + GPWeight(npt)*JacobiDet * matmul(transpose(Matrix_B),matmul(Ct,Matrix_B))
            
            ! internal force vector
            fint = fint - GPWeight(npt)*JacobiDet * matmul(transpose(Matrix_B),stress)
			
          ELSE
          ! distorted element
            ! stop loop
            EXIT
          END IF
        END DO



!~ 	   IF (is_print_elem .AND. is_print_inc .AND. npt .EQ. NGP) THEN
!~ 	   !
!~ 			WRITE(6,*) ' === FINALE ELEMENT-MATRIZEN ==='
!~ 			WRITE(6,*) ' Steifigkeitsmatrix amatrx (mechanisch/elastisch):'

!~ 			DO i1 = 1, 8
!~ 			   WRITE(6,'(A,I2,A,8(ES12.4,1X))') '  Zeile', i1, ':', &
!~ 				  (amatrx(row_u(i1), col_u(i2)), i2=1,8)
!~ 			END DO

!~ 			WRITE(6,*) ''
!~ 			WRITE(6,*) ' Steifigkeitsmatrix amatrx (Phase field):'

!~ 			DO i1 = 1, 4
!~ 			   WRITE(6,'(A,I2,A,4(ES12.4,1X))') '  Zeile', i1+8, ':', &
!~ 				  (amatrx(row_phi(i1), col_phi(i2)), i2=1,4)
!~ 			END DO

!~ 			WRITE(6,*) ''
!~ 		  WRITE(6,*) ' Residuum rhs (mechanisch/elastisch):'
!~ 		  ! rhs 1, 2, 4, 5, 7, 8, 10, 11
!~ 		  DO i1 = 1, 11, 3
!~ 			 WRITE(6,'(A,I2,A,ES12.4)') '  rhs(', i1, ') =', fint(i1)
!~ 			 WRITE(6,'(A,I2,A,ES12.4)') '  rhs(', i1+1, ') =', fint(i1+1)
!~ 		  END DO
!~ 		  WRITE(6,*) ''
!~ 		  WRITE(6,*) ' Residuum rhs (Phase field):'
!~ 		  ! rhs 3, 6, 9, 12
!~ 		  DO i1 = 3, 12, 3
!~ 			 WRITE(6,'(A,I2,A,ES12.4)') '  rhs(', i1, ') =', fint(i1)
!~ 		  END DO
!~ 		  WRITE(6,*) ''
!~ 		  WRITE(6,*) '-----------------------------------------------'
!~ 	   END IF




        ! zero matrices for distorted elements
        IF (.NOT. positivDetJ) THEN
          WRITE(7,*) 'non-positive JACOBIAN det(J)<=0, distorted Element: ', jelem, ' matrices set to zero'
          amatrx = zero; fint = zero; energy = zero
        END IF

        rhs(1:ndofel,1) = fint


      IF (k_M) THEN
      ! mass matrix
        amatrx = MassM
      ELSE IF (.NOT. k_K) THEN
      ! no matrix required
        amatrx = zero
      END IF

      IF (.NOT. k_RHS) THEN
      ! no rhs required
        rhs(1:ndofel,1) = zero
      END IF
      

      DEALLOCATE(Matrix_B, ddsddt, drplde, stran, dstran, Ct, stress)
      DEALLOCATE(MassM, Matrix_N, vel_gp) ! DYNAMIC STEPS

      
      RETURN

!      CONTAINS

!------------------------------------------------------------------------------------

    END SUBROUTINE UEL

!------------------------------------------------------------------------------------
SUBROUTINE SORT_SHAPEFUNC(N_sorted, A, nnode, nPp, D)

! ======================================================================
!     SORT_SHAPEFUNC - Create Matrix with Shapefunctions evaluated at npt (a single GP)
! ======================================================================
!     Creates Shape Function Matrix with 
!     Input:
!       A(nnode)   - Array with Values Shape_Function N_Matrix = [N1, N2, ..., Nnnode]
!       nnode      - number of Shapefunctions (Number of Nodes)
!       nPp        - Add. DOF per Node (phasen)
!       D          - mech. DOF per Node
!     Output:
!       N_Matrix(DOF, DOF*nnode) - sorted N-Matrix (DOF = D+nPp)
!     Example 	 D=2, nPp=1, nnode = 4:
!     N_Matrix = [ N1  0  0  N2  0  0  N3  0  0  N4  0  0 ]
!         		 [  0 N1  0   0 N2  0   0 N3  0   0 N4  0 ]
!         		 [  0  0  0   0  0  0   0  0  0   0  0  0 ]
!     Example 	 D=3, nPp=0, nnode = 4:
!     N_Matrix = [ N1  0  0  N2  0  0  N3  0  0  N4  0  0 ]
!         		 [  0 N1  0   0 N2  0   0 N3  0   0 N4  0 ]
!         		 [  0  0 N1   0  0 N2   0  0 N3   0  0 N4 ]
! ======================================================================

  IMPLICIT NONE
  INTEGER(kind=AbqIK), INTENT(IN) :: nnode, nPp, D
  
  REAL(kind=AbqRK), INTENT(IN) :: A(nnode)
  INTEGER(kind=AbqIK) :: DOF
  REAL(kind=AbqRK), INTENT(OUT) :: N_sorted(D+nPp, (D+nPp)*nnode)
  
  INTEGER :: i1, i2, idx_node, row_idx, col_offset
  
  ! DOF per node
  DOF = D + nPp
  
  ! Init
  DO i1 = 1, DOF
    DO i2 = 1, DOF*nnode
      N_sorted(i1, i2) = zero
    END DO
  END DO
  
  ! Values insert
  DO idx_node = 1, nnode
    DO row_idx = 1, D
      col_offset = (idx_node - 1) * DOF + row_idx
      N_sorted(row_idx, col_offset) = A(idx_node)
    END DO
  END DO
  
  RETURN
END SUBROUTINE SORT_SHAPEFUNC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

