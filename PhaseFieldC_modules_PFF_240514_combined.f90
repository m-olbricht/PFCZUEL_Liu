!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! UMAT material subroutine of multi-phase multi-component phase field model
!
! Stephan Roth, TU Bergakademie Freiberg, 30.07.2020
!
! 30.07.2020: Multi-phase multi-component
! 23.11.2020: all derivatives tested (totalPotential --> stresses --> tangent)
! 15.06.2021: concentrations independent on phases
! 19.07.2021: phase-dependent parameters
! 04.03.2022: chemical reactions
! 03.01.2024: without NCP, DOF: displacements, 1xdamage variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Martin Olbricht 14.05.2024
!
! 26.07.2024: implementation of Amor split
!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SAVED VARIABLES AND INPUT PARAMETERS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SVR: svr(1)     :

!   1:  phase Parameter
!   2-7: stress tensor 
!   8-13: Verzerrungstensor 
!   16: bulkEnergy
!   17: Crack Surface Energy
!   18: Viscous Dissipation Energy
!   19: total Energy
!   21-23: gradient phaseparameter
!   25-27: stiffness in principal direction
!   29: -
!   30: H_sph
!   31: H_dev
!   32: 
!   33: 
!   34: 
!   35: 
!   36: 
!   37: 
!   38:
!   39: 
!   40:  
!
! material parameters:  
!
!   1:  ?
!   2:  number Phases
!   3:  ?
!   4:  numerical Tangent
!   5:  Young's modulus
!   6:  Poisson's ratio
!   7:  gamma star
!   8:  internal length
!   9:  fracture toughness ! USELESS FOR MIXEDMODE-FORMULATION
!   10: Solver (1 -> monolithic, 2 -> Quasi-Newton "Monolithic BFGS" 3 -> staggered "one way partitioned")
!	11: coupling terms (0. -> keine Koppelterme, 1. -> Koppelterme)
!	12: thickness
!   13: Write Flag (0. -> no Prints; 1. -> Print info Elements)
!	14: check Stress (0. -> no stress; 1. -> compare numerical / analytical stress)
!	15: check Tangent (0. -> no check; 1. -> compare numerical / analytical tangent)
!   16: numerical crackresistence eta
!	17: GcI 
!	18: GcII
!	19: dynamic Flag
!	20: density
!
!
! integer properties:
!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PhaseField_module

  USE ABQINTERFACE_PF
  USE SharedValues
  

  IMPLICIT NONE

  PUBLIC :: CheckMaterialParameters, umatPF
  
  
  CONTAINS

!------------------------------------------------------------------------------------

    SUBROUTINE CheckMaterialParameters(props)
      ! Check of all material parameters

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE ABAModul
      USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: props(numMatPar)
      REAL(kind=AbqRK) :: prop_E, prop_nu, prop_gstar, prop_l0, prop_Gc

      ! read material parameters from props
      prop_E      = props( 5) ! Young's modulus
      prop_nu     = props( 6) ! Poisson's ratio
      prop_gstar  = props( 7) ! gamma star, Parameter Star Split Alg
      prop_l0     = props( 8) ! internal length
      prop_Gc     = props( 9) ! fracture toughness

      ! parameter check
      IF (prop_E .LE. zero) THEN
        WRITE(7,*) "Young's modulus should exceed zero. EXIT", prop_E
        CALL XEXIT()
      END IF
      IF (prop_nu .LT. zero) THEN
        WRITE(7,*) "Poisson's ratio must be positive. EXIT", prop_nu
        CALL XEXIT()
      END IF
      IF (prop_gstar .LT. -one) THEN
        WRITE(7,*) "Star Split Parameter must be > -1, EXIT", prop_gstar
        CALL XEXIT()
      END IF
      !---------------------------------------!
      ! INSERT PARAMETER CHECK FOR l0 and Gc--!
      !---------------------------------------!
      IF (prop_l0 .LT. zero) THEN
        WRITE(7,*) "internal length must be positive. EXIT", prop_l0
        CALL XEXIT()
      END IF
      IF (prop_Gc .LT. zero) THEN
        WRITE(7,*) "fracture toughness must be positive. EXIT", prop_Gc
        CALL XEXIT()
      END IF
      !---------------------------------------!
    END SUBROUTINE CheckMaterialParameters

!------------------------------------------------------------------------------------

    SUBROUTINE umatPF(stress,svr,Ct,energy_elast,dissipat_plast,dissipat_creep,rpl,ddsddt,drplde,drpldt, &
                      stran,dstran,time,dtime,Temp,dTemp,predef_umat,dpred,cmname,ndi,nshr,ntens, &
                      num_svr,props_mat,num_matpar,coords_gp,Trafomat,pnewdt,intlen,F0,F1,jelem, &
                      npt,layer,kspt,kstep,kinc)
      ! compute generalised stresses and tangent

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE ABAModul
      USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
      USE FreeEnergyModule      
      USE BulkEnergyModule
      USE CrackSurfaceEnergyModule
      USE ViscousDissipationModule, ONLY: VDED => VDED_Liu, &
                                               d_VDED_d_phase => d_VDED_Liu_d_phase, &
                                               d_d_VDED_d_phase_d_phase => d_d_VDED_Liu_d_phase_d_phase
      USE DegradationFunctionModule, ONLY: degF => quad_degF, &
                                               d_degF_d_phase => d_quad_degF_d_phase, &
                                               d_d_degF_d_phase_d_phase => d_d_quad_degF_d_phase_d_phase
      !USE SplitEnergyModule
      
      USE TensorModule
	  USE SharedValues

      IMPLICIT NONE

      ! UMAT variables
      INTEGER(kind=AbqIK), INTENT(IN) :: ndi, nshr, ntens, num_svr, num_matpar, &
                                         jelem, npt, layer, kspt , kstep, kinc
      REAL(kind=AbqRK), INTENT(IN) :: stran(ntens), dstran(ntens), time(2), &
                                      dtime, Temp, dTemp, predef_umat(*), dpred(*), &
                                      coords_gp(3), Trafomat(3,3), intlen, F0(3,3), &
                                      F1(3,3), props_mat(num_matpar)
      REAL(kind=AbqRK), INTENT(INOUT) :: stress(ntens), svr(num_svr), &
                                         energy_elast, dissipat_plast, dissipat_creep, pnewdt
      REAL(kind=AbqRK), INTENT(OUT) :: Ct(ntens,ntens), rpl, ddsddt(ntens), &
                                       drplde(ntens), drpldt
      CHARACTER*45, INTENT(INOUT) :: cmname

      ! further variables
      INTEGER(kind=AbqIK) :: D, nphase, i1, i2, prop_solver
      INTEGER(kind=AbqIK) :: nHFEDpar
      INTEGER(kind=AbqIK) :: nCSEDpar
      INTEGER(kind=AbqIK) :: nVDEDpar
      INTEGER(kind=AbqIK) :: nMixedModepar
      INTEGER(kind=AbqIK) :: first_HFED_par_index
      INTEGER(kind=AbqIK) :: first_CSED_par_index
      INTEGER(kind=AbqIK) :: first_VDED_par_index
      INTEGER(kind=AbqIK) :: first_MixedMode_par_index
      
      ! pure debug Variables
      INTEGER(kind=AbqIK), PARAMETER :: n_debug_elem = 1
      INTEGER(kind=AbqIK), PARAMETER :: n_debug_inc = 10
      INTEGER(kind=AbqIK), PARAMETER :: n_debug_ip = 1
      
      INTEGER(kind=AbqIK) :: debug_elements(n_debug_elem)
      INTEGER(kind=AbqIK) :: debug_increments(n_debug_inc)
      INTEGER(kind=AbqIK) :: debug_ipoints(n_debug_ip)
      
      ! Debug Werte
      DATA debug_elements /1/
      DATA debug_increments /1,2,3,4,5,96,97,98,99,100/
      DATA debug_ipoints /1/
      ! end pure debug Variables

      REAL(kind=AbqRK) :: eps(3,3), eps_old(3,3) !eps old for staggered MP, H(eps_n)
      REAL(kind=AbqRK) :: delta, stran_per(ntens), stress_per(ntens), Ct_temp1(ntens,ntens), Ct_temp2(ntens,ntens), totalPotential,bulkedVar, csedVar, vdedVar, degVar, per, stress_num(ntens), H_sph, H_dev, Hn_sph, Hn_dev, ElasticEnergytens_sph, ElasticEnergytens_dev
      INTEGER(kind=AbqIK), ALLOCATABLE :: pos_p(:)
      REAL(kind=AbqRK), ALLOCATABLE :: phase(:), phase_old(:), grad_phase(:,:)
      REAL(kind=AbqRK), ALLOCATABLE :: parHFEDMatrix(:)
      REAL(kind=AbqRK), ALLOCATABLE :: parCSEDMatrix(:)
      REAL(kind=AbqRK), ALLOCATABLE :: parVDEDMatrix(:)
      REAL(kind=AbqRK), ALLOCATABLE :: parMixedModeMatrix(:) ! Contains GcI and GcII
      REAL(kind=AbqRK) :: Gc_MM, H_MM, GcI, GcII
      REAL(kind=AbqRK) :: prop_debug, prop_checkstress, prop_checktangent, prop_coupling
      LOGICAL :: numericalTangent, isSpherisymmetric, isNaN, NAN_Check, debug_Flag, check_stress, check_tangent, is_print_elem, is_print_inc, is_print_ip
      LOGICAL :: isMonolithic, isBFGS, isStaggered, isCouplingterms
      LOGICAL, SAVE :: info_printed = .FALSE.
      ! Initialization
      numericalTangent = .FALSE.; isSpherisymmetric = .FALSE.; isNaN = .FALSE.
      
      ! NAN CHECK
      NAN_Check = .TRUE.
      
      
      ! Solver choice
      isMonolithic = .FALSE.
      isBFGS = .FALSE.
      isStaggered = .FALSE.
      prop_solver     = INT(props_mat(10)) ! Solvervariable
	  IF (prop_solver .EQ. 1) THEN
		  isMonolithic = .TRUE.
	  ELSE IF (prop_solver .EQ. 2) THEN
		  isBFGS = .TRUE.
	  ELSE IF (prop_solver .EQ. 3) THEN
		  isStaggered = .TRUE.
	  ELSE 
		  WRITE(6,*) ' NO SOLVER FOUND, check solver variable (1 -> mon, 2 -> BFGS, 3 -> stag) '
		  CALL XEXIT()
	  END IF 
      
      ! Solver coupling terms?
      isCouplingterms = .FALSE. 
      prop_coupling = props_mat( 11 )
      IF (prop_coupling .GT. zero) THEN
		isCouplingterms = .TRUE.
      END IF
      
      ! Write Output for Debug Elements
      debug_Flag = .FALSE. 
      prop_debug = props_mat( 13 )
      IF (prop_debug .GT. zero) THEN
		debug_Flag = .TRUE.
      END IF
      
      ! compare num stress / analyt. stress
      check_stress = .FALSE.
      prop_checkstress = props_mat( 14 )
      IF (prop_checkstress .GT. zero) THEN
		check_stress = .TRUE.
      END IF
      
      ! compare num tang / analyt. tang
      check_tangent = .FALSE.
      prop_checktangent = props_mat( 15 )
      IF (prop_checktangent .GT. zero) THEN
		check_tangent = .TRUE.
      END IF

      ! global Variable to print selectively
      Incrementnumber = kinc
      Elementnumber = jelem
      Integrationpointnumber = npt
      
      ! Flags for Debug
      is_print_inc = .FALSE.
      is_print_elem = .FALSE.
      is_print_ip = .FALSE.
      
      IF (debug_Flag) THEN
		  ! See which printouts
		  DO i1=1,n_debug_elem
		   IF (jelem.EQ.debug_elements(i1)) THEN
			is_print_elem = .TRUE.
			EXIT
		   ENDIF
		  END DO
		  
		  DO i1=1,n_debug_inc
		   IF (kinc.EQ.debug_increments(i1)) THEN
			is_print_inc = .TRUE.
			EXIT
		   ENDIF
		  END DO
		  
		  DO i1=1,n_debug_ip
		   IF (npt.EQ.debug_ipoints(i1)) THEN
			is_print_ip = .TRUE.
			EXIT
		   ENDIF
		  END DO
		  
		  ! Global Variable to use the print flag in functions
		  is_print_elem_sub = is_print_elem
		  is_print_inc_sub = is_print_inc
		  is_print_ip_sub = is_print_ip
	  END IF
      
	  IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
	   WRITE(6,*) "---------------------------------------"
	   WRITE(6,*) " Debug-Ausgabe:"
	   WRITE(6,*) " Increment =", Incrementnumber
	   WRITE(6,*) " Element   =", Elementnumber
	   WRITE(6,*) " IP        =", Integrationpointnumber
	   WRITE(6,*) "---------------------------------------"
	  END IF

      
      ! dimension
      D = 0 ! dummy
      SELECT CASE(ndi+nshr)
      CASE(1)
        ! 1D
        D = 1
      CASE(3)
        ! 1D, spheri-symmetric
        D = 1
        isSpherisymmetric = .TRUE.
      CASE(4)
        ! 2D
        D = 2
      CASE(6)
        ! 3D
        D = 3
      END SELECT

      ! numerical tangent
      IF (props_mat(4) .NE. zero) numericalTangent = .TRUE.
       
      ! number of phases
      nphase = 1 ! damage variable
      !
      ! number of parameters in free energy
      nHFEDpar = 3
      !
      ! number of parameters in cracksurface energy density
      nCSEDpar = 2
      !
      ! number of parameters in viscous dissipation energy density
      nVDEDpar = 1
      !
      ! number of parameters in Mixed Mode Array
      nMixedModepar = 2
      !
      ALLOCATE(pos_p(nphase), phase(nphase), phase_old(nphase), grad_phase(nphase,3))
      ALLOCATE(parHFEDMatrix(nHFEDpar))
      ALLOCATE(parCSEDMatrix(nCSEDpar))
      ALLOCATE(parVDEDMatrix(nVDEDpar))
      ALLOCATE(parMixedModeMatrix(nMixedModepar))
      
      pos_p = 0
      phase = zero; phase_old = zero; grad_phase = zero
      parHFEDMatrix = zero
      parCSEDMatrix = zero
      parVDEDMatrix = zero
      parMixedModeMatrix = zero
	  !
      !
      ! Print Out Solver information
      !
      !
      IF (.NOT. info_printed) THEN
          WRITE(6,*) ''
          WRITE(6,*) '********* Solver *********'
          IF (isMonolithic) THEN
              WRITE(6,*) 'monolithic solver'
          ELSE IF (isBFGS) THEN
              WRITE(6,*) 'BFGS solver'
          ELSE IF (isStaggered) THEN
              WRITE(6,*) 'staggered solver'
          END IF 
          WRITE(6,*) ''
          IF (isCouplingterms) THEN
			WRITE(6,*) 'Coupling terms are active'
		  ELSE
			WRITE(6,*) 'No Coupling terms'
		  END IF
          WRITE(6,*) ''
          WRITE(6,*) '--------------------------'
          IF (debug_Flag) THEN
              WRITE(6,*) ''
              WRITE(6,*) 'PRINT ACTIVE'
              WRITE(6,*) ''
          ELSE 
              WRITE(6,*) ''
              WRITE(6,*) 'PRINT NOT ACTIVE'
              WRITE(6,*) ''
          END IF
          IF (check_stress) THEN
              WRITE(6,*) ''
              WRITE(6,*) 'Check Stress num/analyt'
              WRITE(6,*) ''
          END IF
          IF (check_tangent) THEN
              WRITE(6,*) ''
              WRITE(6,*) 'Check Tangent num/analyt'
              WRITE(6,*) ''
          END IF
          
          info_printed = .TRUE.
      END IF


      ! coordinate positions
      DO i1=1,nphase
        pos_p(i1) = ndi+nshr+(i1-1)*(D+1)+1
      END DO
      

      ! generalized kinematic measures
      eps = zero; phase = zero; grad_phase = zero
      ! strain tensor
      eps = strainCoordinates(D,ntens,stran,isSpherisymmetric)
      ! old strain tensor
      eps_old = strainCoordinates(D,ntens,stran-dstran,isSpherisymmetric)
      ! order parameter -- damage variable
      phase = phaseParameter(nphase,pos_p,ntens,stran)
      
      ! order parameter last increment for VDED
      phase_old(1) = svr(1)
      
      ! gradient of order parameter = gradient of damage variable
      grad_phase(:,:) = gradientOfPhaseParameter(nphase,D,pos_p,ntens,stran)    
        
      ! read material parameters from props_mat
      first_HFED_par_index = 5
      parHFEDMatrix = props_mat(first_HFED_par_index:first_HFED_par_index+nHFEDpar-1)
      first_CSED_par_index = 8
      parCSEDMatrix = props_mat(first_CSED_par_index:first_CSED_par_index+nCSEDpar-1)
      first_VDED_par_index = 16
      parVDEDMatrix = props_mat(first_VDED_par_index:first_VDED_par_index+nVDEDpar-1)
      first_MixedMode_par_index = 17
      parMixedModeMatrix = props_mat(first_MixedMode_par_index:first_MixedMode_par_index+nMixedModepar-1)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
	  ! NAN CHECKS
	  !
	  ! --- Prüfe eps (angenommen eps ist ein Array) ---
	  
   	  DO i1 = 1, SIZE(eps, 1)
		  DO i2 = 1, SIZE(eps, 2)
			  IF (eps(i1, i2) .NE. eps(i1, i2)) THEN
				  WRITE(7,*) 'eps UELstart:'
				  WRITE(7, '(3F15.8)') eps(1, :)
				  WRITE(7, '(3F15.8)') eps(2, :)
				  WRITE(7, '(3F15.8)') eps(3, :)
				  WRITE(7,*) 'NaN in eps bei Index UELstart (', i1, ',', i2, ')', ' Wert: ', eps(i1, i2)
				  WRITE(7,*) 'Increment: ', Incrementnumber
				  WRITE(7,*) 'Element: ', Elementnumber
				  WRITE(7,*) 'IP: ', Integrationpointnumber
				  !
			  END IF
		  END DO
	  END DO
	  
	  
	  ! --- Prüfe phase (kann Skalar oder Array sein) ---
	  DO i1 = 1, SIZE(phase)
		  IF (phase(i1) .NE. phase(i1)) THEN
			  WRITE(7,*) 'phase UELstart: ', phase(1)
			  WRITE(7,*) 'NaN in phase bei Index  UELstart', i1, ' Wert: ', phase(i1)
			  WRITE(7,*) 'Increment: ', Incrementnumber
			  WRITE(7,*) 'Element: ', Elementnumber
   			  WRITE(7,*) 'IP: ', Integrationpointnumber
   			  !
		  END IF
	  END DO
      ! Historyvariable 
      !
      ! H_n+1 (aus aktuellem eps_n+1)  -> BFGS
      ! H_n (aus letztem Konvergiertem eps_n)  -> staggered nach MP
      !
      IF (isMonolithic) THEN
		ElasticEnergytens_sph = HFEDtens_H_sph(eps,nHFEDpar,parHFEDMatrix) ! Monolithic
		ElasticEnergytens_dev = HFEDtens_H_dev(eps,nHFEDpar,parHFEDMatrix) ! Monolithic
	  ELSEIF (isBFGS) THEN
		ElasticEnergytens_sph = HFEDtens_H_sph(eps,nHFEDpar,parHFEDMatrix) ! BFGS
		ElasticEnergytens_dev = HFEDtens_H_dev(eps,nHFEDpar,parHFEDMatrix) ! BFGS
	  ELSEIF (isStaggered) THEN
		ElasticEnergytens_sph = HFEDtens_H_sph(eps_old,nHFEDpar,parHFEDMatrix) !staggered MP
		ElasticEnergytens_dev = HFEDtens_H_dev(eps_old,nHFEDpar,parHFEDMatrix) !staggered MP
      END IF
      !History Update
      Hn_sph = svr(30)
      Hn_dev = svr(31)
      H_sph = zero
      H_dev = zero
      
      
      IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
         write(6,*) '    === History - Feld ==='
         WRITE(6,*) "---------------------------------------"
	     WRITE(6,*) 'ElasticEnergytens_sph : ', ElasticEnergytens_sph 
	     WRITE(6,*) 'ElasticEnergytens_dev : ', ElasticEnergytens_dev
	     WRITE(6,*) 'ElasticEnergytens: ', ElasticEnergytens_sph + ElasticEnergytens_dev
         WRITE(6,*) 'History_n_sph: ', Hn_sph
         WRITE(6,*) 'History_n_dev: ', Hn_dev
		 WRITE(6,*) 'History_n: ', Hn_sph + Hn_dev
         WRITE(6,*) "----------------------------------------"
         WRITE(6,*) 'phase: ', phase
         WRITE(6,*) "----------------------------------------"
      END IF 
      
      IF (ElasticEnergytens_sph .GT. Hn_sph) THEN
        H_sph = ElasticEnergytens_sph
        H_dev = ElasticEnergytens_dev
      ELSE
        H_sph = Hn_sph
        H_dev = Hn_dev
      END IF 
      
      svr(30) = H_sph
      svr(31) = H_dev
      !
      ! compute all first derivatives
      
!~ 	  IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
!~ 	     WRITE(6,*) "---------------------------------------"
!~ 	     WRITE(6,*) "-------VOR STRESSES------"
!~ 	     WRITE(6,*) "---------------------------------------"
!~ 	     WRITE(6,*) "---------------------------------------"
!~ 	     WRITE(6,*) 'elastic stresses:'
!~ 	     WRITE(6,'(F20.15)') (stress(i1), i1=1,pos_p(1)-1)
!~ 	     WRITE(6,*) ''
!~ 	     WRITE(6,*) 'phase stresses:'
!~ 	     WRITE(6,'(F20.15)') (stress(pos_p(1)))
!~ 	     WRITE(6,*) 'gradphase stresses:'
!~ 	     WRITE(6,'(F20.15)') (stress(i1), i1=pos_p(1)+1, ntens)
!~ 	     WRITE(6,*) 'H:', H_sph + H_dev
!~ 	     WRITE(6,*) ''
!~ 	     WRITE(6,'(F20.15)') (stress(pos_p(1)))
!~ 	  END IF
      ! generalized stresses
      stress = stresses(D,ntens,nphase,pos_p,prop_solver,eps,eps_old,phase,phase_old,dtime,grad_phase,nHFEDpar,nCSEDpar,nVDEDpar,nMixedModepar,parHFEDMatrix,parCSEDMatrix,parVDEDMatrix,parMixedModeMatrix,H_sph,H_dev,isMonolithic, isBFGS, isStaggered, isCouplingterms)
      
      ! energetic quantities
      
!~ 	  IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
!~ 	     WRITE(6,*) "---------------------------------------"
!~ 	     WRITE(6,*) 'elastic stresses:'
!~ 	     WRITE(6,'(F20.15)') (stress(i1), i1=1,pos_p(1)-1)
!~ 	     WRITE(6,*) ''
!~ 	     WRITE(6,*) 'phase stresses:'
!~ 	     WRITE(6,'(F20.15)') (stress(pos_p(1)))
!~ 	     WRITE(6,*) 'gradphase stresses:'
!~ 	     WRITE(6,'(F20.15)') (stress(i1), i1=pos_p(1)+1, ntens)
!~ 	  END IF

      ! total potential
      totalPotential = potential(eps,nphase,phase,phase_old,dtime,grad_phase,nHFEDpar,nCSEDpar,nVDEDpar,parHFEDMatrix,parCSEDMatrix,parVDEDMatrix)

      ! Helmholtz free energy density of phase mix and Helmholtz energy density of interface
      !
      !WRITE(7,*) '----------------------'
      !WRITE(7,*) ' Save Bulk Energy SDV '
      !WRITE(7,*) '----------------------'
      !
      bulkedVar = bulkED(eps,phase(1),nHFEDpar,parHFEDMatrix) !ALLSE
      energy_elast = bulkedVar
      
      ! crack surface energy density
      csedVar = CSED(phase(1),grad_phase(1,:),nCSEDpar,parCSEDMatrix)
      dissipat_creep = csedVar
      
      ! degrading Function value
      degVar = degfunction(phase, nphase)

      ! viscous dissipation energy density
      vdedVar = VDED(phase(1),phase_old(1),dtime,nCSEDpar,parCSEDMatrix,nVDEDpar,parVDEDMatrix)
      
      
      ! dissipationEnergy
      dissipat_plast = totalPotential !ALLPD
            

      ! saved variables
      !
      ! stress coordinates
      svr(2:1+ndi+nshr) = stress(1:ndi+nshr)
      !
      ! strain coordinates
      svr(8:7+ndi+nshr) = stran(1:ndi+nshr)
      !
      ! integration point coordinates
      !
      ! phases
      svr(1) = phase(1)
      !
      ! reversible energy density
      svr(16) = energy_elast !bulkE
      svr(17) = dissipat_creep !CSE
      svr(18) = vdedVar !VDE
      svr(19) = dissipat_plast !bulkE + CSE + VDE
      ! Gradient phase
      svr(21:23) = grad_phase(1,:)
      
      
      ! generalised tangent
      IF (numericalTangent) THEN

        ! numerical tangent obtained with perturbation of stresses
        delta = 1.d-6
        Ct = zero
        DO i1=1,ntens
          stran_per = stran
          stran_per(i1) = stran_per(i1) + delta
          ! generalized kinematic measures
          eps = zero; phase = zero; grad_phase = zero
          ! strain tensor
          eps = strainCoordinates(D,ntens,stran_per,isSpherisymmetric)
          ! order parameter -- damage variable
          phase = phaseParameter(nphase,pos_p,ntens,stran_per)
          ! gradient of order parameter = gradient of damage variable
          grad_phase(:,:) = gradientOfPhaseParameter(nphase,D,pos_p,ntens,stran_per)
          !
          ! generalised stresses
          stress_per = stresses(D,ntens,nphase,pos_p,prop_solver,eps,eps_old,phase,phase_old,dtime,grad_phase,nHFEDpar,nCSEDpar,nVDEDpar,nMixedModepar,parHFEDMatrix,parCSEDMatrix,parVDEDMatrix,parMixedModeMatrix,H_sph,H_dev,isMonolithic,isBFGS,isStaggered,isCouplingterms)
          ! tangent
          Ct(i1,1:ntens) = (stress_per(1:ntens)-stress(1:ntens))/delta
        END DO

      ELSE
        ! analytical material tangent
        Ct = tangent(D,ntens,nphase,pos_p,prop_solver,eps,eps_old,phase,phase_old,dtime,grad_phase,nHFEDpar,nCSEDpar,nVDEDpar,nMixedModepar,parHFEDMatrix,parCSEDMatrix,parVDEDMatrix,parMixedModeMatrix,H_sph,H_dev,isMonolithic,isBFGS,isStaggered,isCouplingterms)
      END IF
	  
	  
!~ 	  IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
!~ 		 WRITE(6,*), "---------------------------------------"
!~ 		 WRITE(6,*) 'Ct analytisch:'
!~ 		 DO i1 = 1, ntens
!~ 			WRITE(6,'(7(F16.4,1X))') (Ct(i1,i2), i2=1,ntens)
!~ 		 END DO
!~ 		 WRITE(6,*) ''
!~ 		 WRITE(6,*) "----------------------------------------"
!~ 	  END IF
	  
		IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
		   WRITE(6,*), "---------------------------------------"
		   WRITE(6,*) 'Ct analytisch:'
		   WRITE(6,*) ''
		   
		   ! Ausgabe des 4x4 Blocks
		   WRITE(6,*) 'Block 1:4, 1:4:'
		   DO i1 = 1, 4
			  DO i2 = 1, 4
				 WRITE(6,'(A,I1,A,I1,A,F16.8)') 'Ct(', i1, ',', i2, ') = ', Ct(i1,i2)
			  END DO
		   END DO
		   WRITE(6,*) ''
		   WRITE(6,*) '---------------------------------------'
		   WRITE(6,*) 'phase Ct analytisch: '
		   WRITE(6,'(F16.10)') Ct(pos_p(1),pos_p(1))
		   WRITE(6,*) 'gradphase Ct: '
		   WRITE(6,'(F16.10)') Ct(pos_p(1)+1:ntens,pos_p(1)+1:ntens)
		   WRITE(6,*) '----------------------------------------'
		   WRITE(6,*) ''
!~ 		   ! Ausgabe des 3x3 Blocks (falls ntens >= 7)
!~ 		   IF (ntens >= 7) THEN
!~ 			  WRITE(6,*) 'Block 5:7, 5:7:'
!~ 			  DO i1 = 5, 7
!~ 				 DO i2 = 5, 7
!~ 					WRITE(6,'(A,I1,A,I1,A,F16.10)') 'Ct(', i1, ',', i2, ') = ', Ct(i1,i2)
!~ 				 END DO
!~ 			  END DO
!~ 			  WRITE(6,*) ''
!~ 		   END IF
		END IF
	  
	  !!!!!!!!!!!!
	  !
	  ! Check Comparison num stress / analytical stress
	  ! Mit central diff quotient
	  !!!!!!!!!!!!
	  !
	  IF (check_stress) THEN
	    CALL Check_Stress_Comparison(stress, stran, D, ntens, isSpherisymmetric, &
                                nphase, pos_p, phase_old, dtime, &
                                nHFEDpar, nCSEDpar, nVDEDpar, nMixedModepar, parHFEDMatrix, &
                                parCSEDMatrix, parVDEDMatrix, parMixedModeMatrix, H_sph, H_dev, isMonolithic, isBFGS, isStaggered, isCouplingterms)
	  END IF
	  !!!!!!!!!!!!
	  !
	  ! Check Comparison num tangent / analytical tangent
	  !
	  !!!!!!!!!!!!
	  IF (check_tangent) THEN
		CALL Check_Tangent_Comparison(Ct, stress, stran, D, ntens, isSpherisymmetric, &
								   nphase, pos_p, prop_solver, phase_old, dtime, &
								   nHFEDpar, nCSEDpar, nVDEDpar, nMixedModepar, parHFEDMatrix, &
								   parCSEDMatrix, parVDEDMatrix, parMixedModeMatrix, H_sph, H_dev, isMonolithic, isBFGS, isStaggered, isCouplingterms)
	  END IF
	  !!!!!!!!!

      svr(25) = Ct(1,1)
      svr(26) = Ct(2,2)
      svr(27) = Ct(3,3)
      svr(28) = degVar

      ! dummy output
      rpl = zero; ddsddt = zero; drplde = zero; drpldt = zero
      
      !
	  ! NAN CHECKS
	  !
	  IF (NAN_Check) THEN
		  IF (energy_elast .NE. energy_elast) THEN
			WRITE(7,*) 'energy_elast NAN: ', energy_elast
			!
			! Function Variables
			!
			WRITE(7,*) 'eps: ', eps
			WRITE(7,*) 'phase: ', phase(1)
			WRITE(7,*) 'degF: ', degVar
			!
			WRITE(7,*) 'tangent pure mechanics: '
			DO i1 = 1, 4
				WRITE(7,'(6F12.1)') (Ct(i1,i2), i2=1,4)
			END DO
			
			
			WRITE(7,*) "Increment: ", Incrementnumber
			WRITE(7,*) "Element: ", Elementnumber
			WRITE(7,*) "IP: ", Integrationpointnumber
			!
			CALL XEXIT()
		  END IF
		  IF (dissipat_creep .NE. dissipat_creep) THEN
			WRITE(7,*) 'dissipat_creep NAN: ', dissipat_creep 
			!
			! Function Variables
			!
			WRITE(7,*) 'phase: ', phase(1)
			WRITE(7,*) 'Gradient phase: ', grad_phase(1,:)
			!
			WRITE(7,*) "Increment: ", Incrementnumber
			WRITE(7,*) "Element: ", Elementnumber
			WRITE(7,*) "IP: ", Integrationpointnumber
			!
			CALL XEXIT()
		  END IF
		  IF (vdedVar .NE. vdedVar) THEN
			WRITE(7,*) 'vdedVar NAN: ', vdedVar
			!
			! Function Variables
			!
			WRITE(7,*) 'phase: ', phase(1)
			WRITE(7,*) 'phase_old: ', phase_old(1)
			WRITE(7,*) 'dtime: ', dtime
			!
			WRITE(7,*) "Increment: ", Incrementnumber
			WRITE(7,*) "Element: ", Elementnumber
			WRITE(7,*) "IP: ", Integrationpointnumber
			!
			CALL XEXIT()
		  END IF

      END IF
      

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DEALLOCATE(pos_p, phase, phase_old, grad_phase)
      DEALLOCATE(parHFEDMatrix)
      DEALLOCATE(parCSEDMatrix)
      DEALLOCATE(parVDEDMatrix)

      CONTAINS
      
      END SUBROUTINE umatPF

!------------------------------------------------------------------------------------

      PURE REAL(kind=AbqRK) FUNCTION degfunction(phase, nphase)
      ! computes crack surface energy density

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS
        USE TensorModule
        USE DegradationFunctionModule, ONLY: degF => quad_degF, &
                                       d_degF_d_phase => d_quad_degF_d_phase, &
                                       d_d_degF_d_phase_d_phase => d_d_quad_degF_d_phase_d_phase

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: nphase
        REAL(kind=AbqRK), INTENT(IN) :: phase(nphase)

        degfunction = degF(phase(1))

      END FUNCTION degfunction

!------------------------------------------------------------------------------------


      REAL(kind=AbqRK) FUNCTION potential(eps,nphase,phase,phase_old,dtime,grad_phase,nHFEDpar,nCSEDpar,nVDEDpar,parHFEDMatrix,parCSEDMatrix,parVDEDMatrix)
      ! computes total potential

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS
        USE BulkEnergyModule
        USE CrackSurfaceEnergyModule
        USE ViscousDissipationModule, ONLY: VDED => VDED_Liu, &
                                               d_VDED_d_phase => d_VDED_Liu_d_phase, &
                                               d_d_VDED_d_phase_d_phase => d_d_VDED_Liu_d_phase_d_phase

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
        INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
        INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
        INTEGER(kind=AbqIK), INTENT(IN) :: nphase
        REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrix(nHFEDpar)
        REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrix(nCSEDpar)
        REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrix(nVDEDpar)
        REAL(kind=AbqRK), INTENT(IN) :: eps(3,3), phase(nphase), phase_old(nphase), dtime, grad_phase(nphase,3)
		!
		!WRITE(7,*) '-----------------'
		!WRITE(7,*) '- Potential     -'
		!WRITE(7,*) '-----------------'
		!
        potential = bulkED(eps,phase(1),nHFEDpar,parHFEDMatrix) + CSED(phase(1),grad_phase(1,:),nCSEDpar,parCSEDMatrix) + VDED(phase(1),phase_old(1),dtime,nCSEDpar,parCSEDMatrix,nVDEDpar,parVDEDMatrix)

      END FUNCTION potential

!------------------------------------------------------------------------------------

      PURE REAL(kind=AbqRK) FUNCTION strainCoordinates(D,ntens,s,isSpherisymmetric)
      ! extract strain coordinates

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: D, ntens
        REAL(kind=AbqRK), INTENT(IN) :: s(ntens)
        LOGICAL, INTENT(IN) :: isSpherisymmetric
        DIMENSION strainCoordinates(3,3)

        strainCoordinates = zero
        !
        ! strain tensor
        strainCoordinates(1,1) = s(1)
        IF ( (D .GE. 2) .OR. ( (D .EQ. 1) .AND. (isSpherisymmetric) ) ) THEN
          strainCoordinates(2,2) = s(2)
          strainCoordinates(3,3) = s(3)
        END IF
        IF (D .GE. 2) THEN
          strainCoordinates(1,2) = s(4)/two
          strainCoordinates(2,1) = strainCoordinates(1,2)
        END IF
        IF (D .EQ. 3) THEN
          strainCoordinates(1,3) = s(5)/two
          strainCoordinates(2,3) = s(6)/two
          strainCoordinates(3,1:2) = strainCoordinates(1:2,3)
        END IF

      END FUNCTION strainCoordinates

!------------------------------------------------------------------------------------

      PURE REAL(kind=AbqRK) FUNCTION strainCoordinates_backward(D,ntens,S,isSpherisymmetric)
      ! extract strain coordinates

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: D, ntens
        REAL(kind=AbqRK), INTENT(IN) :: S(3,3)
        LOGICAL, INTENT(IN) :: isSpherisymmetric
        DIMENSION strainCoordinates_backward(ntens)

        strainCoordinates_backward = zero
        !
        ! strain tensor
        strainCoordinates_backward(1) = S(1,1)
        IF ( (D .GE. 2) .OR. ( (D .EQ. 1) .AND. (isSpherisymmetric) ) ) THEN
          strainCoordinates_backward(2) = S(2,2)
          strainCoordinates_backward(3) = S(3,3)
        END IF
        IF (D .GE. 2) THEN
          strainCoordinates_backward(4) = S(1,2) * 2.0
        END IF
        IF (D .EQ. 3) THEN
          strainCoordinates_backward(5) = S(1,3) * 2.0
          strainCoordinates_backward(6)= S(2,3) * 2.0
        END IF

      END FUNCTION strainCoordinates_backward

!------------------------------------------------------------------------------------

      PURE REAL(kind=AbqRK) FUNCTION phaseParameter(nphase,pos_p,ntens,s)
      ! extract phase parameter

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: nphase, pos_p(nphase), ntens
        REAL(kind=AbqRK), INTENT(IN) :: s(ntens)
        INTEGER(kind=AbqIK) :: i1
        DIMENSION phaseParameter(nphase)

        phaseParameter = zero
        FORALL (i1=1:nphase)
          phaseParameter(i1) = s(pos_p(i1))
        END FORALL

      END FUNCTION phaseParameter

!------------------------------------------------------------------------------------

      PURE REAL(kind=AbqRK) FUNCTION gradientOfPhaseParameter(nphase,D,pos_p,ntens,s)
      ! extract gradient of phase parameter

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: nphase, D, pos_p(nphase), ntens
        REAL(kind=AbqRK), INTENT(IN) :: s(ntens)
        INTEGER(kind=AbqIK) :: i1, i2
        DIMENSION gradientOfPhaseParameter(nphase,3)

        gradientOfPhaseParameter = zero
        FORALL (i1=1:nphase)
          gradientOfPhaseParameter(i1,1:D) = s(pos_p(i1)+1:pos_p(i1)+D)
        END FORALL

      END FUNCTION gradientOfPhaseParameter

!------------------------------------------------------------------------------------

      REAL(kind=AbqRK) FUNCTION stresses(D,ntens,nphase,pos_p, prop_solver,eps,eps_old,phase,phase_old,dtime,grad_phase,nHFEDpar,nCSEDpar,nVDEDpar,nMixedModepar,parHFEDMatrix,parCSEDMatrix,parVDEDMatrix,parMixedModeMatrix,H_sph,H_dev,isMonolithic, isBFGS, isStaggered, isCouplingterms)
      ! generalised stresses, modified Voigt notation

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS
        USE SharedValues
        USE BulkEnergyModule
        USE SplitEnergyModule
        USE CrackSurfaceEnergyModule
        USE ViscousDissipationModule, ONLY: VDED => VDED_Liu, &
                                               d_VDED_d_phase => d_VDED_Liu_d_phase, &
                                               d_d_VDED_d_phase_d_phase => d_d_VDED_Liu_d_phase_d_phase
        
        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: D, ntens, nphase, pos_p(nphase), prop_solver
        INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
        INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
        INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
        INTEGER(kind=AbqIK), INTENT(IN) :: nMixedModepar
        REAL(kind=AbqRK), INTENT(IN) :: eps(3,3), eps_old(3,3), phase(nphase), phase_old(nphase), grad_phase(nphase,3), dtime, H_sph, H_dev
        REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrix(nHFEDpar) 
        REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrix(nCSEDpar) 
        REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrix(nVDEDpar)
        REAL(kind=AbqRK), INTENT(IN) :: parMixedModeMatrix(nMixedModepar)
        REAL(kind=AbqRK), ALLOCATABLE :: parCSED_MM_Matrix(:)
        REAL(kind=AbqRK) :: Gc_MM, H_MM
        LOGICAL, INTENT(IN) :: isMonolithic, isBFGS, isStaggered, isCouplingterms
        REAL(kind=AbqRK) :: d_bulkEDMat_d_eps(3,3)
        REAL(kind=AbqRK) :: temp(3)
        INTEGER(kind=AbqIK) :: i1, pos_i1
        DIMENSION stresses(ntens)


		ALLOCATE(parCSED_MM_Matrix(nCSEDpar))
		!
		!WRITE(7,*) '-----------------'
		!WRITE(7,*) ' Stresses	inc: ', temp_inc
		!WRITE(7,*) '-----------------'
	    !
		
		
        ! first derivatives of bulk energy density
        d_bulkEDMat_d_eps = zero
        IF (isMonolithic) THEN
            !monolithic
            d_bulkEDMat_d_eps = d_bulkED_d_eps(eps,phase(1),nHFEDpar,parHFEDMatrix)
        ELSEIF (isBFGS) THEN
            !BFGS
            d_bulkEDMat_d_eps = d_bulkED_d_eps(eps,phase(1),nHFEDpar,parHFEDMatrix)
        ELSEIF (isStaggered) THEN
            !staggered Martinez-Paneda
            d_bulkEDMat_d_eps = d_bulkED_d_eps(eps,phase_old(1),nHFEDpar,parHFEDMatrix)
        ENDIF
        
        !
        stresses=zero
        !
        ! eps: contributions of bulkED
        stresses(1) = d_bulkEDMat_d_eps(1,1)
        stresses(2) = d_bulkEDMat_d_eps(2,2)
        stresses(3) = d_bulkEDMat_d_eps(3,3)
        IF (D .GE. 2) THEN
          stresses(4) = d_bulkEDMat_d_eps(1,2)
        END IF
        IF (D .EQ. 3) THEN
          stresses(5) = d_bulkEDMat_d_eps(1,3)
          stresses(6) = d_bulkEDMat_d_eps(2,3)
        END IF
        !
        !
        ! For Stress Function to change Gc
		parCSED_MM_Matrix(:) = parCSEDMatrix(:)
        !
		! Mixed Mode          
		!
		Gc_MM = Gceff(H_sph, H_dev, nMixedModepar,parMixedModeMatrix)
		H_MM = Heff(H_sph, H_dev, nMixedModepar,parMixedModeMatrix)
		!
		IF (is_print_elem_sub .AND. is_print_inc_sub .AND. is_print_ip_sub) THEN
			WRITE(6,*)
			WRITE(6,*) '=============================================='
			WRITE(6,*) ' DEBUG: parCSED_MM_Matrix BEFORE Mixed Mode'
			WRITE(6,*) '=============================================='
			WRITE(6,*) 'parCSED_MM_Matrix(1) = ', parCSED_MM_Matrix(1)
			WRITE(6,*) 'parCSED_MM_Matrix(2) = ', parCSED_MM_Matrix(2)

			WRITE(6,*)
			WRITE(6,*) '----------------------------------------------'
			WRITE(6,*) ' DEBUG: Arguments for Gceff / Heff'
			WRITE(6,*) '----------------------------------------------'
			WRITE(6,*) 'H_sph            = ', H_sph
			WRITE(6,*) 'H_dev            = ', H_dev
			WRITE(6,*) 'nMixedModepar    = ', nMixedModepar
			WRITE(6,*) 'parMixedModeMatrix = ', parMixedModeMatrix(:)

			WRITE(6,*)
			WRITE(6,*) '----------------------------------------------'
			WRITE(6,*) ' DEBUG: Results from Gceff / Heff'
			WRITE(6,*) '----------------------------------------------'
			WRITE(6,*) 'Gc_MM = ', Gc_MM
			WRITE(6,*) 'H_MM  = ', H_MM

			WRITE(6,*)
			WRITE(6,*) '----------------------------------------------'
			WRITE(6,*) ' DEBUG: Updating parCSED_MM_Matrix(2)'
			WRITE(6,*) '----------------------------------------------'
			WRITE(6,*) 'Old parCSED_MM_Matrix(2) = ', parCSED_MM_Matrix(2)
		END IF	
			
		! Updating the parCSEDMatrix Vector for Energy Calculations
		parCSED_MM_Matrix(2) = Gc_MM
		
		IF (is_print_elem_sub .AND. is_print_inc_sub .AND. is_print_ip_sub) THEN
			WRITE(6,*) 'New parCSED_MM_Matrix(2) = ', parCSED_MM_Matrix(2)

			WRITE(6,*)
			WRITE(6,*) '=============================================='
			WRITE(6,*) ' DEBUG: FINAL parCSED_MM_Matrix ARRAY'
			WRITE(6,*) '=============================================='
			WRITE(6,*) parCSED_MM_Matrix(:)
			WRITE(6,*)
		END IF
        !
        ! phase parameters
        !
        DO i1=1,nphase
          pos_i1 = pos_p(i1)
          !         
          ! phase: contributions of bulk energy + crack surface energy
          IF (isMonolithic) THEN
            !monolithic, staggered is the same, due to H defined earlier, eps in d_bulkED_d_phase_H is not used in the function
            !stresses(pos_i1) = d_bulkED_d_phase_H(eps,phase(1),nHFEDpar,parHFEDMatrix,H_sph,H_dev) + d_CSED_d_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSEDMatrix) +                   d_VDED_d_phase(phase(1),phase_old(1),dtime,nCSEDpar,parCSEDMatrix,nVDEDpar,parVDEDMatrix)
			stresses(pos_i1) = d_bulkED_d_phase_H(eps,phase(1),nHFEDpar,parHFEDMatrix,H_MM) + d_CSED_d_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSED_MM_Matrix) + d_VDED_d_phase(phase(1),phase_old(1),dtime,nCSEDpar,parCSED_MM_Matrix,nVDEDpar,parVDEDMatrix)
          ELSE IF (isBFGS) THEN
            !BFGS
            !stresses(pos_i1) = d_bulkED_d_phase_H(eps,phase(1),nHFEDpar,parHFEDMatrix,H_sph,H_dev) + d_CSED_d_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSEDMatrix) +                   d_VDED_d_phase(phase(1),phase_old(1),dtime,nCSEDpar,parCSEDMatrix,nVDEDpar,parVDEDMatrix)
			stresses(pos_i1) = d_bulkED_d_phase_H(eps,phase(1),nHFEDpar,parHFEDMatrix,H_MM) + d_CSED_d_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSED_MM_Matrix) + d_VDED_d_phase(phase(1),phase_old(1),dtime,nCSEDpar,parCSED_MM_Matrix,nVDEDpar,parVDEDMatrix)
			
			IF (is_print_elem_sub .AND. is_print_inc_sub .AND. is_print_ip_sub) THEN
				WRITE(6,*)
				WRITE(6,*) '=============================================='
				WRITE(6,*) ' DEBUG: STRESS CONTRIBUTIONS (no temp vars)'
				WRITE(6,*) '=============================================='

				!======================
				! BULK ED
				!======================

				WRITE(6,*)
				WRITE(6,*) '----------------------------------------------'
				WRITE(6,*) ' BULK ED TERM'
				WRITE(6,*) '----------------------------------------------'
				WRITE(6,*) 'eps           = ', eps
				WRITE(6,*) 'phase(1)      = ', phase(1)
				WRITE(6,*) 'nHFEDpar      = ', nHFEDpar
				WRITE(6,*) 'parHFEDMatrix = ', parHFEDMatrix(:)
				WRITE(6,*) 'H_MM          = ', H_MM

				WRITE(6,*) 'd_bulkED_d_phase_H = ', &
					 d_bulkED_d_phase_H(eps,phase(1),nHFEDpar,parHFEDMatrix,H_MM)

				!======================
				! CSED
				!======================

				WRITE(6,*)
				WRITE(6,*) '----------------------------------------------'
				WRITE(6,*) ' CSED TERM'
				WRITE(6,*) '----------------------------------------------'
				WRITE(6,*) 'phase(1)          = ', phase(1)
				WRITE(6,*) 'grad_phase(1,:)   = ', grad_phase(1,:)
				WRITE(6,*) 'nCSEDpar          = ', nCSEDpar
				WRITE(6,*) 'parCSED_MM_Matrix = ', parCSED_MM_Matrix(:)

				WRITE(6,*) 'd_CSED_d_phase = ', &
					 d_CSED_d_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSED_MM_Matrix)

				!======================
				! VDED
				!======================

				WRITE(6,*)
				WRITE(6,*) '----------------------------------------------'
				WRITE(6,*) ' VDED TERM'
				WRITE(6,*) '----------------------------------------------'
				WRITE(6,*) 'phase(1)      = ', phase(1)
				WRITE(6,*) 'phase_old(1)  = ', phase_old(1)
				WRITE(6,*) 'dtime         = ', dtime
				WRITE(6,*) 'nCSEDpar      = ', nCSEDpar
				WRITE(6,*) 'parCSED_MM_Matrix = ', parCSED_MM_Matrix(:)
				WRITE(6,*) 'nVDEDpar      = ', nVDEDpar
				WRITE(6,*) 'parVDEDMatrix = ', parVDEDMatrix(:)

				WRITE(6,*) 'd_VDED_d_phase = ', &
					 d_VDED_d_phase(phase(1),phase_old(1),dtime,nCSEDpar,parCSED_MM_Matrix,nVDEDpar,parVDEDMatrix)

				!======================
				! TOTAL
				!======================

				WRITE(6,*)
				WRITE(6,*) '=============================================='
				WRITE(6,*) ' TOTAL STRESS'
				WRITE(6,*) '=============================================='
				WRITE(6,*) 'stresses(pos_i1) = ', stresses(pos_i1)
				WRITE(6,*)
			END IF
          
          
          ELSE IF (isStaggered) THEN
            !staggered Martinez-Paneda
            !stresses(pos_i1) = d_bulkED_d_phase_H(eps_old,phase(1),nHFEDpar,parHFEDMatrix,H_sph,H_dev) + d_CSED_d_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSEDMatrix) +                   d_VDED_d_phase(phase(1),phase_old(1),dtime,nCSEDpar,parCSEDMatrix,nVDEDpar,parVDEDMatrix)
			stresses(pos_i1) = d_bulkED_d_phase_H(eps_old,phase(1),nHFEDpar,parHFEDMatrix,H_MM) + d_CSED_d_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSED_MM_Matrix) + d_VDED_d_phase(phase(1),phase_old(1),dtime,nCSEDpar,parCSED_MM_Matrix,nVDEDpar,parVDEDMatrix)
		  END IF
          
          !
          ! grad_phase: contributions of crack surface energy
          temp(:) = d_CSED_d_grad_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSEDMatrix)
          stresses(pos_i1+1:pos_i1+D) = temp(1:D)
        END DO
        
        
        DEALLOCATE(parCSED_MM_Matrix)

        

      END FUNCTION stresses

!------------------------------------------------------------------------------------

      REAL(kind=AbqRK) FUNCTION tangent(D,ntens,nphase,pos_p,prop_solver,eps,eps_old,phase,phase_old,dtime,grad_phase,nHFEDpar,nCSEDpar,nVDEDpar,nMixedModepar,parHFEDMatrix,parCSEDMatrix,parVDEDMatrix,parMixedModeMatrix,H_sph,H_dev,isMonolithic, isBFGS, isStaggered, isCouplingterms)
      ! generalised tangent, modified Voigt notation

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS
        USE BulkEnergyModule
        USE CrackSurfaceEnergyModule
        USE SplitEnergyModule
        USE ViscousDissipationModule, ONLY: VDED => VDED_Liu, &
                                               d_VDED_d_phase => d_VDED_Liu_d_phase, &
                                               d_d_VDED_d_phase_d_phase => d_d_VDED_Liu_d_phase_d_phase

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: D, ntens, nphase, pos_p(nphase), prop_solver
        INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar, nCSEDpar, nVDEDpar, nMixedModepar
        REAL(kind=AbqRK), INTENT(IN) :: eps(3,3),eps_old(3,3), phase(nphase), phase_old(nphase), grad_phase(nphase,3), dtime, H_sph, H_dev
        REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrix(nHFEDpar) 
        REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrix(nCSEDpar) 
        REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrix(nVDEDpar) 
        REAL(kind=AbqRK), INTENT(IN) :: parMixedModeMatrix(nMixedModepar)
        REAL(kind=AbqRK), ALLOCATABLE :: parCSED_MM_Matrix(:)
        REAL(kind=AbqRK) :: Gc_MM, H_MM
        LOGICAL, INTENT(IN) :: isMonolithic, isBFGS, isStaggered, isCouplingterms
        REAL(kind=AbqRK) :: d_d_bulkEDMat_d_eps_d_eps(3,3,3,3)
        REAL(kind=AbqRK) :: d_d_bulkEDMat_d_eps_d_phase(3,3), d_d_bulkEDMat_d_phase_d_eps(3,3)
        REAL(kind=AbqRK) :: d_d_bulkED_CSED_VDED_d_phase_d_phase
        REAL(kind=AbqRK) :: temp1(3), temp2(3,3)
        INTEGER(kind=AbqIK) :: i1, i2, pos_i1, pos_i2
        DIMENSION tangent(ntens,ntens)

		ALLOCATE(parCSED_MM_Matrix(nCSEDpar))

        ! second derivatives of Helmholtz free energy density
        d_d_bulkEDMat_d_eps_d_eps = zero
        
        IF (isMonolithic) THEN
            !monolithisch
            d_d_bulkEDMat_d_eps_d_eps = d_d_bulkED_d_eps_d_eps(eps,phase(1),nHFEDpar,parHFEDMatrix)
        ELSEIF (isBFGS) THEN
            !
            d_d_bulkEDMat_d_eps_d_eps = d_d_bulkED_d_eps_d_eps(eps,phase(1),nHFEDpar,parHFEDMatrix)
        ELSEIF (isStaggered) THEN
            !staggered Martinez-Paneda
            d_d_bulkEDMat_d_eps_d_eps = d_d_bulkED_d_eps_d_eps(eps,phase_old(1),nHFEDpar,parHFEDMatrix)
        ENDIF

        tangent=zero
        !
        ! eps-eps: contributions of HFED
        
        tangent(1,1) = d_d_bulkEDMat_d_eps_d_eps(1,1,1,1)
        tangent(1,2) = d_d_bulkEDMat_d_eps_d_eps(1,1,2,2)
        tangent(1,3) = d_d_bulkEDMat_d_eps_d_eps(1,1,3,3)
        tangent(2,1) = d_d_bulkEDMat_d_eps_d_eps(2,2,1,1)
        tangent(2,2) = d_d_bulkEDMat_d_eps_d_eps(2,2,2,2)
        tangent(2,3) = d_d_bulkEDMat_d_eps_d_eps(2,2,3,3)
        tangent(3,1) = d_d_bulkEDMat_d_eps_d_eps(3,3,1,1)
        tangent(3,2) = d_d_bulkEDMat_d_eps_d_eps(3,3,2,2)
        tangent(3,3) = d_d_bulkEDMat_d_eps_d_eps(3,3,3,3)
        IF (D .GE. 2) THEN
          tangent(1,4) = d_d_bulkEDMat_d_eps_d_eps(1,1,1,2)
          tangent(2,4) = d_d_bulkEDMat_d_eps_d_eps(2,2,1,2)
          tangent(3,4) = d_d_bulkEDMat_d_eps_d_eps(3,3,1,2)
          tangent(4,1) = d_d_bulkEDMat_d_eps_d_eps(1,2,1,1)
          tangent(4,2) = d_d_bulkEDMat_d_eps_d_eps(1,2,2,2)
          tangent(4,3) = d_d_bulkEDMat_d_eps_d_eps(1,2,3,3)
          tangent(4,4) = d_d_bulkEDMat_d_eps_d_eps(1,2,1,2)
        END IF
        IF (D .EQ. 3) THEN
          tangent(1,5) = d_d_bulkEDMat_d_eps_d_eps(1,1,1,3)
          tangent(1,6) = d_d_bulkEDMat_d_eps_d_eps(1,1,2,3)
          tangent(2,5) = d_d_bulkEDMat_d_eps_d_eps(2,2,1,3)
          tangent(2,6) = d_d_bulkEDMat_d_eps_d_eps(2,2,2,3)
          tangent(3,5) = d_d_bulkEDMat_d_eps_d_eps(3,3,1,3)
          tangent(3,6) = d_d_bulkEDMat_d_eps_d_eps(3,3,2,3)
          tangent(4,5) = d_d_bulkEDMat_d_eps_d_eps(1,2,1,3)
          tangent(4,6) = d_d_bulkEDMat_d_eps_d_eps(1,2,2,3)
          tangent(5,1) = d_d_bulkEDMat_d_eps_d_eps(1,3,1,1)
          tangent(5,2) = d_d_bulkEDMat_d_eps_d_eps(1,3,2,2)
          tangent(5,3) = d_d_bulkEDMat_d_eps_d_eps(1,3,3,3)
          tangent(5,4) = d_d_bulkEDMat_d_eps_d_eps(1,3,1,2)
          tangent(5,5) = d_d_bulkEDMat_d_eps_d_eps(1,3,1,3)
          tangent(5,6) = d_d_bulkEDMat_d_eps_d_eps(1,3,2,3)
          tangent(6,1) = d_d_bulkEDMat_d_eps_d_eps(2,3,1,1)
          tangent(6,2) = d_d_bulkEDMat_d_eps_d_eps(2,3,2,2)
          tangent(6,3) = d_d_bulkEDMat_d_eps_d_eps(2,3,3,3)
          tangent(6,4) = d_d_bulkEDMat_d_eps_d_eps(2,3,1,2)
          tangent(6,5) = d_d_bulkEDMat_d_eps_d_eps(2,3,1,3)
          tangent(6,6) = d_d_bulkEDMat_d_eps_d_eps(2,3,2,3)
        END IF
        !
        ! phase parameters
        !
        !
        !
        ! For Stress Function to change Gc
		parCSED_MM_Matrix(:) = parCSEDMatrix(:)
        !
		! Mixed Mode          
		!
		Gc_MM = Gceff(H_sph, H_dev, nMixedModepar,parMixedModeMatrix)
		H_MM = Heff(H_sph, H_dev, nMixedModepar,parMixedModeMatrix)
		!
		! Updating the parCSEDMatrix Vector for Energy Calculations
		parCSED_MM_Matrix(2) = Gc_MM
		!
        DO i1=1,nphase
          pos_i1 = pos_p(i1)
          !
          ! eps-phase: contributions of bulkEnergy
          ! Initialization
          d_d_bulkEDMat_d_eps_d_phase = zero
          

          IF (isCouplingterms) THEN
			  IF (isMonolithic) THEN
				!monolithisch
				d_d_bulkEDMat_d_eps_d_phase = d_d_bulkED_d_eps_d_phase(eps,phase(1),nHFEDpar,parHFEDMatrix)
			  ELSEIF (isBFGS) THEN
				!BFGS
				d_d_bulkEDMat_d_eps_d_phase = d_d_bulkED_d_eps_d_phase(eps,phase(1),nHFEDpar,parHFEDMatrix)
			  ELSEIF (isStaggered) THEN
				!staggered Martinez-Paneda
				d_d_bulkEDMat_d_eps_d_phase = d_d_bulkED_d_eps_d_phase(eps,phase_old(1),nHFEDpar,parHFEDMatrix)
			  ENDIF
		  ELSE	
			  d_d_bulkEDMat_d_eps_d_phase = zero
	      END IF

          !
          tangent(pos_i1,1) = d_d_bulkEDMat_d_eps_d_phase(1,1) ! dummy
          tangent(pos_i1,2) = d_d_bulkEDMat_d_eps_d_phase(2,2) ! dummy
          tangent(pos_i1,3) = d_d_bulkEDMat_d_eps_d_phase(3,3) ! dummy
          IF (D .GE. 2) THEN
            tangent(pos_i1,4) = d_d_bulkEDMat_d_eps_d_phase(1,2) ! dummy
          END IF
          !IF (D .EQ. 2) THEN
          !  tangent(1:4,pos_i1) = tangent(pos_i1,1:4)
          !ELSE IF (D .EQ. 3) THEN
          IF (D .EQ. 3) THEN
            tangent(pos_i1,5) = d_d_bulkEDMat_d_eps_d_phase(1,3) ! dummy
            tangent(pos_i1,6) = d_d_bulkEDMat_d_eps_d_phase(2,3) ! dummy
            !tangent(1:6,pos_i1) = tangent(pos_i1,1:6)
          END IF
          !
          !
          !
          ! phase-eps: contributions of bulkEnergy
          ! Initialization
          d_d_bulkEDMat_d_phase_d_eps = zero
		  
		  IF (isCouplingterms) THEN
			  IF (isMonolithic) THEN
				!monolithisch
				d_d_bulkEDMat_d_phase_d_eps = zero ! d_d_bulkED_d_eps_d_phase_H(eps,phase(1),nHFEDpar,parHFEDMatrix)
			  ELSEIF (isBFGS) THEN
				!BFGS
				d_d_bulkEDMat_d_phase_d_eps = zero
			  ELSEIF (isStaggered) THEN
				!staggered Martinez-Paneda
				d_d_bulkEDMat_d_phase_d_eps = zero
			  ENDIF 
          ELSE
			  d_d_bulkEDMat_d_phase_d_eps = zero
		  END IF
          
          tangent(1,pos_i1) = d_d_bulkEDMat_d_phase_d_eps(1,1) ! dummy
          tangent(2,pos_i1) = d_d_bulkEDMat_d_phase_d_eps(2,2) ! dummy
          tangent(3,pos_i1) = d_d_bulkEDMat_d_phase_d_eps(3,3) ! dummy
          IF (D .GE. 2) THEN
            tangent(4,pos_i1) = d_d_bulkEDMat_d_phase_d_eps(1,2) ! dummy
          END IF
          !IF (D .EQ. 2) THEN
          !  tangent(1:4,pos_i1) = tangent(pos_i1,1:4)
          !ELSE IF (D .EQ. 3) THEN
          IF (D .EQ. 3) THEN
            tangent(5,pos_i1) = d_d_bulkEDMat_d_phase_d_eps(1,3) ! dummy
            tangent(6,pos_i1) = d_d_bulkEDMat_d_phase_d_eps(2,3) ! dummy
            !tangent(1:6,pos_i1) = tangent(pos_i1,1:6)
          END IF
          !
          !        
          ! eps-grad_phase: contributions of ----
          !
          !
          temp1(:) = zero ! dummy
          tangent(pos_i1+1:pos_i1+D,1) = temp1(1:D)
          temp1(:) = zero ! dummy
          tangent(pos_i1+1:pos_i1+D,2) = temp1(1:D)
          temp1(:) = zero ! dummy
          tangent(pos_i1+1:pos_i1+D,3) = temp1(1:D)
          IF (D .GE. 2) THEN
            temp1(:) = zero ! dummy
            tangent(pos_i1+1:pos_i1+D,4) = temp1(1:D)
          END IF
          IF (D .EQ. 2) THEN
            tangent(1:4,pos_i1+1:pos_i1+D) = TRANSPOSE(tangent(pos_i1+1:pos_i1+D,1:4))
          ELSE IF (D .EQ. 3) THEN
            temp1(:) = zero ! dummy
            tangent(pos_i1+1:pos_i1+D,5) = temp1(1:D)
            temp1(:) = zero ! dummy
            tangent(pos_i1+1:pos_i1+D,6) = temp1(1:D)
            tangent(1:6,pos_i1+1:pos_i1+D) = TRANSPOSE(tangent(pos_i1+1:pos_i1+D,1:6))
          END IF
          !
          DO i2=1,nphase
            pos_i2 = pos_p(i2)
            !
            ! phase-phase: contributions of bulk energy + surface crack energy
            !
            IF (isMonolithic) THEN
              !monolithisch
              !d_d_bulkED_CSED_VDED_d_phase_d_phase = d_d_bulkED_d_phase_d_phase_H(eps,phase(1),nHFEDpar,parHFEDMatrix,H_sph,H_dev) + d_d_CSED_d_phase_d_phase(phase(1),                  grad_phase(1,:),nCSEDpar,parCSEDMatrix) + d_d_VDED_d_phase_d_phase(phase(1),phase_old(1),dtime,nCSEDpar,parCSEDMatrix,nVDEDpar,parVDEDMatrix)
			  d_d_bulkED_CSED_VDED_d_phase_d_phase = d_d_bulkED_d_phase_d_phase_H(eps,phase(1),nHFEDpar,parHFEDMatrix,H_MM) + d_d_CSED_d_phase_d_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSED_MM_Matrix) + d_d_VDED_d_phase_d_phase(phase(1),phase_old(1),dtime,nCSEDpar,parCSED_MM_Matrix,nVDEDpar,parVDEDMatrix)
            ELSEIF (isBFGS) THEN
              !BFGS
              d_d_bulkED_CSED_VDED_d_phase_d_phase = d_d_bulkED_d_phase_d_phase_H(eps,phase(1),nHFEDpar,parHFEDMatrix,H_MM) + d_d_CSED_d_phase_d_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSED_MM_Matrix) + d_d_VDED_d_phase_d_phase(phase(1),phase_old(1),dtime,nCSEDpar,parCSED_MM_Matrix,nVDEDpar,parVDEDMatrix)
            ELSEIF (isStaggered) THEN
              !staggered Martinez-Paneda
              d_d_bulkED_CSED_VDED_d_phase_d_phase = d_d_bulkED_d_phase_d_phase_H(eps_old,phase(1),nHFEDpar,parHFEDMatrix,H_MM) + d_d_CSED_d_phase_d_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSED_MM_Matrix) + d_d_VDED_d_phase_d_phase(phase(1),phase_old(1),dtime,nCSEDpar,parCSED_MM_Matrix,nVDEDpar,parVDEDMatrix)
            ENDIF
            
            tangent(pos_i1,pos_i2) = d_d_bulkED_CSED_VDED_d_phase_d_phase
            !tangent(pos_i1,pos_i2) = zero ! dummy
            !IF (pos_i1 .EQ. pos_i2) THEN 
            !  tangent(pos_i1,pos_i1) = one ! dummy
            !  tangent(pos_i1,pos_i1) = zero ! dummy
            !END IF
            tangent(pos_i2,pos_i1) = tangent(pos_i1,pos_i2)
            !
            
            ! phase-grad_phase: contributions of surface crack energy
            temp1(:) = zero
            temp1 = d_d_CSED_d_phase_d_grad_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSEDMatrix)
            tangent(pos_i1,pos_i2+1:pos_i2+D) = temp1(1:D)
            tangent(pos_i2+1:pos_i2+D,pos_i1) = tangent(pos_i1,pos_i2+1:pos_i2+D)
            !
            
            ! grad_phase-grad_phase: contributions of surface crack energy
            temp2(:,:) = zero
            temp2(:,:) = d_d_CSED_d_grad_phase_d_grad_phase(phase(1),grad_phase(1,:),nCSEDpar,parCSEDMatrix)
            tangent(pos_i1+1:pos_i1+D,pos_i2+1:pos_i2+D) = temp2(1:D,1:D)
            tangent(pos_i2+1:pos_i2+D,pos_i1+1:pos_i1+D) = TRANSPOSE(tangent(pos_i1+1:pos_i1+D,pos_i2+1:pos_i2+D))
            !
          END DO
          !
        END DO
        !
        DEALLOCATE(parCSED_MM_Matrix)
        !
        !
      END FUNCTION tangent

!------------------------------------------------------------------------------------
!
!							SUBROUTINES
!
!------------------------------------------------------------------------------------

SUBROUTINE Check_Stress_Comparison(stress, stran, D, ntens, isSpherisymmetric, &
                                    nphase, pos_p, phase_old, dtime, &
                                    nHFEDpar, nCSEDpar, nVDEDpar, nMixedModepar, parHFEDMatrix, &
                                    parCSEDMatrix, parVDEDMatrix, parMixedModeMatrix, H_sph, H_dev, isMonolithic, isBFGS, isStaggered, isCouplingterms)
  !
  ! Subroutine to check analytical stress against numerical stress
  ! Numerical stress is computed as derivative of total potential
  !
  IMPLICIT NONE

  ! Input/Output variables
  REAL(kind=AbqRK), INTENT(IN) :: stress(:)
  REAL(kind=AbqRK), INTENT(IN) :: stran(:)
  INTEGER(kind=AbqIK), INTENT(IN) :: D
  INTEGER(kind=AbqIK), INTENT(IN) :: ntens
  LOGICAL, INTENT(IN) :: isSpherisymmetric
  INTEGER(kind=AbqIK), INTENT(IN) :: nphase
  INTEGER(kind=AbqIK), INTENT(IN) :: pos_p(:)
  REAL(kind=AbqRK), INTENT(IN) :: phase_old(:)
  REAL(kind=AbqRK), INTENT(IN) :: dtime
  INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar, nCSEDpar, nVDEDpar, nMixedModepar
  REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrix(nHFEDpar)
  REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrix(nCSEDpar)
  REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrix(nVDEDpar)
  REAL(kind=AbqRK), INTENT(IN) :: parMixedModeMatrix(nMixedModepar)
  REAL(kind=AbqRK), INTENT(IN) :: H_sph, H_dev
  LOGICAL, INTENT(IN) :: isMonolithic, isBFGS, isStaggered, isCouplingterms
  
  ! Local variables
  REAL(kind=AbqRK), ALLOCATABLE :: stress_analytical(:), stress_numerical(:)
  REAL(kind=AbqRK), ALLOCATABLE :: stran_per(:)
  REAL(kind=AbqRK), ALLOCATABLE :: eps(:,:), eps_per(:,:), phase(:), phase_per(:), grad_phase(:,:), grad_phase_per(:,:)
  REAL(kind=AbqRK) :: delta, zero
  REAL(kind=AbqRK) :: pot_plus, pot_minus
  INTEGER(kind=AbqIK) :: i1, i2
  LOGICAL :: diff_stress
  
  ! Initialize
  zero = 0.d0
  diff_stress = .FALSE.
  delta = 1.d-7
  
  ! Allocate arrays with proper dimensions
  ALLOCATE(stress_analytical(ntens))
  ALLOCATE(stress_numerical(ntens))
  ALLOCATE(stran_per(ntens))
  ALLOCATE(eps(3,3))
  ALLOCATE(eps_per(3,3))
  ALLOCATE(phase(nphase))
  ALLOCATE(phase_per(nphase))
  ALLOCATE(grad_phase(3,nphase))
  ALLOCATE(grad_phase_per(3,nphase))
  
  !!!!!!!!!!!!
  !
  ! Check Comparison num stress / analytical stress
  !
  !!!!!!!!!!!!
  stress_analytical = stress
          
  ! numerical stress via derivative of potential
  stress_numerical = zero
  DO i1=1,ntens
    ! Forward perturbation
    stran_per = stran
    stran_per(i1) = stran_per(i1) + delta
    eps_per = zero; phase_per = zero; grad_phase_per = zero
    eps_per = strainCoordinates(D,ntens,stran_per,isSpherisymmetric)
    phase_per = phaseParameter(nphase,pos_p,ntens,stran_per)
    grad_phase_per(:,:) = gradientOfPhaseParameter(nphase,D,pos_p,ntens,stran_per)
    pot_plus = potential(eps_per,nphase,phase_per,phase_old,dtime,grad_phase_per,nHFEDpar,nCSEDpar,nVDEDpar,parHFEDMatrix,parCSEDMatrix,parVDEDMatrix)
    
    ! Backward perturbation
    stran_per = stran
    stran_per(i1) = stran_per(i1) - delta
    eps_per = zero; phase_per = zero; grad_phase_per = zero
    eps_per = strainCoordinates(D,ntens,stran_per,isSpherisymmetric)
    phase_per = phaseParameter(nphase,pos_p,ntens,stran_per)
    grad_phase_per(:,:) = gradientOfPhaseParameter(nphase,D,pos_p,ntens,stran_per)
    pot_minus = potential(eps_per,nphase,phase_per,phase_old,dtime,grad_phase_per,nHFEDpar,nCSEDpar,nVDEDpar,parHFEDMatrix,parCSEDMatrix,parVDEDMatrix)
    
    ! Central difference
    stress_numerical(i1) = (pot_plus - pot_minus)/(2.d0*delta)
  END DO
  !
  !
  ! Print Detailed Information
  !
!~   WRITE(6,*) 'Print Test Stress Subroutine'
!~   WRITE(6,*) 'Stress analytical : '
!~   WRITE(6,*) stress_analytical(1)
!~   WRITE(6,*) stress_analytical(2)
!~   WRITE(6,*) stress_analytical(3)
!~   WRITE(6,*) stress_analytical(4)
!~   WRITE(6,*) stress_analytical(5)
!~   WRITE(6,*) stress_analytical(6)
!~   WRITE(6,*) 'Stress numerical : '
!~   WRITE(6,*) stress_numerical(1)
!~   WRITE(6,*) stress_numerical(2)
!~   WRITE(6,*) stress_numerical(3)
!~   WRITE(6,*) stress_numerical(4)
!~   WRITE(6,*) stress_numerical(5)
!~   WRITE(6,*) stress_numerical(6)
!~   !
!~   WRITE(6,*) (Incrementnumber .GT. -1)
  !
  IF (Incrementnumber .GT. 1 .AND. pot_plus*pot_minus .GT. zero) THEN ! Zusätzliche Bedingung um Stop bei num Spannung über Sprung zu verhindern
      DO i1 = 1, ntens
          IF (ABS(stress_analytical(i1) - stress_numerical(i1)) .GT. 0.1) THEN
             diff_stress = .TRUE.
             
             WRITE(6,*) "Increment: ", Incrementnumber
             WRITE(6,*) "Element: ", Elementnumber
             WRITE(6,*) "IP: ", Integrationpointnumber
             WRITE(6,'("Abweichung bei Komponente ",I0,": ",2(E16.8,1X))') i1, stress_analytical(i1), stress_numerical(i1)
             
             ! Lokalisierung
             IF (i1 >= 1 .AND. i1 <= 4) THEN
                WRITE(6,*) '  --> Abweichung in mechanischer Spannung'
             ELSE IF (i1 == 5) THEN
                WRITE(6,*) '  --> Abweichung in Phasenspannung'
             ELSE IF (i1 >= 6 .AND. i1 <= 7) THEN
                WRITE(6,*) '  --> Abweichung in GradPhase Spannung'
             END IF
          END IF
      END DO 
      
      IF (diff_stress) THEN
         ! Get current state for output
         eps = zero; phase = zero; grad_phase = zero
         eps = strainCoordinates(D,ntens,stran,isSpherisymmetric)
         phase = phaseParameter(nphase,pos_p,ntens,stran)
         grad_phase(:,:) = gradientOfPhaseParameter(nphase,D,pos_p,ntens,stran)
         
         WRITE(6,*), "---------------------------------------"
         WRITE(6,*) 'Verzerrung: '
         DO i1 = 1, 3
            WRITE(6,'(3(F10.8,1X))') (eps(i1,i2), i2=1,3)
         END DO
         WRITE(6,*) 'Phasenwert: ', phase(1)				 
         WRITE(6,*), "---------------------------------------"
         WRITE(6,*) 'Stress analytisch:'
         DO i1 = 1, ntens
            WRITE(6,'(E16.8)') stress_analytical(i1)
         END DO
         WRITE(6,*) ''
         WRITE(6,*) 'Stress numerisch:'
         DO i1 = 1, ntens
            WRITE(6,'(E16.8)') stress_numerical(i1)
         END DO
         WRITE(6,*), "---------------------------------------"
         CALL XEXIT()
      END IF
  END IF
  
  ! Deallocate arrays
  DEALLOCATE(stress_analytical, stress_numerical, stran_per, eps, eps_per, phase, phase_per, grad_phase, grad_phase_per)
  
END SUBROUTINE Check_Stress_Comparison

!------------------------------------------------------------------------------------

SUBROUTINE Check_Tangent_Comparison(Ct, stress, stran, D, ntens, isSpherisymmetric, &
                                     nphase, pos_p, prop_solver, phase_old, dtime, &
                                     nHFEDpar, nCSEDpar, nVDEDpar, nMixedModepar, parHFEDMatrix, &
                                     parCSEDMatrix, parVDEDMatrix, parMixedModeMatrix, H_sph, H_dev, isMonolithic, isBFGS, isStaggered, isCouplingterms)
  !
  ! Subroutine to check analytical tangent against numerical tangent
  !
  IMPLICIT NONE
  !
  ! Input/Output variables
  REAL(kind=AbqRK), INTENT(IN) :: Ct(:,:)
  REAL(kind=AbqRK), INTENT(IN) :: stress(:)
  REAL(kind=AbqRK), INTENT(IN) :: stran(:)
  INTEGER(kind=AbqIK), INTENT(IN) :: D
  INTEGER(kind=AbqIK), INTENT(IN) :: ntens
  LOGICAL, INTENT(IN) :: isSpherisymmetric
  INTEGER(kind=AbqIK), INTENT(IN) :: nphase
  INTEGER(kind=AbqIK), INTENT(IN) :: pos_p(:)
  INTEGER(kind=AbqIK), INTENT(IN) :: prop_solver
  REAL(kind=AbqRK), INTENT(IN) :: phase_old(nphase)
  REAL(kind=AbqRK), INTENT(IN) :: dtime
  INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar, nCSEDpar, nVDEDpar, nMixedModepar
  REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrix(nHFEDpar)
  REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrix(nCSEDpar)
  REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrix(nVDEDpar)
  REAL(kind=AbqRK), INTENT(IN) :: parMixedModeMatrix(nMixedModepar)
  REAL(kind=AbqRK), INTENT(IN) :: H_sph, H_dev
  LOGICAL, INTENT(IN) :: isMonolithic, isBFGS, isStaggered, isCouplingterms
  
  ! Local variables
  REAL(kind=AbqRK), ALLOCATABLE :: Ct_temp1(:,:), Ct_temp2(:,:)
  REAL(kind=AbqRK), ALLOCATABLE :: stran_per(:), stress_per(:)
  REAL(kind=AbqRK), ALLOCATABLE :: eps(:,:), eps_old(:,:), eps_per(:,:), phase(:), grad_phase(:,:)
  REAL(kind=AbqRK) :: delta, zero, one, trace_eps, trace_eps_per
  INTEGER(kind=AbqIK) :: i1, i2
  LOGICAL :: diff_tangent, valid_traces

  ! Initialize
  zero = 0.d0
  one = 1.d0
  diff_tangent = .FALSE.
  valid_traces = .TRUE.
  delta = 1.d-7

  ! Allocate arrays
  ALLOCATE(Ct_temp1(ntens,ntens))
  ALLOCATE(Ct_temp2(ntens,ntens))
  ALLOCATE(stran_per(ntens))
  ALLOCATE(stress_per(ntens))
  ALLOCATE(eps_old(3,3))
  ALLOCATE(eps(3,3))
  ALLOCATE(eps_per(3,3))
  ALLOCATE(phase(nphase))
  ALLOCATE(grad_phase(3,nphase))

  !!!!!!!!!!!!
  !
  ! Check Comparison num tangent / analytical tangent
  !
  !!!!!!!!!!!!
  Ct_temp1 = Ct

  ! Calculate current eps and its trace
  eps = zero
  eps = strainCoordinates(D,ntens,stran,isSpherisymmetric)
  trace_eps = eps(1,1) + eps(2,2) + eps(3,3)
	  
  ! numerical tangent
  Ct_temp2 = zero
  DO i1=1,ntens
	stran_per = stran
	stran_per(i1) = stran_per(i1) + delta
	eps_per = zero; phase = zero; grad_phase = zero
	eps_per = strainCoordinates(D,ntens,stran_per,isSpherisymmetric)
	trace_eps_per = eps_per(1,1) + eps_per(2,2) + eps_per(3,3)
  
	! Check if traces have same sign
	IF (ABS(trace_eps) .LT. 1.d-8 .OR. ABS(trace_eps_per) .LT. 1.d-8 .OR. sign(one,trace_eps) .NE. sign(one,trace_eps_per)) THEN
	  valid_traces = .FALSE.
	END IF
  
	phase = phaseParameter(nphase,pos_p,ntens,stran_per)
	grad_phase(:,:) = gradientOfPhaseParameter(nphase,D,pos_p,ntens,stran_per)
	stress_per = stresses(D,ntens,nphase,pos_p,prop_solver,eps_per, eps_old, phase,phase_old,dtime,grad_phase,nHFEDpar,nCSEDpar,nVDEDpar,nMixedModepar,parHFEDMatrix,parCSEDMatrix,parVDEDMatrix,parMixedModeMatrix,H_sph,H_dev,isMonolithic, isBFGS, isStaggered, isCouplingterms)
	Ct_temp2(i1,1:ntens) = (stress_per(1:ntens)-stress(1:ntens))/delta
  END DO
  !
  ! Print Detailed Information
  !
  IF (Incrementnumber .GT. 1 .AND. valid_traces) THEN
      DO i1 = 1, ntens
           DO i2 = 1, ntens
              IF (ABS(Ct_temp1(i1,i2) - Ct_temp2(i1,i2)) .GT. 1.) THEN
                 diff_tangent = .TRUE.
                 
                 WRITE(6,*) "Increment: ", Incrementnumber
                 WRITE(6,*) "Element: ", Elementnumber
                 WRITE(6,*) "IP: ", Integrationpointnumber
                 WRITE(6,'("Abweichung bei (",I0,",",I0,"): ",2(F10.4,1X))') i1, i2, Ct_temp1(i1,i2), Ct_temp2(i1,i2)
                 
                 ! Lokalisierung fuer 2D
                 ! Bereich 1: Ct_temp(1,1) - Ct_temp(4,4)
                 IF (i1 >= 1 .AND. i2 >= 1 .AND. i1 <= 4 .AND. i2 <= 4) THEN
                    WRITE(6,*) '  --> Abweichung in Elast. Tangente'
                 ! Bereich 2: Ct_temp(5,5) - Ct_temp(5,5)
                 ELSE IF (i1 >= 5 .AND. i2 >= 5 .AND. i1 <= 5 .AND. i2 <= 5) THEN
                    WRITE(6,*) '  --> Abweichung in Phase Tangente'
                    
                 ! Bereich 3: Ct_temp(1,5) - Ct_temp(4,5)
                 ELSE IF (i1 >= 1 .AND. i2 >= 5 .AND. i1 <= 4 .AND. i2 <= 5) THEN
                    WRITE(6,*) '  --> Abweichung in phase-eps Tangente'
                    
                 ! Bereich 4: Ct_temp(5,1) - Ct_temp(5,4)
                 ELSE IF (i1 >= 5 .AND. i2 >= 1 .AND. i1 <= 5 .AND. i2 <= 4) THEN
                    WRITE(6,*) '  --> Abweichung in eps-phase Tangente'
                    
                 ! Bereich 7: Ct_temp(6,6) - Ct_temp(7,7)
                 ELSE IF (i1 >= 6 .AND. i2 >= 6 .AND. i1 <= 7 .AND. i2 <= 7) THEN
                    WRITE(6,*) '  --> Abweichung in GradPhase Tangente'
                 END IF
              END IF
           END DO
        END DO 
      IF (diff_tangent) THEN
             WRITE(6,*), "---------------------------------------"
             WRITE(6,*) 'Verzerrung: '
             DO i1 = 1, 3
                WRITE(6,'(3(F10.8,1X))') (eps(i1,i2), i2=1,3)
             END DO
             WRITE(6,*) 'Phasenwert: ', phase(1)				 
             WRITE(6,*), "---------------------------------------"
             WRITE(6,*) 'Ct analytisch:'
             DO i1 = 1, ntens
                WRITE(6,'(7(F16.4,1X))') (Ct_temp1(i1,i2), i2=1,ntens)
             END DO
             WRITE(6,*) ''
             WRITE(6,*) 'Ct numerisch:'
             DO i1 = 1, ntens
                WRITE(6,'(7(F16.4,1X))') (Ct_temp2(i1,i2), i2=1,ntens)
             END DO
             WRITE(6,*), "---------------------------------------"
        CALL XEXIT()
      END IF
  END IF
  
  ! Deallocate arrays
  DEALLOCATE(Ct_temp1, Ct_temp2, stran_per, stress_per, eps, phase, grad_phase)
  
END SUBROUTINE Check_Tangent_Comparison

END MODULE PhaseField_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
