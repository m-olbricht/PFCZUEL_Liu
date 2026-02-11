!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Stephan Roth, TU Bergakademie Freiberg, 33.01.2024
!
! 06.01.2016: Modifikationen fuer UELlib-Kompatibilitaet
! 25.01.2024: Strafenergie um Sprung in Schaedigungsvariable zu minimieren,
!             konstante Kontaktsteifigkeit (linear),
!             Grenzflaechenenergie
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SAVED VARIABLES AND INPUT PARAMETERS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SVR: svr(1)     : damage variable 0<=D<=1
!      svr(2)     : 
!      svr(3)     : damage jump
!      svr(4)     : -
!      svr(5)     : total potential density
!      svr(6)     : stored energy density
!      svr(7)     : interface energy density
!      svr(8)     : contact energy density
!      svr(9)     : penalty energy density
!      svr(10-12) : separation (local frame)
!      svr(13-15) : separation (global frame)
!      svr(16-18) : damage gradient (local frame)
!      svr(19-21) : damage gradient (global frame)
!      svr(22-24) : traction vector (local frame)
!      svr(25-30) : traction vector (global frame)
!	   svr(31)	  : degradation Function interface
!	   svr(32)	  : viscous energy density
!	   svr(33)	  : Hi_sphn (alter history Wert für sphaerischen separationsanteil)
!	   svr(34)	  : Hi_devn (alter history Wert für deviatorischen separationsanteil)
!	   svr(35)	  :
!
! material parameters: thicknessIndex -- plane strain thickness
!                       1 -- internal length l_i
!                       2 -- normal stiffness kn
!                       3 -- tangential stiffness kt
!                       4 -- flag: numerical tangent (inactive)
!                       5 -- penalty factor
!                       6 -- critical interface fracture energy density
!                       7 -- plane strain thickness 
!                       8 -- factor to scale stiffness for negative normal separation
!                       9 -- parameter 1 of degradation function
!                      10 -- parameter 2 of degradation function
!                      11 -- viscous dissipation numerical parameter
!                      12 -- GiI
!                      13 -- GiII
!                      14 -- 
!                      15 -- 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUMAT_module

  USE SharedValues

  IMPLICIT NONE

  PUBLIC :: CheckMaterialParameters, umatCZM

  CONTAINS

!------------------------------------------------------------------------------------

    SUBROUTINE CheckMaterialParameters(props)
      ! Check of all material parameters

      USE ABQINTERFACE
      USE FLOATNUMBERS
      USE ABAModul
      USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: props(numMatPar)
      INTEGER(kind=AbqIK) :: num_counter
      REAL(kind=AbqRK) :: prop_kn, prop_kt, prop_thickness

      prop_thickness = one

      ! read material parameters from props
      prop_kn      = props( 2) ! stiffness normal
	  prop_kt	   = props( 3) ! stiffness tangential
      IF (props(thicknessIndex) .NE. zero) prop_thickness = props(thicknessIndex) ! thickness of 2D plane strain element

      ! parameter check
      IF (prop_kn .LE. zero) THEN
        WRITE(7,*) 'Stiffness kn should exceed zero. EXIT', prop_kn
        CALL XEXIT()
      END IF
	  IF (prop_kt .LE. zero) THEN
        WRITE(7,*) 'Stiffness kt should exceed zero. EXIT', prop_kt
        CALL XEXIT()
      END IF
      IF (prop_thickness .LT. zero) THEN
        WRITE(7,*) 'Thickness of 2D element should exceed zero. EXIT', prop_thickness
        CALL XEXIT()
      END IF
!~       IF (prop_s0 .LT. zero) THEN
!~         WRITE(7,*) 'cohesive length must exceed zero. EXIT', prop_thickness
!~         CALL XEXIT()
!~       END IF
!      IF (prop_thickness .LT. zero) THEN
!        WRITE(7,*) 'Thickness of 2D element should exceed zero. EXIT', prop_thickness
!        CALL XEXIT()
!      END IF

    END SUBROUTINE CheckMaterialParameters

!------------------------------------------------------------------------------------

    SUBROUTINE umatCZM(stress,svr,Ct,energy_elast,dissipat_plast,dissipat_creep,rpl,ddsddt,drplde,drpldt, &
                       stran,dstran,time,dtime,Temp,dTemp,predef_umat,dpred,cmname,ndi,nshr,ntens, &
                       num_svr,props_mat,num_matpar,coords_gp,Trafomat,pnewdt,intlen,F0,F1,jelem, &
                       npt,layer,kspt,kstep,kinc)
      ! compute cohesive tractions and cohesive stiffness

      USE ABQINTERFACE
      USE FLOATNUMBERS
      USE ABAModul
      USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
      USE CohesiveEnergyModule
      USE SplitIntEnergyModule
      USE ContactEnergyModule
	  USE AliasModuleCZ



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

	  ! pure debug Variables
      INTEGER(kind=AbqIK), PARAMETER :: n_debug_elem = 1
      INTEGER(kind=AbqIK), PARAMETER :: n_debug_inc = 5
      INTEGER(kind=AbqIK), PARAMETER :: n_debug_ip = 1
      
      
      INTEGER(kind=AbqIK) :: debug_elements(n_debug_elem)
      INTEGER(kind=AbqIK) :: debug_increments(n_debug_inc)
      INTEGER(kind=AbqIK) :: debug_ipoints(n_debug_ip)
      
      ! Debug Werte
      DATA debug_elements /2/
      DATA debug_increments /1,2,3,4,5/
      DATA debug_ipoints /1/


      ! further variables
      INTEGER(kind=AbqIK) :: D, i1, i2, pos_damage, pos_damageJump, pos_damageGradient, nIEDpar, nCEDpar, nVDEDpar, nMixedModeCZpar, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK) :: damage, damage_old, damageJump
      REAL(kind=AbqRK) :: lambda
!~       REAL(kind=AbqRK) :: Hin, Hi, Himin
	  REAL(kind=AbqRK) :: Hi_sph, Hi_dev, Hi_sphn, Hi_devn
	  REAL(kind=AbqRK) :: Gi_MM, Hi_MM
      REAL(kind=AbqRK) :: degFi ! Nur für svr Deg Funktion
      REAL(kind=AbqRK) :: prop_l0, prop_kn, prop_kt, prop_G0, prop_pen, prop_compressionstiffness, prop_df_one, prop_df_two, prop_visCZ, prop_GiI, prop_GiII
      REAL(kind=AbqRK) :: totalPotential, cohesiveEnergy, cohesive_energy_sph, cohesive_energy_dev, interfaceEnergy, penaltyEnergy, contactEnergy, viscousEnergy, transformedStress(6)
      REAL(kind=AbqRK), ALLOCATABLE :: sep(:), sepGlobal(:), trac(:), damageGradient(:), damageGradientLocal(:), damageGradientGlobal(:), &
                                       K(:,:), tmp(:), tmp2(:,:), parIEDMatrix(:), parIED_MM_Matrix(:), parCEDMatrix(:), parVDEDMatrix(:), parMixedModeCZMatrix(:)
      LOGICAL :: NAN_Check, debug_Flag, is_print_inc, is_print_elem, is_print_ip

      ! initialization
      ! add, if needed
      NAN_Check = .TRUE.
      
      !debug Flag
      debug_Flag = .TRUE.
      
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
	  END IF

      ! read material parameters from props_mat
      prop_l0      = props_mat( 1) ! internal length
      prop_kn      = props_mat( 2) ! stiffness n
      prop_kt	   = props_mat( 3) ! stiffness t
      prop_pen     = props_mat( 5) ! penalty factor
      prop_G0      = props_mat( 6) ! critical interface fracture energy density
      prop_compressionstiffness = props_mat( 8) ! penalty factor for negative normal separation
      prop_df_one  = props_mat( 9) ! parameter 1 of degradation function
      prop_df_two  = props_mat(10) ! parameter 2 of degradation function
      prop_visCZ   = props_mat(11) ! parameter viscous dissipation
      prop_GiI	   = props_mat(12) ! GiI
      prop_GiII	   = props_mat(13) ! GiII
	  ! arrays of material parameters
      !
      ! interface energy density -- see InterfaceEnergyModule.f90
      nIEDpar = 2
      ALLOCATE(parIEDMatrix(nIEDpar))
      parIEDMatrix(1) = prop_G0
      parIEDMatrix(2) = prop_l0
      ! interface energy with new effective Gi (updated in the code)
      ALLOCATE(parIED_MM_Matrix(nIEDpar))
      parIED_MM_Matrix(1) = zero
      parIED_MM_Matrix(2) = zero
      ! contact energy density -- see ContactEnergyModule.f90
      nCEDpar = 3
      ALLOCATE(parCEDMatrix(nCEDpar))
      parCEDMatrix(1) = prop_kn
      parCEDMatrix(2) = prop_kt
      parCEDMatrix(3) = prop_compressionstiffness
      ! viscous dissipation energy density -- see ViscousDissipationEnergyModule.f90
      nVDEDpar = 1
      ALLOCATE(parVDEDMatrix(nVDEDpar))
      parVDEDMatrix(1) = prop_visCZ
      ! MixedMode Parameter
      nMixedModeCZpar = 2
      ALLOCATE(parMixedModeCZMatrix(nMixedModeCZpar))
      parMixedModeCZMatrix(1) = prop_GiI
      parMixedModeCZMatrix(2) = prop_GiII
      
	  
!~ 	  IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
!~ 		  WRITE(6,*) '___________________________________'
!~ 		  WRITE(6,*) '======IED======='
!~ 		  WRITE(6,*) 'parIEDMatrix(1)', parIEDMatrix(1)
!~ 		  WRITE(6,*) 'parIEDMatrix(2)', parIEDMatrix(2)
!~ 		  WRITE(6,*) '======CED======='
!~ 		  WRITE(6,*) 'parCEDMatrix(1)', parCEDMatrix(1)
!~ 		  WRITE(6,*) 'parCEDMatrix(2)', parCEDMatrix(2)
!~ 		  WRITE(6,*) 'parCEDMatrix(3)', parCEDMatrix(3)
!~ 		  WRITE(6,*) '======VDED======='
!~ 		  WRITE(6,*) 'parVDEDMatrix(1)', parVDEDMatrix(1)
!~ 		  WRITE(6,*) '=====Mixed======='
!~ 		  WRITE(6,*) 'parMixedModeCZMatrix(1)', parMixedModeCZMatrix(1)
!~ 		  WRITE(6,*) 'parMixedModeCZMatrix(2)', parMixedModeCZMatrix(2)
!~ 		  WRITE(6,*) '_________________________________________'
!~ 	  END IF
      ! model dimension
      D = ndi+nshr
      ! index of normal direction
      normalDirectionIndex = D
      tangentialDirectionIndex = D-1

      ALLOCATE(sep(D), sepGlobal(D), trac(D), K(D,D), tmp(D), tmp2(D,D), damageGradient(D-1), damageGradientLocal(D), damageGradientGlobal(D))
      sep=zero; sepGlobal=zero; trac=zero; K=zero; tmp=zero; tmp2=zero; damageGradient=zero; damageGradientLocal=zero; damageGradientGlobal=zero

      ! position in generalized arrays
      pos_damage         = ntens-(D+1)+1
      pos_damageJump     = ntens-(D+1)+2
      pos_damageGradient = ntens-(D+1)+3


      ! generalized kinematic measures
      !
      ! separation vector coordinates
      sep(1:D) = stran(1:D)
      ! damage variable
      damage = stran(pos_damage)
      ! damage jump
      damageJump = stran(pos_damageJump)
      ! damage gradient
      damageGradient = stran(pos_damageGradient:pos_damageGradient-2+D)
      ! damage old
      damage_old = svr(1)
!~       ! normalised effective separation
!~       lambda = effSep(D,normalDirectionIndex,sep,prop_delta0)
      ! damage gradient array of size D
      damageGradientLocal(1:D-1) = damageGradient(1:D-1)
	  !
	  !
	  !
	  IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
	   WRITE(6,*) "---------------------------------------"
	   WRITE(6,*) " Debug-Ausgabe:"
	   WRITE(6,*) " Increment =", Incrementnumber
	   WRITE(6,*) " Element   =", Elementnumber
	   WRITE(6,*) " IP        =", Integrationpointnumber
	   WRITE(6,*) "---------------------------------------"
	  END IF
	  !
	  ! Historyvariable
	  !
	  cohesive_energy_sph = CZEDtens_H_sph(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrix)
	  cohesive_energy_dev = CZEDtens_H_dev(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrix)
	  !
	  !
	  Hi_sphn = svr(33)
	  Hi_devn = svr(34)
	  !
	  Hi_sph = zero
	  Hi_dev = zero
	  !	  
	  !
	  IF (cohesive_energy_sph .GT. Hi_sphn) THEN
	    Hi_sph = cohesive_energy_sph
	  ELSE
		Hi_sph = Hi_sphn
	  ENDIF
	  
	  IF (cohesive_energy_dev .GT. Hi_devn) THEN
	    Hi_dev = cohesive_energy_dev
	  ELSE
		Hi_dev = Hi_devn
	  ENDIF
	  
	  svr(33) = Hi_sph
	  svr(34) = Hi_dev
	  
	  !
!~ 	  Himin = three*prop_G0/(16*prop_l0)
!~ 	  IF (Hi .LT. Himin) THEN
!~ 	    Hi = Himin
!~ 	  END IF 
	  !
	  !
!~ 	  IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
!~ 		WRITE(6,*) 'Hi: ', Hi
!~ 	  END IF
		
		
      ! energetic density quantities
      !
      cohesiveEnergy  = CZED(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrix,damage,prop_df_one,prop_df_two)
      interfaceEnergy = IED(D-1,damage,damageGradient,nIEDpar,parIEDMatrix) 
      penaltyEnergy   = PenED(damageJump,prop_pen)
      contactEnergy   = ContactED(D,normalDirectionIndex,sep,nCEDpar,parCEDMatrix)
      viscousEnergy   = VDED(damage,damage_old,dtime,nIEDpar,parIEDMatrix,nVDEDpar,parVDEDMatrix)
      !
      totalPotential  = cohesiveEnergy + interfaceEnergy + penaltyEnergy + contactEnergy + viscousEnergy
      !
!~       IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
!~         WRITE(6,*) 'OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO'
!~ 		WRITE(6,*) 'viscous Energy: ', viscousEnergy
!~ 		WRITE(6,*) ''
!~ 		WRITE(6,*) 'prop visnum: ', prop_visCZ
!~ 		WRITE(6,*) 'damage: ', damage
!~ 		WRITE(6,*) 'damage_old: ', damage_old
!~ 		WRITE(6,*) 'dtime: ', dtime
!~ 		WRITE(6,*) 'OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO'
!~ 	  END IF
      !
      ! assignments
      !energy_elast = penalty Energy (further processed in CZUEL, not added to bulkEnergy)
      energy_elast   = penaltyEnergy
      !dissipat_plast = total Potential CZ
      dissipat_plast = totalPotential
      !dissipat_creep = crack surface energy + interface Energy(CZ Variant of CSE)
      dissipat_creep = interfaceEnergy		
		
!~ 	  IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
!~ 		WRITE(6,*) 'cohesive Energy density: ', cohesiveEnergy
!~ 		WRITE(6,*) ''
!~ 		WRITE(6,*) '---------------------------'
!~ 	  END IF
      ! compute generalized stresses
      !
      stress = zero
      !
      ! stress(1) = T_t (Tangentialtraction)
      ! stress(2) = T_n (Normaltraction)
      ! stress(3) = T_phi (Driving Force Phase / ONLY for 2D)
      ! stress(4) = T_damage_Jump (Damage Jump / ONLY for 2D)
      ! stress(5) = T_phi_grad (Damage Gradient / ONLY for 2D)
      !
      !
      ! local tractions
      trac = d_CZED_d_sep(sep,D,normalDirectionIndex,tangentialDirectionIndex,nCEDpar,parCEDMatrix,damage,prop_df_one,prop_df_two)
      ! consider negative normal separation / compression
!~       trac = tracWithContact(D,normalDirectionIndex,trac,sep,nCEDpar,parCEDMatrix)
      ! arrange in stress array
      stress(1:D) = trac(1:D)
      !
      ! damage conjugate -- damage driving force
      ! contribution from elastic free energy density and surface energy density
      !
      ! For Stress Function to change Gc
	  parIED_MM_Matrix(:) = parIEDMatrix(:)
      !
	  ! Mixed Mode          
	  !
	  Gi_MM = Gieff(Hi_sph, Hi_dev, nMixedModeCZpar,parMixedModeCZMatrix)
	  Hi_MM = Hieff(Hi_sph, Hi_dev, nMixedModeCZpar,parMixedModeCZMatrix)
	  !
	  IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
	  		WRITE(6,*)
			WRITE(6,*) '=============================================='
			WRITE(6,*) ' DEBUG: Separation Damage'
			WRITE(6,*) '=============================================='
			WRITE(6,*) 'separation normal = ', sep(normalDirectionIndex)
			WRITE(6,*) 'separation tangential = ', sep(tangentialDirectionIndex)
			WRITE(6,*) 'damage: ', damage

	  
			WRITE(6,*)
			WRITE(6,*) '=============================================='
			WRITE(6,*) ' DEBUG: parIED_MM_Matrix BEFORE Mixed Mode'
			WRITE(6,*) '=============================================='
			WRITE(6,*) 'parIED_MM_Matrix(1) = ', parIED_MM_Matrix(1)
			WRITE(6,*) 'parIED_MM_Matrix(2) = ', parIED_MM_Matrix(2)

			WRITE(6,*)
			WRITE(6,*) '----------------------------------------------'
			WRITE(6,*) ' DEBUG: Arguments for Gieff / Hieff'
			WRITE(6,*) '----------------------------------------------'
			WRITE(6,*) 'Hi_sph            = ', Hi_sph
			WRITE(6,*) 'Hi_dev            = ', Hi_dev
			WRITE(6,*) 'nMixedModepar    = ', nMixedModeCZpar
			WRITE(6,*) 'parMixedModeMatrix = ', parMixedModeCZMatrix(:)

			WRITE(6,*)
			WRITE(6,*) '----------------------------------------------'
			WRITE(6,*) ' DEBUG: Results from Gieff / Hieff'
			WRITE(6,*) '----------------------------------------------'
			WRITE(6,*) 'Gi_MM = ', Gi_MM
			WRITE(6,*) 'Hi_MM  = ', Hi_MM

			WRITE(6,*)
			WRITE(6,*) '----------------------------------------------'
			WRITE(6,*) ' DEBUG: Updating parIED_MM_Matrix(2)'
			WRITE(6,*) '----------------------------------------------'
			WRITE(6,*) 'Old parIED_MM_Matrix(2) = ', parIED_MM_Matrix(2)
	  END IF	
	  !
	  parIED_MM_Matrix(1) = Gi_MM
      !
      stress(pos_damage) = d_CZED_d_damage_Hi(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrix,damage,prop_df_one,prop_df_two, Hi_MM) &
                         + d_IED_d_phase(D-1,damage,damageGradient,nIEDpar,parIED_MM_Matrix) &
                         + d_VDED_d_phase(damage,damage_old,dtime,nIEDpar,parIED_MM_Matrix,nVDEDpar,parVDEDMatrix)


      !DEBUGDEBUG		
!~ 	  IF (Incrementnumber .GE. 1 .AND. Elementnumber .EQ. 1 .AND. Integrationpointnumber .EQ. 1) THEN
!~ 	    WRITE(7,*) "========================================"
!~ 	    WRITE(7,*) "Incrementnumber", Incrementnumber
!~ 	    WRITE(7,*) "Elementnumber", Elementnumber
!~ 	    WRITE(7,*) "IPTnumber", Integrationpointnumber
!~ 	    WRITE(7,*) "========================================"
!~ 	    WRITE(7,*) "sep: ", sep(:)
!~ 	    WRITE(7,*) "phase ", damage
!~ 	    WRITE(7,*) "trac: ", trac(:)
!~ 	    WRITE(7,*) "T0: ", prop_t0
!~ 	    WRITE(7,*) "s0: ", prop_delta0
!~ 	    WRITE(7,*) "CZED: ", cohesiveEnergy
!~ 	    WRITE(7,*) "driving_Force_d_damage: ", stress(pos_damage)
!~ 	    WRITE(7,*) "d_CZED_d_damage: ", dbg_A
!~ 	    WRITE(7,*) "d_INTERF_d_damage: ", dbg_B
!~ 	    WRITE(7,*) "========================================"
!~ 	  END IF
      
      
      
      !
      ! damage jump conjugate
      ! contribution from penalty energy density: reduce damage jump to a minimum (approx. zero)
      stress(pos_damageJump) = d_PenED_d_damageJump(damageJump,prop_pen)
      !
      ! damage gradient conjugate
      ! contribution from surface energy density
      stress(pos_damageGradient:pos_damageGradient+D-2) = d_IED_d_grad_phase(D-1,damage,damageGradient,nIEDpar,parIED_MM_Matrix)

	
      ! generalised tangent
      !
      Ct = zero
      !
      ! separation
      ! contribution from elastic free energy density
      K = d_CZED_d_sep_d_sep(sep,D,normalDirectionIndex,tangentialDirectionIndex,nCEDpar,parCEDMatrix,damage,prop_df_one,prop_df_two)
      ! consider negative normal separation / compression
!~       K = tangentWithContact(D,normalDirectionIndex,K,sep,nCEDpar,parCEDMatrix)
      ! arrange  in tangent array
      Ct(1:D,1:D) = K(1:D,1:D)
      !
      ! damage
      ! contribution from elastic free energy density and surface energy density
      Ct(pos_damage,pos_damage) = d_CZED_d_damage_d_damage_Hi(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrix,damage,prop_df_one,prop_df_two, Hi_MM) &
                                + d_IED_d_phase_d_phase(D-1,damage,damageGradient,nIEDpar,parIED_MM_Matrix) &
                                + d_VDED_d_phase_d_phase(damage,damage_old,dtime,nIEDpar,parIED_MM_Matrix,nVDEDpar,parVDEDMatrix)
      !
      ! mixed submatrix
      ! contribution from elastic free energy density and surface energy density
      tmp(1:D) = d_CZED_d_sep_d_damage(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrix,damage,prop_df_one,prop_df_two)
      Ct(1:D,pos_damage) = tmp(1:D)
      ! d_CZED_d_damage_d_sep Beitrag ist 0
      Ct(pos_damage,1:D) = zero ! d_CZED_d_damage_d_sep -> 0 History Feld
      !
      ! damage gradient
      ! contribution from surface energy density
      tmp2(1:D-1,1:D-1) = d_IED_d_grad_phase_d_grad_phase(D-1,damage,damageGradient,nIEDpar,parIED_MM_Matrix)
      Ct(pos_damageGradient:pos_damageGradient+D-2,pos_damageGradient:pos_damageGradient+D-2) = tmp2(1:D-1,1:D-1)
      !
      ! damage jump
      ! contribution from penalty energy: reduce damage jump to a minimum (approx. zero)
      Ct(pos_damageJump,pos_damageJump) = d_PenED_d_damageJump_d_damageJump(damageJump,prop_pen)
      !
      ! no coupling between damage and damage gradient
      ! no coupling between damage jump and separation, damage, or damage gradient
      !
      
      
      
      ! ====================================================================================================================
	  !
	  ! SUBROUTINE CHECKS
	  !
	  ! ====================================================================================================================
	  ! Stress deviations
	  IF (is_print_elem .AND. is_print_inc .AND. is_print_ip) THEN
		  CALL Check_Stress_Comparison_CZ(stress, stran, D, ntens, normalDirectionIndex, &
										 tangentialDirectionIndex, prop_df_one, prop_df_two, &
										 Hi_MM, prop_pen, damage_old, dtime, nCEDpar, parCEDMatrix, &
										 nIEDpar, parIED_MM_Matrix, nVDEDpar, parVDEDMatrix)
	  END IF
	  
	  ! NAN in Ct                     
      IF (Incrementnumber .GT. zero) THEN    
		  CALL check_matrix_entries(Ct, D, pos_damage, pos_damageGradient, pos_damageJump, &
									sep, normalDirectionIndex, tangentialDirectionIndex, &
									nCEDpar, parCEDMatrix, damage, prop_df_one, prop_df_two, &
									Hi_MM, damageGradient, damage_old, dtime, nIEDpar, &
									parIED_MM_Matrix, nVDEDpar, parVDEDMatrix, damageJump, prop_pen)
	  END IF
	  ! ====================================================================================================================



      ! Transformations (change of vector basis)
      !
      ! separation in global frame, Q^T*s
      sepGlobal            = TransformedVectorLocalToGlobal(D,sep,TrafoMat)
      ! damage gradient in global frame, Q^T*grad_d
      damageGradientGlobal = TransformedVectorLocalToGlobal(D,damageGradientLocal,TrafoMat)
      ! traction in global frame, Q^T*sigma*Q (stress tensor in Abaqus notation 3D: 11-22-33-12-13-23, 2D: 11-22-33-12)
      transformedStress    = TransformedTractionVectorLocalToStressArrayGlobal(D,trac,TrafoMat)
      
      ! 
      ! Degradationsfunktion
      degFi = degF(prop_df_one, prop_df_two, damage)


      ! save current state variables
      !
      svr(1)  = damage
      svr(2)  = 0.
      svr(3)  = damageJump
      !
      !
      ! energy densities (per surface area)
      svr(5) = totalPotential
      svr(6) = cohesiveEnergy
      svr(7) = interfaceEnergy
      svr(8) = contactEnergy
      svr(9) = penaltyEnergy
      svr(32)= viscousEnergy
      !
      ! separation vector
      ! local frame
      svr(10:9+D)  = sep
      ! global frame
      svr(13:12+D) = sepGlobal
      !
      ! damage gradient vector
      ! local frame
      svr(16:15+D) = damageGradientLocal
      ! global frame
      svr(19:18+D) = damageGradientGlobal
      !
      ! cohesive traction vector
      ! local frame
      svr(22:21+D) = trac(1:D)
      ! global frame (stress tensor in Abaqus notation 3D: 11-22-33-12-13-23, 2D: 11-22-33-12)
      svr(25:30) = transformedStress
      !
      ! Degradationsfunktion
	  svr(31) = degFi
	  
      ! dummy output
      rpl = zero; ddsddt = zero; drplde = zero; drpldt = zero
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! PRINTOUTS
      !
	  !WRITE(7,*) 'INCREMENT: ', kinc
	  !WRITE(7,*) 'ELEMENT: ', jelem
	  !WRITE(7,*) 'time: ', time
      !WRITE(7,*) 'traction: ', svr(21+D)
	  !WRITE(7,*) 'sep: ', svr(9+D)
	  !WRITE(7,*) 'damage: ', svr(1)
	  !WRITE(7,*) 'drivingforce_stored', stress_A
	  !WRITE(7,*) 'drivingforce_crack', stress_B
	  !WRITE(7,*) 'damage driving force ', stress(pos_damage)
	  !WRITE(7,*) 'G0: ', props_mat( 6)
	  !WRITE(7,*) 'Hi: ', svr(4)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !
	  ! NAN CHECKS
	  !
!~ 	  IF (NAN_Check) THEN
!~ 		  IF (cohesiveEnergy .NE. cohesiveEnergy) THEN
!~ 			WRITE(7,*) 'cohesive Energy NAN: ', cohesiveEnergy
!~ 			WRITE(7,*) 'INCREMENT: ', kinc
!~ 			WRITE(7,*) 'ELEMENT: ', jelem
!~ 			!
!~ 			! Function Variables
!~ 			!
!~ 			WRITE(7,*) 'separation: ', sep
!~ 			WRITE(7,*) 'separation: ', sep
!~ 		  END IF
!~ 		  IF (interfaceEnergy .NE. interfaceEnergy) THEN
!~ 			WRITE(7,*) 'interface Energy NAN: ', interfaceEnergy
!~ 			WRITE(7,*) 'INCREMENT: ', kinc
!~ 			WRITE(7,*) 'ELEMENT: ', jelem
!~ 			!
!~ 			! Function Variables
!~ 			!
!~ 			WRITE(7,*) 'damage: ', damage
!~ 			DO i1=1,D
!~ 				WRITE(7,*) 'damage gradient: ', damageGradient(i1)
!~ 			END DO
!~ 		  END IF
!~ 		  IF (penaltyEnergy .NE. penaltyEnergy) THEN
!~ 			WRITE(7,*) 'penalty Energy NAN: ', penaltyEnergy
!~ 			WRITE(7,*) 'INCREMENT: ', kinc
!~ 			WRITE(7,*) 'ELEMENT: ', jelem
!~ 			!
!~ 			! Function Variables
!~ 			!
!~ 			WRITE(7,*) 'damageJump: ', damageJump
!~ 		  END IF
!~ 		  IF (contactEnergy .NE. contactEnergy) THEN
!~ 			WRITE(7,*) 'contact Energy NAN: ', contactEnergy
!~ 			WRITE(7,*) 'INCREMENT: ', kinc
!~ 			WRITE(7,*) 'ELEMENT: ', jelem
!~ 			!
!~ 			! Function Variables
!~ 			!
!~ 			WRITE(7,*) 'separation: ', sep
!~ 		  END IF
!~ 		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 	  END IF
	  
	  
      

      DEALLOCATE(sep, sepGlobal, trac, K, tmp, tmp2, damageGradient, damageGradientLocal, damageGradientGlobal, parIEDMatrix, parCEDMatrix)

      CONTAINS

!------------------------------------------------------------------------------------

    END SUBROUTINE umatCZM

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION TransformedVectorLocalToGlobal(D,vec,TrafoMat)
      ! rotate/transform vector from local to global frame: v^glob = transposed(Q)*v^loc

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D
      REAL(kind=AbqRK), INTENT(IN) :: vec(D), TrafoMat(3,3)
      REAL(kind=AbqRK) :: TrafoMatD(D,D)
      DIMENSION TransformedVectorLocalToGlobal(D)

      TrafoMatD = Trafomat(1:D,1:D)

      TransformedVectorLocalToGlobal = zero
      ! rotate with local transformation matrix
      TransformedVectorLocalToGlobal = MATMUL(TRANSPOSE(TrafoMatD),vec)

    END FUNCTION TransformedVectorLocalToGlobal

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION TransformedTractionVectorLocalToStressArrayGlobal(D,trac,TrafoMat)
      ! rotate/transform (stress valued) traction vector from local to global frame via stress matrix
      ! not verified

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D
      REAL(kind=AbqRK), INTENT(IN) :: trac(D), TrafoMat(3,3)
      REAL(kind=AbqRK) :: stressMatrix(3,3), transformedStress(3,3)
      DIMENSION TransformedTractionVectorLocalToStressArrayGlobal(6)

      stressMatrix = zero
      IF (D .EQ. 2) THEN
        stressMatrix(1,2) = trac(1)
        stressMatrix(2,2) = trac(2)
        stressMatrix(2,1) = stressMatrix(1,2)
      ELSE IF (D .EQ. 3) THEN
        stressMatrix(1,3) = trac(1)
        stressMatrix(2,3) = trac(2)
        stressMatrix(3,3) = trac(3)
        stressMatrix(3,1) = stressMatrix(1,3)
        stressMatrix(3,2) = stressMatrix(2,3)
      END IF

      transformedStress = MATMUL(TRANSPOSE(TrafoMat),MATMUL(stressMatrix,TrafoMat))
      
      ! reorder in array according to Abaqus convention: 3D: 11-22-33-12-13-23, 2D: 11-22-33-12
      TransformedTractionVectorLocalToStressArrayGlobal = zero
      TransformedTractionVectorLocalToStressArrayGlobal(1) = transformedStress(1,1)
      TransformedTractionVectorLocalToStressArrayGlobal(2) = transformedStress(2,2)
      TransformedTractionVectorLocalToStressArrayGlobal(3) = transformedStress(3,3)
      TransformedTractionVectorLocalToStressArrayGlobal(4) = transformedStress(1,2)
      IF (D .EQ. 3) THEN
        TransformedTractionVectorLocalToStressArrayGlobal(5) = transformedStress(1,3)
        TransformedTractionVectorLocalToStressArrayGlobal(6) = transformedStress(2,3)
      END IF

    END FUNCTION TransformedTractionVectorLocalToStressArrayGlobal

!------------------------------------------------------------------------------------

REAL(kind=AbqRK) FUNCTION potential_CZ(D, normalDirectionIndex, tangentialDirectionIndex, damage, sep, damageGradient, damageJump, &
                     prop_df_one, prop_df_two, Hi_MM, prop_pen, damage_old, dtime, &
                     nCEDpar, parCEDMatrix, nIEDpar, parIED_MM_Matrix, nVDEDpar, parVDEDMatrix)
  !
  ! Computes total potential energy for cohesive zone with damage
  ! Sum of: CZED + IED + PenED + ContactED
  !
  USE CohesiveEnergyModule
  USE ContactEnergyModule
  USE AliasModuleCZ


  IMPLICIT NONE
  
  ! Input
  INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex, tangentialDirectionIndex
  REAL(kind=AbqRK), INTENT(IN) :: damage, damageJump
  REAL(kind=AbqRK), INTENT(IN) :: sep(D), damageGradient(D-1)
  REAL(kind=AbqRK), INTENT(IN) :: Hi_MM, prop_pen, damage_old, dtime
  INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, nIEDpar, nVDEDpar
  REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrix(nCEDpar), parIED_MM_Matrix(nIEDpar), parVDEDMatrix(nVDEDpar)
  REAL(kind=AbqRK), INTENT(IN) :: prop_df_one, prop_df_two

  ! Output
  
  ! Local
  REAL(kind=AbqRK) :: CZED_energy, IED_energy, PenED_energy, ContactED_energy, VDED_energy
  
  ! Cohesive Zone Energy Density
  CZED_energy = CZED(sep,D,normalDirectionIndex,tangentialDirectionIndex,nCEDpar,parCEDMatrix,damage,prop_df_one,prop_df_two)
  
  ! Interface Energy Density (gradient damage)
  IED_energy = IED(D-1, damage, damageGradient, nIEDpar, parIED_MM_Matrix)
  
  ! Penalty Energy Density (damage jump)
  PenED_energy = PenED(damageJump, prop_pen)
  
  ! Contact Energy Density (penalty for penetration)
  ContactED_energy = ContactED(D, normalDirectionIndex, sep, nCEDpar, parCEDMatrix)
  
  ! Viscous Energy Dis
  VDED_energy = VDED(damage,damage_old,dtime,nIEDpar,parIED_MM_Matrix,nVDEDpar,parVDEDMatrix)
  
  ! Total potential
  potential_CZ = CZED_energy + IED_energy + PenED_energy + ContactED_energy + VDED_energy
  
END FUNCTION potential_CZ

!------------------------------------------------------------------------------------
!
!							SUBROUTINES
!
!------------------------------------------------------------------------------------

SUBROUTINE Check_Stress_Comparison_CZ(stress, stran, D, ntens, normalDirectionIndex, &
										 tangentialDirectionIndex, prop_df_one, prop_df_two, &
										 Hi_MM, prop_pen, damage_old, dtime, nCEDpar, parCEDMatrix, &
										 nIEDpar, parIED_MM_Matrix, nVDEDpar, parVDEDMatrix)
                                       
  
  !
  ! Subroutine to check analytical stress against numerical stress
  ! Numerical stress is computed as derivative of total potential
  ! For Cohesive Zone with Damage model
  !
  !
  IMPLICIT NONE

  ! Input/Output variables
  REAL(kind=AbqRK), INTENT(IN) :: stress(:)
  REAL(kind=AbqRK), INTENT(IN) :: stran(:)
  INTEGER(kind=AbqIK), INTENT(IN) :: D
  INTEGER(kind=AbqIK), INTENT(IN) :: ntens
  INTEGER(kind=AbqIK), INTENT(IN) :: normalDirectionIndex, tangentialDirectionIndex
  REAL(kind=AbqRK), INTENT(IN) :: prop_df_one, prop_df_two
  REAL(kind=AbqRK), INTENT(IN) :: Hi_MM, prop_pen, damage_old, dtime
  INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, nIEDpar, nVDEDpar
  REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrix(nCEDpar)
  REAL(kind=AbqRK), INTENT(IN) :: parIED_MM_Matrix(nIEDpar)
  REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrix(nVDEDpar)
  
  ! Local variables
  REAL(kind=AbqRK), ALLOCATABLE :: stress_analytical(:), stress_numerical(:)
  REAL(kind=AbqRK), ALLOCATABLE :: stran_per(:)
  REAL(kind=AbqRK), ALLOCATABLE :: sep(:), sep_per(:)
  REAL(kind=AbqRK), ALLOCATABLE :: damageGradient(:), damageGradient_per(:)
  REAL(kind=AbqRK) :: damage, damage_per
  REAL(kind=AbqRK) :: damageJump, damageJump_per
  REAL(kind=AbqRK) :: delta, zero
  REAL(kind=AbqRK) :: pot_plus, pot_minus
  REAL(kind=AbqRK) :: sep_normal_unperturbed, sep_normal_plus, sep_normal_minus
  INTEGER(kind=AbqIK) :: i1, pos_damage, pos_damageJump, pos_damageGradient
  LOGICAL :: diff_stress, contact_jump
  
  ! Initialize
  zero = 0.d0
  diff_stress = .FALSE.
  contact_jump = .FALSE.
  delta = 1.d-7
  
  ! Allocate arrays with proper dimensions
  ALLOCATE(stress_analytical(ntens))
  ALLOCATE(stress_numerical(ntens))
  ALLOCATE(stran_per(ntens))
  ALLOCATE(sep(D))
  ALLOCATE(sep_per(D))
  ALLOCATE(damageGradient(D-1))
  ALLOCATE(damageGradient_per(D-1))
  
  ! Position in generalized arrays
  pos_damage         = ntens-(D+1)+1
  pos_damageJump     = ntens-(D+1)+2
  pos_damageGradient = ntens-(D+1)+3
  
  ! Extract current state
  sep(1:D) = stran(1:D)
  damage = stran(pos_damage)
  damageJump = stran(pos_damageJump)
  damageGradient = stran(pos_damageGradient:pos_damageGradient-2+D)
  
  ! Store unperturbed normal separation for contact jump check
  sep_normal_unperturbed = sep(normalDirectionIndex)
  
  !!!!!!!!!!!!
  !
  ! Check Comparison numerical stress / analytical stress
  !
  !!!!!!!!!!!!
  stress_analytical = stress
          
  ! Numerical stress via derivative of potential
  stress_numerical = zero
  DO i1=1,ntens
    ! Forward perturbation
    stran_per = stran
    stran_per(i1) = stran_per(i1) + delta
    
    ! Extract perturbed state
    sep_per(1:D) = stran_per(1:D)
    damage_per = stran_per(pos_damage)
    damageJump_per = stran_per(pos_damageJump)
    damageGradient_per = stran_per(pos_damageGradient:pos_damageGradient-2+D)
    
    ! Store perturbed normal separation
    sep_normal_plus = sep_per(normalDirectionIndex)
    
    ! Compute total potential (forward)
    pot_plus = potential_CZ(D, normalDirectionIndex, tangentialDirectionIndex, damage_per, sep_per, damageGradient_per, &
                           damageJump_per, prop_df_one, prop_df_two, &
                           Hi_MM, prop_pen, damage_old, dtime, &
                           nCEDpar, parCEDMatrix, nIEDpar, parIED_MM_Matrix, nVDEDpar, parVDEDMatrix)
    
    ! Backward perturbation
    stran_per = stran
    stran_per(i1) = stran_per(i1) - delta
    
    ! Extract perturbed state
    sep_per(1:D) = stran_per(1:D)
    damage_per = stran_per(pos_damage)
    damageJump_per = stran_per(pos_damageJump)
    damageGradient_per = stran_per(pos_damageGradient:pos_damageGradient-2+D)
    
    ! Store perturbed normal separation
    sep_normal_minus = sep_per(normalDirectionIndex)
    
    ! Compute total potential (backward)
    pot_minus = potential_CZ(D, normalDirectionIndex, tangentialDirectionIndex, damage_per, sep_per, damageGradient_per, &
                           damageJump_per, prop_df_one, prop_df_two, &
                           Hi_MM, prop_pen, damage_old, dtime, &
                           nCEDpar, parCEDMatrix, nIEDpar, parIED_MM_Matrix, nVDEDpar, parVDEDMatrix)
    
    ! Check for contact jump (sign change in normal separation due to perturbation)
    IF (i1 == normalDirectionIndex) THEN
      IF (sep_normal_plus * sep_normal_minus .LT. zero) THEN
        contact_jump = .TRUE.
      END IF
    END IF
    
    ! Central difference
    stress_numerical(i1) = (pot_plus - pot_minus)/(2.d0*delta)
  END DO
  
!~   IF (Incrementnumber .GT. 1) THEN
!~ 	  WRITE(6,*) '==================================='
!~ 	  WRITE(6,*) 'Print Test Stress Subroutine'
!~ 	  WRITE(6,*) 'Stress analytical : '
!~ 	  WRITE(6,*) stress_analytical(1)
!~ 	  WRITE(6,*) stress_analytical(2)
!~ 	  WRITE(6,*) stress_analytical(3)
!~ 	  WRITE(6,*) stress_analytical(4)
!~ 	  WRITE(6,*) stress_analytical(5)
!~ 	  WRITE(6,*) 'Stress numerical : '
!~ 	  WRITE(6,*) stress_numerical(1)
!~ 	  WRITE(6,*) stress_numerical(2)
!~ 	  WRITE(6,*) stress_numerical(3)
!~ 	  WRITE(6,*) stress_numerical(4)
!~ 	  WRITE(6,*) stress_numerical(5)
!~   END IF
  !
  !
  ! Print Detailed Information if deviations found
  !
  IF (Incrementnumber .GT. 1 .AND. .NOT. contact_jump) THEN
      DO i1 = 1, ntens
          IF (i1 .NE. 4 .AND. ABS(stress_analytical(i1) - stress_numerical(i1)) .GT. 0.001) THEN
          !
          ! NE 4 => 4 is driving force of the penalty energy
          ! small increment delta already changes derivative a lot
          ! No respective delta for pen-energy coded, so no check to avoid stops
          !
             diff_stress = .TRUE.
             
             WRITE(6,*) "Increment: ", Incrementnumber
             WRITE(6,*) "Element: ", Elementnumber
             WRITE(6,*) "IP: ", Integrationpointnumber
             WRITE(6,'("Abweichung bei Komponente ",I0,": ",2(E16.8,1X))') i1, stress_analytical(i1), stress_numerical(i1)
             
             ! Localization of deviation
             IF (i1 >= 1 .AND. i1 <= D) THEN
                WRITE(6,*) '  --> Abweichung in Traktionsvektor (Separation)'
                IF (i1 == normalDirectionIndex) THEN
                   WRITE(6,*) '      (Normalkomponente)'
                ELSE
                   WRITE(6,*) '      (Tangentialkomponente)'
                END IF
             ELSE IF (i1 == pos_damage) THEN
                WRITE(6,*) '  --> Abweichung in Damage Driving Force'
             ELSE IF (i1 == pos_damageJump) THEN
                WRITE(6,*) '  --> Abweichung in Damage Jump Konjugat'
             ELSE IF (i1 >= pos_damageGradient .AND. i1 <= pos_damageGradient+D-2) THEN
                WRITE(6,*) '  --> Abweichung in Damage Gradient Konjugat'
             END IF
          END IF
      END DO 
      
      IF (diff_stress) THEN
         WRITE(6,*), "---------------------------------------"
         WRITE(6,*) 'Separationsvektor: '
         DO i1 = 1, D
            WRITE(6,'(E16.8)') sep(i1)
         END DO
         WRITE(6,*) 'Damage: ', damage
         WRITE(6,*) 'Damage Jump: ', damageJump
         WRITE(6,*) 'Damage Gradient: '
         DO i1 = 1, D-1
            WRITE(6,'(E16.8)') damageGradient(i1)
         END DO
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
  ELSE IF (contact_jump) THEN
      WRITE(6,*) "Contact jump detected during perturbation - skipping stress check"
      WRITE(6,*) "sep_normal: ", sep_normal_unperturbed
  END IF
  
  ! Deallocate arrays
  DEALLOCATE(stress_analytical, stress_numerical, stran_per, sep, sep_per, &
             damageGradient, damageGradient_per)
  
END SUBROUTINE Check_Stress_Comparison_CZ

!------------------------------------------------------------------------------------

SUBROUTINE check_matrix_entries(Ct, D, pos_damage, pos_damageGradient, pos_damageJump, &
                                sep, normalDirectionIndex, tangentialDirectionIndex, &
                                nCEDpar, parCEDMatrix, damage, prop_df_one, prop_df_two, &
                                Hi_MM, damageGradient, damage_old, dtime, nIEDpar, &
                                parIED_MM_Matrix, nVDEDpar, parVDEDMatrix, damageJump, prop_pen)
                                
    ! Überprüft alle Matrixeinträge aus Ct
    ! Für NAN wird zugeordnet 
    
    ! In welcher Submatrix befindet sich NAN?
    ! Welcher Teilbeitrag aus welcher Energie ist dafür verantwortlich?
    ! Welche Werte haben die Parameter dieser Energie?
    USE ABQINTERFACE
    USE FLOATNUMBERS
    USE ABAModul
    
    USE CohesiveEnergyModule                        
    USE AliasModuleCZ    
                                
    IMPLICIT NONE
    
    ! Argumente
    INTEGER(kind=AbqIK), INTENT(IN) :: D
    INTEGER(kind=AbqIK), INTENT(IN) :: pos_damage
    INTEGER(kind=AbqIK), INTENT(IN) :: pos_damageGradient
    INTEGER(kind=AbqIK), INTENT(IN) :: pos_damageJump
    INTEGER(kind=AbqIK), INTENT(IN) :: normalDirectionIndex, tangentialDirectionIndex
    INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, nIEDpar, nVDEDpar
    REAL(kind=AbqRK), INTENT(IN) :: Ct(:,:)
    REAL(kind=AbqRK), INTENT(IN) :: sep(:), damage, damage_old, dtime, damageJump
    REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrix(:), prop_df_one, prop_df_two
    REAL(kind=AbqRK), INTENT(IN) :: Hi_MM, damageGradient(:)
    REAL(kind=AbqRK), INTENT(IN) :: parIED_MM_Matrix(:), parVDEDMatrix(:), prop_pen
    
    ! Lokale Variablen
    INTEGER(kind=AbqIK) :: i1, i2, i, j
    INTEGER(kind=AbqIK) :: matrix_size
    LOGICAL :: found_nan
    
    ! Vorberechnete einzelne Funktionsbeiträge
    REAL(kind=AbqRK), ALLOCATABLE :: contrib_d_CZED_d_sep_d_sep(:,:)
    REAL(kind=AbqRK)  :: contrib_d_CZED_d_damage_d_damage_Hi
    REAL(kind=AbqRK)  :: contrib_d_IED_d_phase_d_phase
    REAL(kind=AbqRK)  :: contrib_d_VDED_d_phase_d_phase
    REAL(kind=AbqRK), ALLOCATABLE :: contrib_d_IED_d_grad_phase_d_grad_phase(:,:)
    REAL(kind=AbqRK) :: contrib_d_PenED_d_damageJump_d_damageJump
    REAL(kind=AbqRK), ALLOCATABLE :: contrib_d_CZED_d_sep_d_damage(:)
    
    matrix_size = SIZE(Ct, 1)
    
    ! Allokiere Arrays für die Beiträge
    ALLOCATE(contrib_d_CZED_d_sep_d_sep(D, D))
    ALLOCATE(contrib_d_IED_d_grad_phase_d_grad_phase(D-1, D-1))
    ALLOCATE(contrib_d_CZED_d_sep_d_damage(D))


    ! ========================================================================
    ! Berechne alle Funktionsbeiträge einmalig zu Beginn
    ! ========================================================================
    
    contrib_d_CZED_d_sep_d_sep = d_CZED_d_sep_d_sep(sep, D, normalDirectionIndex, &
                                                      tangentialDirectionIndex, nCEDpar, &
                                                      parCEDMatrix, damage, prop_df_one, prop_df_two)
    
    contrib_d_CZED_d_damage_d_damage_Hi = d_CZED_d_damage_d_damage_Hi(sep, D, normalDirectionIndex, &
                                                                        tangentialDirectionIndex, nCEDpar, &
                                                                        parCEDMatrix, damage, prop_df_one, &
                                                                        prop_df_two, Hi_MM)
    
    contrib_d_IED_d_phase_d_phase = d_IED_d_phase_d_phase(D-1, damage, damageGradient, &
                                                           nIEDpar, parIED_MM_Matrix)
    
    contrib_d_VDED_d_phase_d_phase = d_VDED_d_phase_d_phase(damage, damage_old, dtime, &
                                                             nIEDpar, parIED_MM_Matrix, &
                                                             nVDEDpar, parVDEDMatrix)
    
    contrib_d_IED_d_grad_phase_d_grad_phase = d_IED_d_grad_phase_d_grad_phase(D-1, damage, &
                                                                               damageGradient, nIEDpar, &
                                                                               parIED_MM_Matrix)
    
    contrib_d_PenED_d_damageJump_d_damageJump = d_PenED_d_damageJump_d_damageJump(damageJump, prop_pen)
    
    contrib_d_CZED_d_sep_d_damage = d_CZED_d_sep_d_damage(sep, D, normalDirectionIndex, &
                                                           tangentialDirectionIndex, nCEDpar, &
                                                           parCEDMatrix, damage, prop_df_one, prop_df_two)
    
    ! ========================================================================
    ! Durchlaufe alle Einträge der Matrix
    ! ========================================================================
    
    DO i1 = 1, matrix_size
        DO i2 = 1, matrix_size
            
            ! (NAN Check)
            IF (Ct(i1,i2) .NE. Ct(i1,i2)) THEN
                
                DO i = 1, SIZE(Ct,1)
				   WRITE(6,'(100(ES12.4,1X))') (Ct(i,j), j = 1, SIZE(Ct,2))
				END DO
                
                
                
                ! Fallunterscheidung basierend auf Indizes
                SELECT CASE (get_case_number(i1, i2, D, pos_damage, pos_damageGradient, pos_damageJump))
                    
                    CASE (1)
                        ! Case 1: Ct(1:D, 1:D)
                        WRITE(6,*) '=== Case 1, d_sep_sep: i1=', i1, ', i2=', i2, ' ==='
                        
                        ! Prüfe den entsprechenden Beitrag
                        IF (contrib_d_CZED_d_sep_d_sep(i1, i2) .NE. contrib_d_CZED_d_sep_d_sep(i1, i2)) THEN
                            WRITE(6,*) 'contrib_d_CZED_d_sep_d_sep(', i1, ',', i2, ') ist NAN!'
                            WRITE(6,*) 'Prüfe Parameter von d_CZED_d_sep_d_sep:'
                            
                            ! Parameterprüfung für d_CZED_d_sep_d_sep
                            DO i = 1, SIZE(sep)
                                IF (sep(i) .NE. sep(i)) THEN
                                    WRITE(6,*) '  FEHLER: sep(', i, ') ist NAN!'
                                END IF
                            END DO
                            IF (damage .NE. damage) THEN
                                WRITE(6,*) '  FEHLER: damage ist NAN!'
                            END IF
                        END IF
                        !
                        !
                        !
                        CALL XEXIT()
                    !
                    !:::::::::::::::::::::::
                    !
                    CASE (2)
                        ! Case 2: Ct(pos_damage, pos_damage)
                        WRITE(6,*) '=== Case 2, d_phase_phase: i1=', i1, ', i2=', i2, ' ==='
                        
                        ! Prüfe ersten Beitrag: d_CZED_d_damage_d_damage_Hi
                        IF (contrib_d_CZED_d_damage_d_damage_Hi .NE. contrib_d_CZED_d_damage_d_damage_Hi) THEN
                            WRITE(6,*) 'contrib_d_CZED_d_damage_d_damage_Hi ist NAN!'
                            WRITE(6,*) 'Prüfe Parameter von d_CZED_d_damage_d_damage_Hi:'
                            
                            ! Parameterprüfung (Dummy-Beispiel)
                            IF (Hi_MM .NE. Hi_MM) THEN
                                WRITE(6,*) '  FEHLER: Hi_MM ist NAN!'
                            END IF
                        END IF
                        
                        ! Prüfe zweiten Beitrag: d_IED_d_phase_d_phase
                        IF (contrib_d_IED_d_phase_d_phase .NE. contrib_d_IED_d_phase_d_phase) THEN
                            WRITE(6,*) 'contrib_d_IED_d_phase_d_phase ist NAN!'
                            WRITE(6,*) 'Prüfe Parameter von d_IED_d_phase_d_phase:'
                            
                            ! Parameterprüfung (Dummy-Beispiel)
                            IF (damage .NE. damage) THEN
                                WRITE(6,*) '  FEHLER: damage ist NAN!'
                            END IF
                            DO i = 1, SIZE(damageGradient)
                                IF (damageGradient(i) .NE. damageGradient(i)) THEN
                                    WRITE(6,*) '  FEHLER: damageGradient(', i, ') ist NAN!'
                                END IF
                            END DO
                            IF (parIED_MM_Matrix(1) .NE. parIED_MM_Matrix(1)) THEN
								WRITE(6,*) 'FEHLER: Gieff ist NAN'
							END IF
                        END IF
                        
                        ! Prüfe dritten Beitrag: d_VDED_d_phase_d_phase
                        IF (contrib_d_VDED_d_phase_d_phase .NE. contrib_d_VDED_d_phase_d_phase) THEN
                            WRITE(6,*) 'contrib_d_VDED_d_phase_d_phase ist NAN!'
                            WRITE(6,*) 'Prüfe Parameter von d_VDED_d_phase_d_phase:'
                            
                            ! Parameterprüfung (Dummy-Beispiel)
                            IF (damage .NE. damage) THEN
                                WRITE(6,*) '  FEHLER: damage ist NAN!'
                            END IF
                            IF (damage_old .NE. damage_old) THEN
                                WRITE(6,*) '  FEHLER: damage_old ist NAN!'
                            END IF
                            IF (dtime .LE. zero) THEN
                                WRITE(6,*) '  FEHLER: dtime ist LE zero!'
                            END IF
                            IF (parIED_MM_Matrix(1) .NE. parIED_MM_Matrix(1)) THEN
								WRITE(6,*) 'FEHLER: Gieff ist NAN'
							END IF
                        END IF
                        !
                        !
                        !
                        CALL XEXIT()
                    !
                    !:::::::::::::::::::::::
                    !
                    CASE (3)
                        ! Case 3: Ct(pos_damageGradient:pos_damageGradient+D-2, pos_damageGradient:pos_damageGradient+D-2)
                        WRITE(6,*) '=== Case 3 d_gradphase_gradphase: i1=', i1, ', i2=', i2, ' ==='
                        
                        ! Prüfe den entsprechenden Beitrag
                        IF (contrib_d_IED_d_grad_phase_d_grad_phase(i1-pos_damageGradient+1, i2-pos_damageGradient+1) .NE. &
                            contrib_d_IED_d_grad_phase_d_grad_phase(i1-pos_damageGradient+1, i2-pos_damageGradient+1)) THEN
                            WRITE(6,*) 'contrib_d_IED_d_grad_phase_d_grad_phase(', &
                                       i1-pos_damageGradient+1, ',', i2-pos_damageGradient+1, ') ist NAN!'
                            WRITE(6,*) 'Prüfe Parameter von d_IED_d_grad_phase_d_grad_phase:'
                            
                            ! Parameterprüfung (Dummy-Beispiel)
                            IF (damage .NE. damage) THEN
                                WRITE(6,*) '  FEHLER: damage ist NAN!'
                            END IF
                            DO i = 1, SIZE(damageGradient)
                                IF (damageGradient(i) .NE. damageGradient(i)) THEN
                                    WRITE(6,*) '  FEHLER: damageGradient(', i, ') ist NAN!'
                                END IF
                            END DO
                            IF (parIED_MM_Matrix(1) .NE. parIED_MM_Matrix(1)) THEN
								WRITE(6,*) 'FEHLER: Gieff ist NAN'
							END IF
                        END IF                        
                        !
                        !
                        !
                        CALL XEXIT()
                    !
                    !:::::::::::::::::::::::
                    !
                    CASE (4)
                        ! Case 4: Ct(pos_damageJump, pos_damageJump)
                        WRITE(6,*) '=== Case 4 d_damageJump_damageJump: i1=', i1, ', i2=', i2, ' ==='
                        
                        ! Prüfe den Beitrag
                        IF (contrib_d_PenED_d_damageJump_d_damageJump .NE. contrib_d_PenED_d_damageJump_d_damageJump) THEN
                            WRITE(6,*) 'contrib_d_PenED_d_damageJump_d_damageJump ist NAN!'
                            WRITE(6,*) 'Prüfe Parameter von d_PenED_d_damageJump_d_damageJump:'
                            
                            ! Parameterprüfung (Dummy-Beispiel)
                            IF (damageJump .NE. damageJump) THEN
                                WRITE(6,*) '  FEHLER: damageJump ist NAN!'
                            END IF
                            IF (prop_pen .NE. prop_pen) THEN
                                WRITE(6,*) '  FEHLER: prop_pen ist NAN!'
                            END IF
                        END IF
                        !
                        !
                        !
                        CALL XEXIT()
                    !
                    !:::::::::::::::::::::::
                    !
                    CASE (5)
                        ! Case 5: Ct(1:D, pos_damage)
                        WRITE(6,*) '=== Case 5 d_sep_damage: i1=', i1, ', i2=', i2, ' ==='
                        WRITE(6,*) 'VALUE: ', Ct(i1,i2)
                        ! Prüfe den entsprechenden Beitrag
                        IF (contrib_d_CZED_d_sep_d_damage(i1) .NE. contrib_d_CZED_d_sep_d_damage(i1)) THEN
                            WRITE(6,*) 'contrib_d_CZED_d_sep_d_damage(', i1, ') ist NAN!'
                            WRITE(6,*) 'Prüfe Parameter von d_CZED_d_sep_d_damage:'
                            
                            ! Parameterprüfung (Dummy-Beispiel)
                            DO i = 1, SIZE(sep)
                                IF (sep(i) .NE. sep(i)) THEN
                                    WRITE(6,*) '  FEHLER: sep(', i, ') ist NAN!'
                                END IF
                            END DO
                            IF (damage .NE. damage) THEN
                                WRITE(6,*) '  FEHLER: damage ist NAN!'
                            END IF
                        END IF
                        
                    CASE DEFAULT
						WRITE(6,*) 'NAN ROUTINE'
                        WRITE(6,*) 'Warning: Kein Case definiert für i1=', i1, ', i2=', i2
                        !
                        !
                        !
                        CALL XEXIT()
                END SELECT
                
            END IF
            
        END DO
    END DO
    
    ! Speicher freigeben
    DEALLOCATE(contrib_d_CZED_d_sep_d_sep)
    DEALLOCATE(contrib_d_IED_d_grad_phase_d_grad_phase)
    DEALLOCATE(contrib_d_CZED_d_sep_d_damage)
    
END SUBROUTINE check_matrix_entries

!--------------------------------------------------------------------------------------------

! Hilfsfunktion zur Bestimmung der Case-Nummer
INTEGER FUNCTION get_case_number(i1, i2, D, pos_damage, pos_damageGradient, pos_damageJump)
    IMPLICIT NONE
    
    INTEGER(kind=AbqIK), INTENT(IN) :: i1, i2, D
    INTEGER(kind=AbqIK), INTENT(IN) :: pos_damage, pos_damageGradient, pos_damageJump
    
    IF (i1 >= 1 .AND. i1 <= D .AND. i2 >= 1 .AND. i2 <= D) THEN
        get_case_number = 1
        RETURN
    END IF
    
    IF (i1 == pos_damage .AND. i2 == pos_damage) THEN
        get_case_number = 2
        RETURN
    END IF
    
    IF (i1 >= pos_damageGradient .AND. i1 <= pos_damageGradient+D-2 .AND. &
        i2 >= pos_damageGradient .AND. i2 <= pos_damageGradient+D-2) THEN
        get_case_number = 3
        RETURN
    END IF
    
    IF (i1 == pos_damageJump .AND. i2 == pos_damageJump) THEN
        get_case_number = 4
        RETURN
    END IF
    
    IF (i1 >= 1 .AND. i1 <= D .AND. i2 == pos_damage) THEN
        get_case_number = 5
        RETURN
    END IF
    
    get_case_number = 0
    
END FUNCTION get_case_number

!------------------------------------------------------------------------------------

END MODULE CZUMAT_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
