!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! bulk energy module
!
! Martin Olbricht, TU Bergakademie Freiberg, 12.01.2024
!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BulkEnergyModule

  USE SplitEnergyModule, ONLY: HFEDpos => HFEDposNoSplit, &
                                               d_HFEDpos_d_eps_e => d_HFEDposNoSplit_d_eps_e, &
                                               d_d_HFEDpos_d_eps_e_d_eps_e => d_d_HFEDposNoSplit_d_eps_e_d_eps_e, &
                                               HFEDneg => HFEDnegNoSplit, &
                                               d_HFEDneg_d_eps_e => d_HFEDnegNoSplit_d_eps_e, &
                                               d_d_HFEDneg_d_eps_e_d_eps_e => d_d_HFEDnegNoSplit_d_eps_e_d_eps_e
                                               

!~   USE SplitEnergyModule, ONLY: HFEDpos => HFEDposAmorSplit, &
!~                                                d_HFEDpos_d_eps_e => d_HFEDposAmorSplit_d_eps_e, &
!~                                                d_d_HFEDpos_d_eps_e_d_eps_e => d_d_HFEDposAmorSplit_d_eps_e_d_eps_e, &
!~                                                HFEDneg => HFEDnegAmorSplit, &
!~                                                d_HFEDneg_d_eps_e => d_HFEDnegAmorSplit_d_eps_e, &
!~                                                d_d_HFEDneg_d_eps_e_d_eps_e => d_d_HFEDnegAmorSplit_d_eps_e_d_eps_e
                                               
                                               
!~   USE SplitEnergyModule, ONLY: HFEDpos => HFEDposMieheSplit, &
!~  											  d_HFEDpos_d_eps_e => d_HFEDposMieheSplit_NOTRANSFORMATION_d_eps_e, &
!~                                                d_d_HFEDpos_d_eps_e_d_eps_e => d_d_HFEDposAmorSplit_d_eps_e_d_eps_e, &
!~                                                HFEDneg => HFEDnegMieheSplit, &
!~ 											      d_HFEDneg_d_eps_e => d_HFEDnegMieheSplit_NOTRANSFORMATION_d_eps_e, &
!~                                                d_d_HFEDneg_d_eps_e_d_eps_e => d_d_HFEDnegAmorSplit_d_eps_e_d_eps_e
                                               
                                               
!~   USE SplitEnergyModule, ONLY: HFEDpos => HFEDposStarSplit, &
!~                                                d_HFEDpos_d_eps_e => d_HFEDposStarSplit_d_eps_e, &
!~                                                d_d_HFEDpos_d_eps_e_d_eps_e => d_d_HFEDposStarSplit_d_eps_e_d_eps_e, &
!~                                                HFEDneg => HFEDnegStarSplit, &
!~                                                d_HFEDneg_d_eps_e => d_HFEDnegStarSplit_d_eps_e, &
!~                                                d_d_HFEDneg_d_eps_e_d_eps_e => d_d_HFEDnegStarSplit_d_eps_e_d_eps_e

                                               
  USE DegradationFunctionModule, ONLY: degF => quad_degF, &
                                           d_degF_d_phase => d_quad_degF_d_phase, &
                                           d_d_degF_d_phase_d_phase => d_d_quad_degF_d_phase_d_phase                 

!~   USE DegradationFunctionModule, ONLY: degF => No_degF, &
!~                                            d_degF_d_phase => d_No_degF_d_phase, &
!~                                            d_d_degF_d_phase_d_phase => d_d_No_degF_d_phase_d_phase    


  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS

!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION HFEDtens_H(eps,nHFEDpar,parHFEDMatrixPhase)
    ! bulk energy density

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE FreeEnergyModule
      USE SplitEnergyModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      !
      HFEDtens_H = HFEDpos(eps,nHFEDpar,parHFEDMatrixPhase)
    END FUNCTION HFEDtens_H
    
!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION bulkED(eps,phase,nHFEDpar,parHFEDMatrixPhase)
    ! bulk energy density

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE FreeEnergyModule
      USE SplitEnergyModule
      USE DegradationFunctionModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: HFEDtens, HFEDcomp
      REAL(kind=AbqRK) :: degD
      !
	  !
	  !
      HFEDtens = HFEDpos(eps,nHFEDpar,parHFEDMatrixPhase)
      HFEDcomp = HFEDneg(eps,nHFEDpar,parHFEDMatrixPhase)
      degD = degF(phase)
      !
      bulkED = degD*HFEDtens + HFEDcomp
      !
!~       WRITE(7,*) '------------------------------------------------- '
!~       WRITE(7,*) '-------------------bulkED------------------------ '
!~       WRITE(7,*) '------------------------------------------------- '
!~ 	  WRITE(7,*) 'eps: ', eps
!~ 	  WRITE(7,*) 'phase: ', phase
!~ 	  WRITE(7,*) 'HFEDtens: ', HFEDtens
!~ 	  WRITE(7,*) 'HFEDcomp: ', HFEDcomp
!~ 	  WRITE(7,*) 'bulkED: ', degD*HFEDtens + HFEDcomp
!~ 	  WRITE(7,*) '------------------------------------------------- '
      !
      
    END FUNCTION bulkED
    
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!
!                               first derivatives
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION d_bulkED_d_eps(eps,phase,nHFEDpar,parHFEDMatrixPhase)
    ! derivative of bulk energy density w.r.t. strain

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE FreeEnergyModule
      USE SplitEnergyModule
      USE DegradationFunctionModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: degD
      REAL(kind=AbqRK) :: d_HFEDtens_d_eps_e(3,3), d_HFEDcomp_d_eps_e(3,3)
      DIMENSION d_bulkED_d_eps(3,3)
      !
      d_HFEDtens_d_eps_e = d_HFEDpos_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      d_HFEDcomp_d_eps_e = d_HFEDneg_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      !
      degD = degF(phase)
      !
      !
      d_bulkED_d_eps = degD*d_HFEDtens_d_eps_e + d_HFEDcomp_d_eps_e
      !
!~ 	  WRITE(7,*) '-------------------------------------------------'
!~ 	  WRITE(7,*) '---------------d_bulkED_d_eps--------------------'
!~ 	  WRITE(7,*) '-------------------------------------------------'
!~ 	  WRITE(7,*) 'eps: ', eps
!~ 	  WRITE(7,*) 'phase: ', phase
!~ 	  WRITE(7,*) 'd_HFEDtens_d_eps_e: ', d_HFEDtens_d_eps_e
!~ 	  WRITE(7,*) 'd_HFEDcomp_d_eps_e: ', d_HFEDcomp_d_eps_e
!~ 	  WRITE(7,*) 'd_bulkED_d_eps: ', degD*d_HFEDtens_d_eps_e + d_HFEDcomp_d_eps_e
!~ 	  WRITE(7,*) '------------------------------------------------- '
      !

    END FUNCTION d_bulkED_d_eps

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_bulkED_d_phase(eps,phase,nHFEDpar,parHFEDMatrixPhase)
    ! derivative of bulk energy density w.r.t. phase

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE FreeEnergyModule
      USE DegradationFunctionModule
      USE SplitEnergyModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: d_degD_d_phase, HFEDtens
      !
      HFEDtens = HFEDpos(eps,nHFEDpar,parHFEDMatrixPhase)
      d_degD_d_phase = d_degF_d_phase(phase)
      !
      d_bulkED_d_phase = d_degD_d_phase*HFEDtens
!      d_bulkED_d_phase = zero
    END FUNCTION d_bulkED_d_phase
    
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_bulkED_d_phase_H(eps,phase,nHFEDpar,parHFEDMatrixPhase,H)
    ! derivative of bulk energy density w.r.t. phase

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE FreeEnergyModule
      USE DegradationFunctionModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: phase,H
      REAL(kind=AbqRK) :: d_degD_d_phase_H, HFEDtens
      !
      HFEDtens = H
      d_degD_d_phase_H = d_degF_d_phase(phase)
      !
      d_bulkED_d_phase_H = d_degD_d_phase_H*HFEDtens
!      d_bulkED_d_phase = zero
    END FUNCTION d_bulkED_d_phase_H
    

!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!
!                               second derivatives
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!


    REAL(kind=AbqRK) FUNCTION d_d_bulkED_d_eps_d_eps(eps,phase,nHFEDpar,parHFEDMatrixPhase)
    ! doubled partial derivative of bulk energy density w.r.t. strain strain

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE FreeEnergyModule
      USE SplitEnergyModule
      USE DegradationFunctionModule
      USE SharedValues ! Debug
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqIK) :: i1, i2 ! Debug Print
      REAL(kind=AbqRK) :: d_d_HFEDtens_d_eps_e_d_eps_e(3,3,3,3), d_d_HFEDcomp_d_eps_e_d_eps_e(3,3,3,3), stiffness_degraded_dummy(3,3,3,3), stiffness_degraded_dummy_voigt(6,6), stiffness_dummy(3,3,3,3), stiffness_dummy_voigt(6,6)! debug dummys
      REAL(kind=AbqRK) :: degD
      DIMENSION d_d_bulkED_d_eps_d_eps(3,3,3,3)
      !
      !
      d_d_HFEDtens_d_eps_e_d_eps_e = d_d_HFEDpos_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      d_d_HFEDcomp_d_eps_e_d_eps_e = d_d_HFEDneg_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      !
!~       IF (Incrementnumber .GT. 2 .AND. Elementnumber .EQ. 1 .AND. Integrationpointnumber .EQ. 1) THEN
!~ 		stiffness_dummy = d_d_HFEDtens_d_eps_e_d_eps_e
!~ 		stiffness_dummy_voigt(1,1) = stiffness_dummy(1,1,1,1)
!~ 		stiffness_dummy_voigt(1,2) = stiffness_dummy(1,1,2,2)
!~ 		stiffness_dummy_voigt(1,3) = stiffness_dummy(1,1,3,3)
!~ 		stiffness_dummy_voigt(1,4) = stiffness_dummy(1,1,1,2)
!~ 		stiffness_dummy_voigt(1,5) = stiffness_dummy(1,1,1,3)
!~ 		stiffness_dummy_voigt(1,6) = stiffness_dummy(1,1,2,3)
!~ 		stiffness_dummy_voigt(2,1) = stiffness_dummy(2,2,1,1)
!~ 		stiffness_dummy_voigt(2,2) = stiffness_dummy(2,2,2,2)
!~ 		stiffness_dummy_voigt(2,3) = stiffness_dummy(2,2,3,3)
!~ 		stiffness_dummy_voigt(2,4) = stiffness_dummy(2,2,1,2)
!~ 		stiffness_dummy_voigt(2,5) = stiffness_dummy(2,2,1,3)
!~ 		stiffness_dummy_voigt(2,6) = stiffness_dummy(2,2,2,3)
!~ 		stiffness_dummy_voigt(3,1) = stiffness_dummy(3,3,1,1)
!~ 		stiffness_dummy_voigt(3,2) = stiffness_dummy(3,3,2,2)
!~ 		stiffness_dummy_voigt(3,3) = stiffness_dummy(3,3,3,3)
!~ 		stiffness_dummy_voigt(3,4) = stiffness_dummy(3,3,1,2)
!~ 		stiffness_dummy_voigt(3,5) = stiffness_dummy(3,3,1,3)
!~ 		stiffness_dummy_voigt(3,6) = stiffness_dummy(3,3,2,3)
!~ 		stiffness_dummy_voigt(4,1) = stiffness_dummy(2,1,1,1)
!~ 		stiffness_dummy_voigt(4,2) = stiffness_dummy(2,1,2,2)
!~ 		stiffness_dummy_voigt(4,3) = stiffness_dummy(2,1,3,3)
!~ 		stiffness_dummy_voigt(4,4) = stiffness_dummy(2,1,1,2)
!~ 		stiffness_dummy_voigt(4,5) = stiffness_dummy(2,1,1,3)
!~ 		stiffness_dummy_voigt(4,6) = stiffness_dummy(2,1,2,3)
!~ 		stiffness_dummy_voigt(5,1) = stiffness_dummy(3,1,1,1)
!~ 		stiffness_dummy_voigt(5,2) = stiffness_dummy(3,1,2,2)
!~ 		stiffness_dummy_voigt(5,3) = stiffness_dummy(3,1,3,3)
!~ 		stiffness_dummy_voigt(5,4) = stiffness_dummy(3,1,1,2)
!~ 		stiffness_dummy_voigt(5,5) = stiffness_dummy(3,1,1,3)
!~ 		stiffness_dummy_voigt(5,6) = stiffness_dummy(3,1,2,3)
!~ 		stiffness_dummy_voigt(6,1) = stiffness_dummy(3,2,1,1)
!~ 		stiffness_dummy_voigt(6,2) = stiffness_dummy(3,2,2,2)
!~ 		stiffness_dummy_voigt(6,3) = stiffness_dummy(3,2,3,3)
!~ 		stiffness_dummy_voigt(6,4) = stiffness_dummy(3,2,1,2)
!~ 		stiffness_dummy_voigt(6,5) = stiffness_dummy(3,2,1,3)
!~ 		stiffness_dummy_voigt(6,6) = stiffness_dummy(3,2,2,3)
!~ 		WRITE(6,*), "---------------------------------------"
!~ 		WRITE(6,*) 'stiffness: '
!~ 		DO i1 = 1, 6
!~ 			WRITE(6,'(7(F16.4,1X))') (stiffness_dummy_voigt(i1,i2), i2=1,6)
!~ 		END DO
!~ 		WRITE(6,*) 'phase: ', phase
!~ 		WRITE(6,*) 'degF: ', degD
!~ 		WRITE(6,*), "---------------------------------------"
!~ 	  END IF

      degD = degF(phase)
      !
      !
!~       WRITE(7,*) '------------------------------------------------- '
!~       WRITE(7,*) '----------d_d_bulkED_d_eps_d_eps----------------- '
!~       WRITE(7,*) '------------------------------------------------- '
!~ 	  WRITE(7,*) 'eps: ', eps
!~ 	  WRITE(7,*) 'phase: ', phase
!~ 	  WRITE(7,*) 'd_d_HFEDtens_d_eps_e_d_eps_e: ', d_d_HFEDtens_d_eps_e_d_eps_e
!~ 	  WRITE(7,*) 'd_d_HFEDcomp_d_eps_e_d_eps_e: ', d_d_HFEDcomp_d_eps_e_d_eps_e
!~ 	  WRITE(7,*) 'degD: ', degD
!~ 	  WRITE(7,*) '------------------------------------------------- '
      !
      !
      d_d_bulkED_d_eps_d_eps = degD*d_d_HFEDtens_d_eps_e_d_eps_e + d_d_HFEDcomp_d_eps_e_d_eps_e
      
      
            !		!
!~ 	  IF (Incrementnumber .GT. 2 .AND. Elementnumber .EQ. 1 .AND. Integrationpointnumber .EQ. 1) THEN
!~ 		stiffness_degraded_dummy = d_d_bulkED_d_eps_d_eps
!~ 		stiffness_degraded_dummy_voigt(1,1) = stiffness_degraded_dummy(1,1,1,1)
!~ 		stiffness_degraded_dummy_voigt(1,2) = stiffness_degraded_dummy(1,1,2,2)
!~ 		stiffness_degraded_dummy_voigt(1,3) = stiffness_degraded_dummy(1,1,3,3)
!~ 		stiffness_degraded_dummy_voigt(1,4) = stiffness_degraded_dummy(1,1,1,2)
!~ 		stiffness_degraded_dummy_voigt(1,5) = stiffness_degraded_dummy(1,1,1,3)
!~ 		stiffness_degraded_dummy_voigt(1,6) = stiffness_degraded_dummy(1,1,2,3)
!~ 		stiffness_degraded_dummy_voigt(2,1) = stiffness_degraded_dummy(2,2,1,1)
!~ 		stiffness_degraded_dummy_voigt(2,2) = stiffness_degraded_dummy(2,2,2,2)
!~ 		stiffness_degraded_dummy_voigt(2,3) = stiffness_degraded_dummy(2,2,3,3)
!~ 		stiffness_degraded_dummy_voigt(2,4) = stiffness_degraded_dummy(2,2,1,2)
!~ 		stiffness_degraded_dummy_voigt(2,5) = stiffness_degraded_dummy(2,2,1,3)
!~ 		stiffness_degraded_dummy_voigt(2,6) = stiffness_degraded_dummy(2,2,2,3)
!~ 		stiffness_degraded_dummy_voigt(3,1) = stiffness_degraded_dummy(3,3,1,1)
!~ 		stiffness_degraded_dummy_voigt(3,2) = stiffness_degraded_dummy(3,3,2,2)
!~ 		stiffness_degraded_dummy_voigt(3,3) = stiffness_degraded_dummy(3,3,3,3)
!~ 		stiffness_degraded_dummy_voigt(3,4) = stiffness_degraded_dummy(3,3,1,2)
!~ 		stiffness_degraded_dummy_voigt(3,5) = stiffness_degraded_dummy(3,3,1,3)
!~ 		stiffness_degraded_dummy_voigt(3,6) = stiffness_degraded_dummy(3,3,2,3)
!~ 		stiffness_degraded_dummy_voigt(4,1) = stiffness_degraded_dummy(2,1,1,1)
!~ 		stiffness_degraded_dummy_voigt(4,2) = stiffness_degraded_dummy(2,1,2,2)
!~ 		stiffness_degraded_dummy_voigt(4,3) = stiffness_degraded_dummy(2,1,3,3)
!~ 		stiffness_degraded_dummy_voigt(4,4) = stiffness_degraded_dummy(2,1,1,2)
!~ 		stiffness_degraded_dummy_voigt(4,5) = stiffness_degraded_dummy(2,1,1,3)
!~ 		stiffness_degraded_dummy_voigt(4,6) = stiffness_degraded_dummy(2,1,2,3)
!~ 		stiffness_degraded_dummy_voigt(5,1) = stiffness_degraded_dummy(3,1,1,1)
!~ 		stiffness_degraded_dummy_voigt(5,2) = stiffness_degraded_dummy(3,1,2,2)
!~ 		stiffness_degraded_dummy_voigt(5,3) = stiffness_degraded_dummy(3,1,3,3)
!~ 		stiffness_degraded_dummy_voigt(5,4) = stiffness_degraded_dummy(3,1,1,2)
!~ 		stiffness_degraded_dummy_voigt(5,5) = stiffness_degraded_dummy(3,1,1,3)
!~ 		stiffness_degraded_dummy_voigt(5,6) = stiffness_degraded_dummy(3,1,2,3)
!~ 		stiffness_degraded_dummy_voigt(6,1) = stiffness_degraded_dummy(3,2,1,1)
!~ 		stiffness_degraded_dummy_voigt(6,2) = stiffness_degraded_dummy(3,2,2,2)
!~ 		stiffness_degraded_dummy_voigt(6,3) = stiffness_degraded_dummy(3,2,3,3)
!~ 		stiffness_degraded_dummy_voigt(6,4) = stiffness_degraded_dummy(3,2,1,2)
!~ 		stiffness_degraded_dummy_voigt(6,5) = stiffness_degraded_dummy(3,2,1,3)
!~ 		stiffness_degraded_dummy_voigt(6,6) = stiffness_degraded_dummy(3,2,2,3)
!~ 		WRITE(6,*), "---------------------------------------"
!~ 		WRITE(6,*) 'stiffness degraded: '
!~ 		DO i1 = 1, 6
!~ 			WRITE(6,'(7(F16.4,1X))') (stiffness_degraded_dummy_voigt(i1,i2), i2=1,6)
!~ 		END DO
!~ 		WRITE(6,*) 'phase: ', phase
!~ 		WRITE(6,*) 'degF: ', degD
!~ 		WRITE(6,*), "---------------------------------------"
!~ 	  END IF

    END FUNCTION d_d_bulkED_d_eps_d_eps

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_d_bulkED_d_eps_d_phase(eps,phase,nHFEDpar,parHFEDMatrixPhase)
    ! partial derivative of the bulk energy w.r.t. strain and phase parameter
    
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE FreeEnergyModule
      USE DegradationFunctionModule
      USE SplitEnergyModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: d_degD_d_phase
      REAL(kind=AbqRK) :: d_HFEDtens_d_eps_e(3,3), d_HFEDcomp_d_eps_e(3,3)
      DIMENSION d_d_bulkED_d_eps_d_phase(3,3)
      !
      d_HFEDtens_d_eps_e = d_HFEDpos_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      d_HFEDcomp_d_eps_e = d_HFEDneg_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      !
      d_degD_d_phase = d_degF_d_phase(phase)
      !
      d_d_bulkED_d_eps_d_phase = d_degD_d_phase*d_HFEDtens_d_eps_e
!      d_d_bulkED_d_eps_d_phase = zero

    END FUNCTION d_d_bulkED_d_eps_d_phase

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_d_bulkED_d_eps_d_phase_H(eps,phase,nHFEDpar,parHFEDMatrixPhase)
    ! partial derivative of the bulk energy w.r.t. strain and phase parameter
    !
    ! Ist äquivalent zu d_d_bulkED_d_eps_d_phase
    !
    
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE FreeEnergyModule
      USE DegradationFunctionModule
      USE SplitEnergyModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: d_degD_d_phase_H
      REAL(kind=AbqRK) :: d_HFEDtens_d_eps_e(3,3), d_HFEDcomp_d_eps_e(3,3)
      DIMENSION d_d_bulkED_d_eps_d_phase_H(3,3)
      !
      d_HFEDtens_d_eps_e = d_HFEDpos_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      d_HFEDcomp_d_eps_e = d_HFEDneg_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      !
      d_degD_d_phase_H = d_degF_d_phase(phase)
      !
      d_d_bulkED_d_eps_d_phase_H = d_degD_d_phase_H*d_HFEDtens_d_eps_e
!      d_d_bulkED_d_eps_d_phase = zero

    END FUNCTION d_d_bulkED_d_eps_d_phase_H

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_d_bulkED_d_phase_d_phase(eps,phase,nHFEDpar,parHFEDMatrixPhase)
    ! partial derivative of the bulk energy w.r.t. phase and phase
    !

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE FreeEnergyModule
      USE DegradationFunctionModule
      USE SplitEnergyModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: d_d_degD_d_phase_d_phase, HFEDtens
      !
      HFEDtens = HFEDpos(eps,nHFEDpar,parHFEDMatrixPhase)
      d_d_degD_d_phase_d_phase = d_d_degF_d_phase_d_phase(phase)
      !
      d_d_bulkED_d_phase_d_phase = d_d_degD_d_phase_d_phase * HFEDtens
!      d_d_bulkED_d_phase_d_phase = zero

    END FUNCTION d_d_bulkED_d_phase_d_phase

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_d_bulkED_d_phase_d_phase_H(eps,phase,nHFEDpar,parHFEDMatrixPhase,H)
    ! partial derivative of the bulk energy w.r.t. phase and phase
    !

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE FreeEnergyModule
      USE DegradationFunctionModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: phase,H
      REAL(kind=AbqRK) :: d_d_degD_d_phase_d_phase_H, HFEDtens
      !
      HFEDtens = H
      d_d_degD_d_phase_d_phase_H = d_d_degF_d_phase_d_phase(phase)
      !
      d_d_bulkED_d_phase_d_phase_H = d_d_degD_d_phase_d_phase_H * HFEDtens

    END FUNCTION d_d_bulkED_d_phase_d_phase_H

!------------------------------------------------------------------------------------!

END MODULE BulkEnergyModule
