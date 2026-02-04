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

!~   USE SplitEnergyModule, ONLY: HFEDpos => HFEDposNoSplit, &
!~                                                d_HFEDpos_d_eps_e => d_HFEDposNoSplit_d_eps_e, &
!~                                                d_d_HFEDpos_d_eps_e_d_eps_e => d_d_HFEDposNoSplit_d_eps_e_d_eps_e, &
!~                                                HFEDneg => HFEDnegNoSplit, &
!~                                                d_HFEDneg_d_eps_e => d_HFEDnegNoSplit_d_eps_e, &
!~                                                d_d_HFEDneg_d_eps_e_d_eps_e => d_d_HFEDnegNoSplit_d_eps_e_d_eps_e
                                               
!~   USE SplitEnergyModule, ONLY: HFEDpos => HFEDposAmorSplit, &
!~                                                d_HFEDpos_d_eps_e => d_HFEDposAmorSplit_d_eps_e, &
!~                                                d_d_HFEDpos_d_eps_e_d_eps_e => d_d_HFEDposAmorSplit_d_eps_e_d_eps_e, &
!~                                                HFEDneg => HFEDnegAmorSplit, &
!~                                                d_HFEDneg_d_eps_e => d_HFEDnegAmorSplit_d_eps_e, &
!~                                                d_d_HFEDneg_d_eps_e_d_eps_e => d_d_HFEDnegAmorSplit_d_eps_e_d_eps_e

  USE SplitEnergyModule, ONLY: HFEDpos => HFEDposLiuMixedMode, &
                                               d_HFEDpos_d_eps_e => d_HFEDposLiuMixedMode_d_eps_e, &
                                               d_d_HFEDpos_d_eps_e_d_eps_e => d_d_HFEDposLiuMixedMode_d_eps_e_d_eps_e, &
                                               HFEDneg => HFEDnegLiuMixedMode, &
                                               d_HFEDneg_d_eps_e => d_HFEDnegLiuMixedMode_d_eps_e, &
                                               d_d_HFEDneg_d_eps_e_d_eps_e => d_d_HFEDnegLiuMixedMode_d_eps_e_d_eps_e

                                               
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

    REAL(kind=AbqRK) FUNCTION HFEDtens_H_sph(eps,nHFEDpar,parHFEDMatrixPhase)
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
      HFEDtens_H_sph = HFEDpos_sph_Liu_Split(eps,nHFEDpar,parHFEDMatrixPhase)
    END FUNCTION HFEDtens_H_sph
    
!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION HFEDtens_H_dev(eps,nHFEDpar,parHFEDMatrixPhase)
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
      HFEDtens_H_dev = HFEDpos_dev_Liu_Split(eps,nHFEDpar,parHFEDMatrixPhase)
    END FUNCTION HFEDtens_H_dev
    
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
      USE SharedValues
      
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
!~       IF (Incrementnumber .LE. 5 .AND. Elementnumber .EQ. 552 .AND. Integrationpointnumber .EQ. 1) THEN
!~ 	WRITE(6,*) 'HFEDtens: ', HFEDtens
!~ 	WRITE(6,*) 'd_degD_d_phase_H', d_degD_d_phase_H
!~ 	WRITE(6,*) 'phase: ', phase
!~       END IF
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
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: d_d_HFEDtens_d_eps_e_d_eps_e(3,3,3,3), d_d_HFEDcomp_d_eps_e_d_eps_e(3,3,3,3), stiffness_degraded_dummy(3,3,3,3), stiffness_degraded_dummy_voigt(6,6), stiffness_dummy(3,3,3,3), stiffness_dummy_voigt(6,6)! debug dummys
      REAL(kind=AbqRK) :: degD
      DIMENSION d_d_bulkED_d_eps_d_eps(3,3,3,3)
      !
      !
      d_d_HFEDtens_d_eps_e_d_eps_e = d_d_HFEDpos_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      d_d_HFEDcomp_d_eps_e_d_eps_e = d_d_HFEDneg_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      !
      degD = degF(phase)
      !
      !
      d_d_bulkED_d_eps_d_eps = degD*d_d_HFEDtens_d_eps_e_d_eps_e + d_d_HFEDcomp_d_eps_e_d_eps_e
      
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

    REAL(kind=AbqRK) FUNCTION d_d_bulkED_d_phase_d_phase_H(eps,phase,nHFEDpar,parHFEDMatrixPhase, H)
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
      REAL(kind=AbqRK), INTENT(IN) :: phase, H
      REAL(kind=AbqRK) :: d_d_degD_d_phase_d_phase_H, HFEDtens
      !
      HFEDtens = H
      d_d_degD_d_phase_d_phase_H = d_d_degF_d_phase_d_phase(phase)
      !
      d_d_bulkED_d_phase_d_phase_H = d_d_degD_d_phase_d_phase_H * HFEDtens

    END FUNCTION d_d_bulkED_d_phase_d_phase_H

!------------------------------------------------------------------------------------!

END MODULE BulkEnergyModule
