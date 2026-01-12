!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! crack surface density functional
!
! Martin Olbricht, TU Bergakademie Freiberg, 12.01.2024
!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CrackSurfaceEnergyModule


  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS
  
  

    PURE REAL(kind=AbqRK) FUNCTION CSED(phase,grad_phase,nCSEDpar,parCSEDMatrixPhase)
    ! CrackSurfaceEnergyDensity
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrixPhase(nCSEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: grad_phase(3)
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      !
      l0               = parCSEDMatrixPhase(1)
      Gc               = parCSEDMatrixPhase(2)
      !
      CSED = one/two*Gc*(one/l0*phase**two+l0*SINGLECONTRACTIONOneOne(grad_phase, grad_phase))
      
      END FUNCTION CSED

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               first derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_CSED_d_phase(phase,grad_phase,nCSEDpar,parCSEDMatrixPhase)
    ! Derivative CrackSurfaceEnergyDensity w. r. t. phase parameter

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrixPhase(nCSEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: grad_phase(3)
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      !
      l0               = parCSEDMatrixPhase(1)
      Gc               = parCSEDMatrixPhase(2)
      !
      d_CSED_d_phase = Gc*one/l0*phase
      
      END FUNCTION d_CSED_d_phase

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_CSED_d_grad_phase(phase,grad_phase,nCSEDpar,parCSEDMatrixPhase)
    ! Derivative CrackSurfaceEnergyDensity w. r. t. grad phase parameter

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrixPhase(nCSEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: grad_phase(3)
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      DIMENSION d_CSED_d_grad_phase(3)
      !
      l0               = parCSEDMatrixPhase(1)
      Gc               = parCSEDMatrixPhase(2)
      !
      d_CSED_d_grad_phase = Gc*l0*grad_phase
      
      END FUNCTION d_CSED_d_grad_phase

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               second derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_d_CSED_d_phase_d_phase(phase,grad_phase,nCSEDpar,parCSEDMatrixPhase)
    ! Derivative CrackSurfaceEnergyDensity w. r. t. phase parameter phase parameter

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrixPhase(nCSEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: grad_phase(3)
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      !
      l0               = parCSEDMatrixPhase(1)
      Gc               = parCSEDMatrixPhase(2)
      !
      d_d_CSED_d_phase_d_phase = one/l0*Gc
      
      END FUNCTION d_d_CSED_d_phase_d_phase

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_d_CSED_d_grad_phase_d_grad_phase(phase,grad_phase,nCSEDpar,parCSEDMatrixPhase)
    ! Derivative CrackSurfaceEnergyDensity w. r. t. grad phase parameter grad phase parameter

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrixPhase(nCSEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: grad_phase(3)
      INTEGER(kind=AbqIK) :: i1
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      DIMENSION d_d_CSED_d_grad_phase_d_grad_phase(3,3)
      !
      l0               = parCSEDMatrixPhase(1)
      Gc               = parCSEDMatrixPhase(2)
      !
      d_d_CSED_d_grad_phase_d_grad_phase = zero
      !
      FORALL (i1=1:3)
        d_d_CSED_d_grad_phase_d_grad_phase(i1,i1) = Gc*l0
      END FORALL
      
      END FUNCTION d_d_CSED_d_grad_phase_d_grad_phase

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_d_CSED_d_phase_d_grad_phase(phase,grad_phase,nCSEDpar,parCSEDMatrixPhase)
    ! Derivative CrackSurfaceEnergyDensity w. r. t. phase parameter grad_phase parameter

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrixPhase(nCSEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: grad_phase(3)
      !
      DIMENSION d_d_CSED_d_phase_d_grad_phase(3)
      !
      d_d_CSED_d_phase_d_grad_phase = zero

      END FUNCTION d_d_CSED_d_phase_d_grad_phase

!------------------------------------------------------------------------------------

END MODULE CrackSurfaceEnergyModule
