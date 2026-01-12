!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! viscous dissipation energy module
!
! Martin Olbricht, TU Bergakademie Freiberg, 14.03.2024
!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



MODULE ViscousDissipationModule

  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS
  
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION VDED_GL(phase,phase_old,dtime,nCSEDpar,parCSEDMatrixPhase,nVDEDpar,parVDEDMatrixPhase)
    ! viscous dissipation energy
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrixPhase(nCSEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixPhase(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: phase_old
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      REAL(kind=AbqRK) :: eta
      !
      l0               = parCSEDMatrixPhase(1)
      Gc               = parCSEDMatrixPhase(2)
      !
      eta              = parVDEDMatrixPhase(1)
      !
      VDED_GL = zero
      
      IF (dtime .NE. zero) THEN
        VDED_GL = half*eta/dtime*(phase-phase_old)**two
      END IF 
      
      END FUNCTION VDED_GL

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               first derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_VDED_GL_d_phase(phase,phase_old,dtime,nCSEDpar,parCSEDMatrixPhase,nVDEDpar,parVDEDMatrixPhase)
    ! Derivative viscous dissipation energy w. r. t. phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrixPhase(nCSEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixPhase(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: phase_old
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      REAL(kind=AbqRK) :: eta
      !
      l0               = parCSEDMatrixPhase(1)
      Gc               = parCSEDMatrixPhase(2)
      !
      eta              = parVDEDMatrixPhase(1)
      !
      d_VDED_GL_d_phase = zero
      IF (dtime .NE. zero) THEN
        d_VDED_GL_d_phase = eta/dtime*(phase-phase_old)
      END IF 
      
      END FUNCTION d_VDED_GL_d_phase

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               second derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_d_VDED_GL_d_phase_d_phase(phase,phase_old,dtime,nCSEDpar,parCSEDMatrixPhase,nVDEDpar,parVDEDMatrixPhase)
    ! Derivative viscous dissipation energy w. r. t. phase parameter phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrixPhase(nCSEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixPhase(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: phase_old
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      REAL(kind=AbqRK) :: eta
      !
      l0               = parCSEDMatrixPhase(1)
      Gc               = parCSEDMatrixPhase(2)
      !
      eta              = parVDEDMatrixPhase(1)
      !
      d_d_VDED_GL_d_phase_d_phase = zero
      IF (dtime .NE. zero) THEN
        d_d_VDED_GL_d_phase_d_phase = one/dtime*eta
      END IF 
      
      END FUNCTION d_d_VDED_GL_d_phase_d_phase

!------------------------------------------------------------------------------------
!
!
!	Liu Variation (A modified phase-field model for cohesive interface failure in quasi-brittle solids, 2023)
!
!
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION VDED_Liu(phase,phase_old,dtime,nCSEDpar,parCSEDMatrixPhase,nVDEDpar,parVDEDMatrixPhase)
    ! viscous dissipation energy
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrixPhase(nCSEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixPhase(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: phase_old
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      REAL(kind=AbqRK) :: eta
      !
      l0               = parCSEDMatrixPhase(1)
      Gc               = parCSEDMatrixPhase(2)
      !
      eta              = parVDEDMatrixPhase(1)
      !
      VDED_Liu = zero
      
      IF (dtime .NE. zero) THEN
        VDED_Liu = Gc/l0*eta/dtime*(phase-phase_old)*phase
      END IF 
      
      END FUNCTION VDED_Liu

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               first derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_VDED_Liu_d_phase(phase,phase_old,dtime,nCSEDpar,parCSEDMatrixPhase,nVDEDpar,parVDEDMatrixPhase)
    ! Derivative viscous dissipation energy w. r. t. phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrixPhase(nCSEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixPhase(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: phase_old
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      REAL(kind=AbqRK) :: eta
      !
      l0               = parCSEDMatrixPhase(1)
      Gc               = parCSEDMatrixPhase(2)
      !
      eta              = parVDEDMatrixPhase(1)
      !
      d_VDED_Liu_d_phase = zero
      IF (dtime .NE. zero) THEN
        d_VDED_Liu_d_phase = Gc/l0*eta/dtime*(phase-phase_old)
      END IF 
      
      END FUNCTION d_VDED_Liu_d_phase

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               second derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_d_VDED_Liu_d_phase_d_phase(phase,phase_old,dtime,nCSEDpar,parCSEDMatrixPhase,nVDEDpar,parVDEDMatrixPhase)
    ! Derivative viscous dissipation energy w. r. t. phase parameter phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCSEDpar
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parCSEDMatrixPhase(nCSEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixPhase(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: phase_old
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      REAL(kind=AbqRK) :: eta
      !
      l0               = parCSEDMatrixPhase(1)
      Gc               = parCSEDMatrixPhase(2)
      !
      eta              = parVDEDMatrixPhase(1)
      !
      d_d_VDED_Liu_d_phase_d_phase = zero
      IF (dtime .NE. zero) THEN
        d_d_VDED_Liu_d_phase_d_phase = zero
      END IF 
      
      END FUNCTION d_d_VDED_Liu_d_phase_d_phase

!------------------------------------------------------------------------------------

END MODULE ViscousDissipationModule
