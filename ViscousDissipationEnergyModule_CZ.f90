!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! viscous dissipation energy module Regularization for CZ
!
! Martin Olbricht, TU Bergakademie Freiberg, 28.11.2025
!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



MODULE ViscousDissipationCZModule

  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS
  
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION VDED_GL(phase,phase_old,dtime,nIEDpar,parIEDMatrixCZ,nVDEDpar,parVDEDMatrixCZ)
    ! viscous dissipation energy
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nIEDpar
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parIEDMatrixCZ(nIEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixCZ(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: phase_old
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      REAL(kind=AbqRK) :: eta
      !
      l0               = parIEDMatrixCZ(2)
      Gc               = parIEDMatrixCZ(1)
      !
      eta              = parVDEDMatrixCZ(1)
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

    PURE REAL(kind=AbqRK) FUNCTION d_VDED_GL_d_phase(phase,phase_old,dtime,nIEDpar,parIEDMatrixCZ,nVDEDpar,parVDEDMatrixCZ)
    ! Derivative viscous dissipation energy w. r. t. phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nIEDpar
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parIEDMatrixCZ(nIEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixCZ(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: phase_old
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      REAL(kind=AbqRK) :: eta
      !
      l0               = parIEDMatrixCZ(2)
      Gc               = parIEDMatrixCZ(1)
      !
      eta              = parVDEDMatrixCZ(1)
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

    PURE REAL(kind=AbqRK) FUNCTION d_d_VDED_GL_d_phase_d_phase(phase,phase_old,dtime,nIEDpar,parIEDMatrixCZ,nVDEDpar,parVDEDMatrixCZ)
    ! Derivative viscous dissipation energy w. r. t. phase parameter phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nIEDpar
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parIEDMatrixCZ(nIEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixCZ(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: phase_old
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      REAL(kind=AbqRK) :: eta
      !
      l0               = parIEDMatrixCZ(2)
      Gc               = parIEDMatrixCZ(1)
      !
      eta              = parVDEDMatrixCZ(1)
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

    REAL(kind=AbqRK) FUNCTION VDED_Liu(phase,phase_old,dtime,nIEDpar,parIEDMatrixCZ,nVDEDpar,parVDEDMatrixCZ)
    ! viscous dissipation energy
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      USE SharedValues
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nIEDpar
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parIEDMatrixCZ(nIEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixCZ(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: phase_old
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      REAL(kind=AbqRK) :: eta
      !
      l0               = parIEDMatrixCZ(2)
      Gc               = parIEDMatrixCZ(1)
      !
      eta              = parVDEDMatrixCZ(1)
      !
      VDED_Liu = zero
      
      IF (dtime .NE. zero) THEN
        VDED_Liu = Gc/l0*eta/dtime*(phase-phase_old)*phase
      END IF 
      
!~       IF (Incrementnumber .EQ. 3 .AND. Elementnumber .EQ. 121 .AND. Integrationpointnumber .EQ. 1) THEN
!~ 		WRITE(6,*) 'VISCOUS ROUTINE: '
!~ 		WRITE(6,*) 'Gc: ', Gc
!~ 		WRITE(6,*) 'l0: ', l0
!~ 		WRITE(6,*) 'eta: ', eta
!~ 		WRITE(6,*) 'VDED: ', VDED_Liu
!~       END IF
      
      END FUNCTION VDED_Liu

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               first derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_VDED_Liu_d_phase(phase,phase_old,dtime,nIEDpar,parIEDMatrixCZ,nVDEDpar,parVDEDMatrixCZ)
    ! Derivative viscous dissipation energy w. r. t. phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nIEDpar
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parIEDMatrixCZ(nIEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixCZ(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: phase_old
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      REAL(kind=AbqRK) :: eta
      !
      l0               = parIEDMatrixCZ(2)
      Gc               = parIEDMatrixCZ(1)
      !
      eta              = parVDEDMatrixCZ(1)
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

    PURE REAL(kind=AbqRK) FUNCTION d_d_VDED_Liu_d_phase_d_phase(phase,phase_old,dtime,nIEDpar,parIEDMatrixCZ,nVDEDpar,parVDEDMatrixCZ)
    ! Derivative viscous dissipation energy w. r. t. phase parameter phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nIEDpar
      INTEGER(kind=AbqIK), INTENT(IN) :: nVDEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parIEDMatrixCZ(nIEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: parVDEDMatrixCZ(nVDEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK), INTENT(IN) :: phase_old
      REAL(kind=AbqRK), INTENT(IN) :: dtime
      REAL(kind=AbqRK) :: l0
      REAL(kind=AbqRK) :: Gc
      REAL(kind=AbqRK) :: eta
      !
      l0               = parIEDMatrixCZ(2)
      Gc               = parIEDMatrixCZ(1)
      !
      eta              = parVDEDMatrixCZ(1)
      !
      d_d_VDED_Liu_d_phase_d_phase = zero
      IF (dtime .NE. zero) THEN
        d_d_VDED_Liu_d_phase_d_phase = zero
      END IF 
      
      END FUNCTION d_d_VDED_Liu_d_phase_d_phase

!------------------------------------------------------------------------------------

END MODULE ViscousDissipationCZModule
