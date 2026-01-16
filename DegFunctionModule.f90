!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! degradation function module
!
! Martin Olbricht, TU Bergakademie Freiberg, 12.01.2024
!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE DegradationFunctionModule

  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS
  
!------------------------------------------------------------------------------------
!
!
!
!       Quadratic
!
!
!
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION quad_degF(phase)
    ! Degradation function
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: kappa
      !
      kappa = 0.0000001
      !

      quad_degF =   (one-phase)**two + kappa

!~       IF (quad_degF .LT. kappa) THEN
!~         quad_degF = kappa
!~       END IF
      
      END FUNCTION quad_degF

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               first derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_quad_degF_d_phase(phase)
    ! Derivative Degradation function w. r. t. phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: kappa
      !
      d_quad_degF_d_phase  = two*(phase-one)
      
      END FUNCTION d_quad_degF_d_phase

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               second derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_d_quad_degF_d_phase_d_phase(phase)
    ! Derivative Degradation function w. r. t. phase parameter phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: kappa
      !
      d_d_quad_degF_d_phase_d_phase = two

      END FUNCTION d_d_quad_degF_d_phase_d_phase

!------------------------------------------------------------------------------------
!
!
!
!       No Degradation function (To recover the linear-Elastic)
!
!
!
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION No_degF(phase)
    ! Degradation function
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: kappa
      !
      No_degF = one
      
      END FUNCTION No_degF

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               first derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_No_degF_d_phase(phase)
    ! Derivative Degradation function w. r. t. phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: kappa
      !
      d_No_degF_d_phase  = zero
      
      END FUNCTION d_No_degF_d_phase

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               second derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_d_No_degF_d_phase_d_phase(phase)
    ! Derivative Degradation function w. r. t. phase parameter phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: phase
      REAL(kind=AbqRK) :: kappa
      !
      d_d_No_degF_d_phase_d_phase  = zero
      
      END FUNCTION d_d_No_degF_d_phase_d_phase

!------------------------------------------------------------------------------------

END MODULE DegradationFunctionModule
