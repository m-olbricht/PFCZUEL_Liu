!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! degradation function module
!
! Martin Olbricht, TU Bergakademie Freiberg, 07.01.2025
!
! Degradationsfunktion standard quadratic function, 07.01.2025
!
! Degradationsfunktion nach Liu (2023), 1. Ableitung VZ Prüfen!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE DegradationFunctionInterfaceModule

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

    PURE REAL(kind=AbqRK) FUNCTION quad_degF(delta0, delta1, phase)
    ! Degradation function
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: delta0, delta1, phase
      REAL(kind=AbqRK) :: kappa
      !
      kappa = 0.000001
      !

      quad_degF =   (one-phase)**two + kappa

      IF (quad_degF .LT. kappa) THEN
        quad_degF = kappa
      END IF
      
      END FUNCTION quad_degF

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               first derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_quad_degF_d_phase(delta0, delta1, phase)
    ! Derivative Degradation function w. r. t. phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: delta0, delta1, phase
      REAL(kind=AbqRK) :: kappa
      !
      d_quad_degF_d_phase  = two*(phase-one)
      
      END FUNCTION d_quad_degF_d_phase

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               second derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_d_quad_degF_d_phase_d_phase(delta0, delta1, phase)
    ! Derivative Degradation function w. r. t. phase parameter phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: delta0, delta1, phase
      REAL(kind=AbqRK) :: kappa
      !
      d_d_quad_degF_d_phase_d_phase = two

      END FUNCTION d_d_quad_degF_d_phase_d_phase

!------------------------------------------------------------------------------------
!
!
!
!       Liu (2023 A modified phase-field model for cohesive interface failure in quasi-brittle solids)
!
!
!
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION quad_degFi_Liu(delta0, delta1, phase)
    ! Degradation function
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: delta0, delta1, phase
      REAL(kind=AbqRK) :: Chi1, Chi2, kappa
      !
      Chi1 = two * ((delta1 - delta0)*delta0**two)/delta1**three
      Chi2 = delta0**two/delta1**two
      !
      kappa = 0.000001
      !
      quad_degFi_Liu = ((delta0/((one-phase)*delta0+phase*delta1))**two + (Chi1 + Chi2)*phase**two - (Chi1 + two * Chi2) * phase) + kappa
      
	  IF (quad_degFi_Liu .LT. kappa) THEN
        quad_degFi_Liu = kappa
      END IF
      
      END FUNCTION quad_degFi_Liu

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               first derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_quad_degFi_Liu_d_phase(delta0, delta1, phase)
    ! Derivative Degradation function w. r. t. phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: delta0, delta1, phase
      REAL(kind=AbqRK) :: Chi1, Chi2, kappa
      !
      Chi1 = two * ((delta1 - delta0)*delta0**two)/delta1**three
      Chi2 = delta0**two/delta1**two
      !
      d_quad_degFi_Liu_d_phase = (-two * delta0**two*(-delta0+delta1))/((one-phase)*delta0 + phase * delta1)**three + two * (Chi1 + Chi2)*phase - (Chi1 + two * Chi2)
      
      END FUNCTION d_quad_degFi_Liu_d_phase

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               second derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_d_quad_degFi_Liu_d_phase_d_phase(delta0, delta1, phase)
    ! Derivative Degradation function w. r. t. phase parameter phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: delta0, delta1, phase
      REAL(kind=AbqRK) :: Chi1, Chi2, kappa
      !
      Chi1 = two * ((delta1 - delta0)*delta0**two)/delta1**three
      Chi2 = delta0**two/delta1**two
      !
      d_d_quad_degFi_Liu_d_phase_d_phase  = (six*delta0**two*(-delta0+delta1)**two)/((one-phase)*delta0 + phase*delta1)**four + two * (Chi1 + Chi2)

      END FUNCTION d_d_quad_degFi_Liu_d_phase_d_phase

!------------------------------------------------------------------------------------
!
!
!
!       No Degradation function (To recover the linear-Elastic)
!
!
!
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION No_degF(delta0, delta1, phase)
    ! Degradation function
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: delta0, delta1, phase
      REAL(kind=AbqRK) :: kappa
      !
      No_degF = one
      
      END FUNCTION No_degF

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               first derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_No_degF_d_phase(delta0, delta1, phase)
    ! Derivative Degradation function w. r. t. phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: delta0, delta1, phase
      REAL(kind=AbqRK) :: kappa
      !
      d_No_degF_d_phase  = zero
      
      END FUNCTION d_No_degF_d_phase

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!                               second derivatives
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_d_No_degF_d_phase_d_phase(delta0, delta1, phase)
    ! Derivative Degradation function w. r. t. phase parameter phase parameter
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: delta0, delta1, phase
      REAL(kind=AbqRK) :: kappa
      !
      d_d_No_degF_d_phase_d_phase  = zero
      
      END FUNCTION d_d_No_degF_d_phase_d_phase

!------------------------------------------------------------------------------------

END MODULE DegradationFunctionInterfaceModule

