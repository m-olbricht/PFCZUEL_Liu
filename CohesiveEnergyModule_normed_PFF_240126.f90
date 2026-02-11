!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! cohesive energy module
!
! Stephan Roth, TU Bergakademie Freiberg, 08.01.2024
!
! Martin Olbricht, Adjusted, to test PFCZ
!
! 08.01.2024: linear cohesive traction-separation law
! 12.01.2024: cubic degradation function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CohesiveEnergyModule

  !MacAulay
  USE MathModul
  USE SharedValues, ONLY: Incrementnumber,Elementnumber,Integrationpointnumber
  !
  USE DegradationFunctionInterfaceModule, ONLY: Degradation => quad_degF, &
                                              d_Degradation_d_damage => d_quad_degF_d_phase, &
                                              d_Degradation_d_damage_d_damage => d_d_quad_degF_d_phase_d_phase
!~   USE DegradationFunctionInterfaceModule, ONLY: Degradation => quad_degFi_Liu, &
!~                                               d_Degradation_d_damage => d_quad_degFi_Liu_d_phase, &
!~                                               d_Degradation_d_damage_d_damage => d_d_quad_degFi_Liu_d_phase_d_phase
!~   USE DegradationFunctionInterfaceModule, ONLY: Degradation => No_degF, &
!~                                            d_Degradation_d_damage => d_No_degF_d_phase, &
!~                                            d_Degradation_d_damage_d_damage => d_d_No_degF_d_phase_d_phase    

  USE SplitIntEnergyModule, ONLY: CZEDpos => CZEDposLiuMixedMode, &
                                               d_CZEDpos_d_sep => d_CZEDposLiuMixedMode_d_sep, &
                                               d_CZEDpos_d_sep_d_sep => d_CZEDposLiuMixedMode_d_sep_d_sep, &
                                               CZEDneg => CZEDnegLiuMixedMode, &
                                               d_CZEDneg_d_sep => d_CZEDnegLiuMixedMode_d_sep, &
                                               d_CZEDneg_d_sep_d_sep => d_CZEDnegLiuMixedMode_d_sep_d_sep


  IMPLICIT NONE

  CONTAINS

!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION CZED(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase,damage,prop_df_one,prop_df_two)
    ! cohesive zone energy density

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      REAL(kind=AbqRK), INTENT(IN) :: damage, prop_df_one, prop_df_two
      REAL(kind=AbqRK) :: CZEDtens, CZEDcomp, degi
      !
      CZEDtens = CZEDpos(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
      !
      CZEDcomp = CZEDneg(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
      !
      degi = Degradation(prop_df_one,prop_df_two,damage)
      !
      CZED = degi*CZEDtens + CZEDcomp

    END FUNCTION CZED

!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!
!                               first derivatives
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION d_CZED_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase,damage,prop_df_one,prop_df_two)
    ! cohesive zone energy density: first derivative w.r.t. separation coordinates

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      REAL(kind=AbqRK), INTENT(IN) :: damage, prop_df_one, prop_df_two
      REAL(kind=AbqRK) :: d_CZEDtens_d_sep(D),d_CZEDcomp_d_sep(D)
      REAL(kind=AbqRK) :: degi
      DIMENSION d_CZED_d_sep(D)
      !
      d_CZEDtens_d_sep = d_CZEDpos_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
      !
      d_CZEDcomp_d_sep = d_CZEDneg_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
      !
      degi = Degradation(prop_df_one,prop_df_two,damage)
      !
      d_CZED_d_sep = degi*d_CZEDtens_d_sep + d_CZEDcomp_d_sep

    END FUNCTION d_CZED_d_sep

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_CZED_d_damage(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase,damage,prop_df_one,prop_df_two)
    ! cohesive zone energy density: mixed second derivative w.r.t. separation coordinates and damage

      USE ABQINTERFACE
      USE FLOATNUMBERS
      USE SplitIntEnergyModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      REAL(kind=AbqRK), INTENT(IN) :: damage, prop_df_one, prop_df_two
      REAL(kind=AbqRK) :: CZEDtens, d_degi_d_damage
      REAL(kind=AbqRK) :: degi
	  !
	  CZEDtens = CZEDpos(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
      !
      d_degi_d_damage = d_Degradation_d_damage(prop_df_one,prop_df_two,damage)
      !
      d_CZED_d_damage = CZEDtens * d_degi_d_damage
	  

    END FUNCTION d_CZED_d_damage
    
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_CZED_d_damage_Hi(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase,damage,prop_df_one,prop_df_two, Hi)
    ! cohesive zone energy density: mixed second derivative w.r.t. separation coordinates and damage

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      REAL(kind=AbqRK), INTENT(IN) :: Hi
      REAL(kind=AbqRK), INTENT(IN) :: damage, prop_df_one, prop_df_two
      REAL(kind=AbqRK) :: degi
	  !
      !
      d_CZED_d_damage_Hi = Hi * d_Degradation_d_damage(prop_df_one,prop_df_two,damage)
	  

    END FUNCTION d_CZED_d_damage_Hi

!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!
!                               second derivatives
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION d_CZED_d_sep_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase,damage,prop_df_one,prop_df_two)
    ! cohesive zone energy density: second derivative w.r.t. separation coordinates

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      REAL(kind=AbqRK), INTENT(IN) :: damage, prop_df_one, prop_df_two
      REAL(kind=AbqRK) :: d_CZEDtens_d_sep_d_sep(D,D), d_CZEDcomp_d_sep_d_sep(D,D)
      REAL(kind=AbqRK) :: degi
      DIMENSION d_CZED_d_sep_d_sep(D,D)
      !
      d_CZED_d_sep_d_sep = zero
      !
      d_CZEDtens_d_sep_d_sep = d_CZEDpos_d_sep_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
      !
      d_CZEDcomp_d_sep_d_sep = d_CZEDneg_d_sep_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
      !
      degi = Degradation(prop_df_one,prop_df_two,damage)
      !
      d_CZED_d_sep_d_sep = degi*d_CZEDtens_d_sep_d_sep + d_CZEDcomp_d_sep_d_sep

    END FUNCTION d_CZED_d_sep_d_sep


!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_CZED_d_sep_d_damage(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase,damage,prop_df_one,prop_df_two)
    ! partial derivative of the bulk energy w.r.t. strain and phase parameter
    
      USE ABQINTERFACE
      USE FLOATNUMBERS
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      REAL(kind=AbqRK), INTENT(IN) :: damage, prop_df_one, prop_df_two
      REAL(kind=AbqRK) :: d_CZEDtens_d_sep(D),d_CZEDcomp_d_sep(D)
      REAL(kind=AbqRK) :: d_degi_d_damage
      DIMENSION d_CZED_d_sep_d_damage(D)
      !
      d_CZEDtens_d_sep = d_CZEDpos_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
      !
      d_CZEDcomp_d_sep = d_CZEDneg_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
      !
      d_degi_d_damage = Degradation(prop_df_one,prop_df_two,damage)
      !
      d_CZED_d_sep_d_damage = d_degi_d_damage*d_CZEDtens_d_sep + d_CZEDcomp_d_sep

    END FUNCTION d_CZED_d_sep_d_damage

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_CZED_d_damage_d_damage(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase,damage,prop_df_one,prop_df_two)
    ! cohesive zone energy density: mixed second derivative w.r.t. separation coordinates and damage

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      REAL(kind=AbqRK), INTENT(IN) :: damage, prop_df_one, prop_df_two
      REAL(kind=AbqRK) :: CZEDtens, d_degi_d_damage_d_damage
	  !
	  CZEDtens = CZEDpos(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
      !
      d_degi_d_damage_d_damage = d_Degradation_d_damage_d_damage(prop_df_one,prop_df_two,damage)
      !
      d_CZED_d_damage_d_damage = CZEDtens * d_degi_d_damage_d_damage
	  

    END FUNCTION d_CZED_d_damage_d_damage
    
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_CZED_d_damage_d_damage_Hi(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase,damage,prop_df_one,prop_df_two, Hi)
    ! cohesive zone energy density: mixed second derivative w.r.t. separation coordinates and damage

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      REAL(kind=AbqRK), INTENT(IN) :: Hi
      REAL(kind=AbqRK), INTENT(IN) :: damage, prop_df_one, prop_df_two
	  !
      !
      d_CZED_d_damage_d_damage_Hi = Hi * d_Degradation_d_damage_d_damage(prop_df_one,prop_df_two,damage)
	  

    END FUNCTION d_CZED_d_damage_d_damage_Hi

!------------------------------------------------------------------------------------

END MODULE CohesiveEnergyModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

