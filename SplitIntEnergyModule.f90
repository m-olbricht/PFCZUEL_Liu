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

MODULE SplitIntEnergyModule



  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS

!------------------------------------------------------------------------------------
!
!
!
! Liu et al MIXED MODE (2023) DOI: 10.1016/j.ijmecsci.2023.108368
!
! Mixed-Mode Formulation for Tension, Amor Split for Compression
!
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION Gieff(Hisph, Hidev, nMixedModeCZpar,parMixedModeCZMatrixPhase)
    ! Psi pos Liu et al MIXED MODE (2023)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nMixedModeCZpar
      REAL(kind=AbqRK), INTENT(IN) :: parMixedModeCZMatrixPhase(nMixedModeCZpar)
      REAL(kind=AbqRK), INTENT(IN) :: Hisph, Hidev
      REAL(kind=AbqRK) :: kappa, m, GiI, GiII
      !
      m = 2.6 ! Benzegaggh - Kenane mit empirischem Exponenten m
      ! 2.6 ist der von Liu gewählte Wert aus dem Referenzpaper (Quelle 51 des Liu Papers)
      !
      !
      kappa = 1e-10 ! numerical Parameter to avoid /0
	  !
      Gieff = zero
      !
      GiI               = parMixedModeCZMatrixPhase(1)
      GiII              = parMixedModeCZMatrixPhase(2)
      !
      Gieff = GiI + (GiII-GiI) * (Hidev/(Hisph+Hidev + kappa)) ** m
      
    END FUNCTION Gieff
    
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION CorrectionFactor(Hisph, Hidev, nMixedModeCZpar,parMixedModeCZMatrixPhase)
    ! Psi pos Liu et al MIXED MODE (2023)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nMixedModeCZpar
      REAL(kind=AbqRK), INTENT(IN) :: parMixedModeCZMatrixPhase(nMixedModeCZpar)
      REAL(kind=AbqRK), INTENT(IN) :: Hisph, Hidev
      REAL(kind=AbqRK) :: Hi, kappa, GiI, GiII, Gi
      !
	  !
      CorrectionFactor = zero
      !
      kappa = zero ! numerical Parameter to avoid /0
      !
      GiI               = parMixedModeCZMatrixPhase(1)
      GiII              = parMixedModeCZMatrixPhase(2)
      !
      Gi = Gieff(Hisph, Hidev, nMixedModeCZpar,parMixedModeCZMatrixPhase)
      !
      Hi = (Hisph+Hidev)
      !
      !
      CorrectionFactor = Hi*(GiI * GiII)/(Gi*(Hidev*GiI+Hisph*GiII) + kappa)
      
    END FUNCTION CorrectionFactor
    
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION Hieff(Hisph, Hidev, nMixedModeCZpar,parMixedModeCZMatrixPhase)
    ! Psi pos Liu et al MIXED MODE (2023)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nMixedModeCZpar
      REAL(kind=AbqRK), INTENT(IN) :: parMixedModeCZMatrixPhase(nMixedModeCZpar)
      REAL(kind=AbqRK), INTENT(IN) :: Hisph, Hidev
      !
      REAL(kind=AbqRK) :: Gi, kappa, GiI, GiII, alpha
      !
	  !
      Hieff = zero
      !
      GiI               = parMixedModeCZMatrixPhase(1)
      GiII              = parMixedModeCZMatrixPhase(2)
      !
      alpha = CorrectionFactor(Hisph, Hidev, nMixedModeCZpar,parMixedModeCZMatrixPhase)
      !
      Gi = Gieff(Hisph, Hidev, nMixedModeCZpar,parMixedModeCZMatrixPhase)
      !
      !
      IF (GiI .GT. zero .AND. GiII .GT. zero) THEN
		Hieff = (Hisph/GiI + Hidev/GiII) * alpha * Gi
	  ELSE
	  
	  END IF
      
    END FUNCTION Hieff
    
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION CZEDtens_H_sph(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
    ! cohesive zone energy density

      USE ABQINTERFACE
      USE FLOATNUMBERS
	  USE MathModul
	  
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      REAL(kind=AbqRK) :: CZEDsph, CZEDtang
      !
      REAL(kind=AbqRK) :: kn, kt, sn, st
      kn               = parCEDMatrixPhase(1)
      kt               = parCEDMatrixPhase(2)
      !
      sn = MacAulay(sep(normalDirectionIndex))
      st = sep(tangentialDirectionIndex)
      !
      !
      CZEDtens_H_sph = half*kn*sn**two
      
    END FUNCTION CZEDtens_H_sph

!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION CZEDtens_H_dev(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
    ! cohesive zone energy density

      USE ABQINTERFACE
      USE FLOATNUMBERS
      USE MathModul

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      !
      REAL(kind=AbqRK) :: kn, kt, sn, st
      kn               = parCEDMatrixPhase(1)
      kt               = parCEDMatrixPhase(2)
      !
      st = sep(tangentialDirectionIndex)
      !
      CZEDtens_H_dev = half*kt*st**two
      
    END FUNCTION CZEDtens_H_dev

!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!
!                               CZED cohesive zone energy density
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION CZEDposLiuMixedMode(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
    ! cohesive zone energy density

      USE ABQINTERFACE
      USE FLOATNUMBERS
      USE MathModul

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      REAL(kind=AbqRK) :: CZEDsph, CZEDtang
      !
      REAL(kind=AbqRK) :: kn, kt, sn, st
      kn               = parCEDMatrixPhase(1)
      kt               = parCEDMatrixPhase(2)
      !
      sn = MacAulay(sep(normalDirectionIndex))
      st = sep(tangentialDirectionIndex)
      !
      CZEDsph = half*kn*sn**two
      CZEDtang = half*kt*st**two
      !
      CZEDposLiuMixedMode = CZEDsph + CZEDtang
      
    END FUNCTION CZEDposLiuMixedMode

!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION CZEDnegLiuMixedMode(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
    ! cohesive zone energy density

      USE ABQINTERFACE
      USE FLOATNUMBERS
      USE MathModul

      IMPLICIT NONE

      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      REAL(kind=AbqRK) :: CZEDsph
      !
      REAL(kind=AbqRK) :: kn, kt, sn, st
      kn               = parCEDMatrixPhase(1)
      !
      sn = MacAulay_neg(sep(normalDirectionIndex))
      !
      CZEDsph = half*kn*sn**two
      !
      CZEDnegLiuMixedMode = CZEDsph
      
    END FUNCTION CZEDnegLiuMixedMode

!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!
!                               first derivatives
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION d_CZEDposLiuMixedMode_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
    ! derivative of bulk energy density w.r.t. strain

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      REAL(kind=AbqRK) :: d_CZEDpos_d_sep(D)
      DIMENSION d_CZEDposLiuMixedMode_d_sep(D)
      !
      REAL(kind=AbqRK) :: kn, kt, sn, st
      kn               = parCEDMatrixPhase(1)
      kt               = parCEDMatrixPhase(2)
      !
      sn = MacAulay(sep(normalDirectionIndex))
      st = sep(tangentialDirectionIndex)
      !
      !
      d_CZEDpos_d_sep(normalDirectionIndex) = kn*sn
      d_CZEDpos_d_sep(tangentialDirectionIndex) = kt*st
      d_CZEDposLiuMixedMode_d_sep = d_CZEDpos_d_sep
      

    END FUNCTION d_CZEDposLiuMixedMode_d_sep

!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION d_CZEDnegLiuMixedMode_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
    ! derivative of bulk energy density w.r.t. strain

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      !
      !
      REAL(kind=AbqRK) :: kn, kt, sn, st
      DIMENSION d_CZEDnegLiuMixedMode_d_sep(D)
      !
      kn               = parCEDMatrixPhase(1)
      kt               = parCEDMatrixPhase(2)
      !
      sn = MacAulay_neg(sep(normalDirectionIndex))
      !
      d_CZEDnegLiuMixedMode_d_sep = zero
      !
      d_CZEDnegLiuMixedMode_d_sep(normalDirectionIndex) = kn*sn

    END FUNCTION d_CZEDnegLiuMixedMode_d_sep

!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!
!                               second derivatives
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION d_CZEDposLiuMixedMode_d_sep_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
    ! derivative of bulk energy density w.r.t. strain

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      !
      REAL(kind=AbqRK) :: en(D), et(D)
      REAL(kind=AbqRK) :: d_CZEDpos_d_sep_d_sep(D,D)
      REAL(kind=AbqRK) :: d_CZEDtang_d_sep_d_sep(D,D)
      !
      REAL(kind=AbqRK) :: kn, kt, sn, st
      DIMENSION d_CZEDposLiuMixedMode_d_sep_d_sep(D,D)
      !
      kn               = parCEDMatrixPhase(1)
      kt               = parCEDMatrixPhase(2)
      !
      d_CZEDpos_d_sep_d_sep = zero
      d_CZEDtang_d_sep_d_sep = zero
      !
      en = zero
      et = zero
      en(normalDirectionIndex) = one
      et(tangentialDirectionIndex) = one
      !
      ! IF cond. as substitute for Heaviside-Func
      IF (sep(normalDirectionIndex) .GE. zero) THEN
		d_CZEDpos_d_sep_d_sep = kn*DYADE(en,en)
	  END IF
	  d_CZEDtang_d_sep_d_sep = kt*DYADE(et,et)
	  !
	  d_CZEDposLiuMixedMode_d_sep_d_sep = d_CZEDpos_d_sep_d_sep + d_CZEDtang_d_sep_d_sep
      

    END FUNCTION d_CZEDposLiuMixedMode_d_sep_d_sep

!------------------------------------------------------------------------------------!

    REAL(kind=AbqRK) FUNCTION d_CZEDnegLiuMixedMode_d_sep_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
    ! derivative of bulk energy density w.r.t. strain

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      !
      REAL(kind=AbqRK) :: en(D)
      REAL(kind=AbqRK) :: d_CZEDneg_d_sep_d_sep(D,D)
      !
      REAL(kind=AbqRK) :: kn, kt, sn, st
      DIMENSION d_CZEDnegLiuMixedMode_d_sep_d_sep(D,D)
      !
      kn               = parCEDMatrixPhase(1)
      kt               = parCEDMatrixPhase(2)
      !
      !
      en = zero
      en(normalDirectionIndex) = one
      !
      ! IF cond. as substitute for Heaviside-Func
      IF (sep(normalDirectionIndex) .LT. zero) THEN
		d_CZEDneg_d_sep_d_sep = kn * DYADE(en,en)
	  ELSE 
		d_CZEDneg_d_sep_d_sep = zero
	  END IF
      !
	  d_CZEDnegLiuMixedMode_d_sep_d_sep = d_CZEDneg_d_sep_d_sep
	  
	  
    END FUNCTION d_CZEDnegLiuMixedMode_d_sep_d_sep

!~ !------------------------------------------------------------------------------------!

!~     REAL(kind=AbqRK) FUNCTION d_CZEDtangLiuMixedMode_d_sep_d_sep(sep, D, normalDirectionIndex, tangentialDirectionIndex, nCEDpar, parCEDMatrixPhase)
!~     ! derivative of bulk energy density w.r.t. strain

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nCEDpar, D, normalDirectionIndex, tangentialDirectionIndex
!~       REAL(kind=AbqRK), INTENT(IN) :: parCEDMatrixPhase(nCEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: sep(D)
!~       !
!~       REAL(kind=AbqRK) :: et(D)
!~       REAL(kind=AbqRK) :: d_CZEDtangLiuMixedMode_d_sep_d_sep(D,D)
!~       !
!~       REAL(kind=AbqRK) :: kn, kt, sn, st
!~       kn               = parCEDMatrixPhase(1)
!~       kt               = parCEDMatrixPhase(2)
!~       !
!~       d_CZEDtangLiuMixedMode_d_sep_d_sep = zero
!~       !
!~       et = zero
!~       et(tangentialDirectionIndex) = one
!~       !
!~       d_CZEDtangLiuMixedMode_d_sep_d_sep = kt*DYADE(et,et)

!~     END FUNCTION d_CZEDtangLiuMixedMode_d_sep_d_sep

!~ !------------------------------------------------------------------------------------!

END MODULE SplitIntEnergyModule
