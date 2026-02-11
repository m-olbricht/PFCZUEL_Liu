!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Split Module of the Helmholtz Free Energy
!
! Martin Olbricht, TU Bergakademie Freiberg, 09.07.2024
!
! 09.07.2024 Isotropic No Split
!
! 26.07.2024 Amor Split
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE SplitEnergyModule

  USE ABQINTERFACE_PF
  USE SharedValues
  
  IMPLICIT NONE


  CONTAINS

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION eps_e(eps)
    ! elastic strain

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      DIMENSION eps_e(3,3)
      !
      eps_e = eps
 
    END FUNCTION eps_e

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION HFEDposNoSplit(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos no split algorithm

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      eps_e_i = eps_e(eps)
      !
      ETensor=ElastTensor(E,nu)
      !
      !
      HFEDposNoSplit = one/two*DOUBLECONTRACTIONSTwoFourTwo(eps_e_i,ETensor,eps_e_i)
      
    END FUNCTION HFEDposNoSplit

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION HFEDnegNoSplit(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi neg no split algorithm

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      !
      HFEDnegNoSplit = zero
 
    END FUNCTION HFEDnegNoSplit

!------------------------------------------------------------------------------------
!
!
! 1. Ableitung
!
!
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_HFEDposNoSplit_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos no split algorithm

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      DIMENSION d_HFEDposNoSplit_d_eps_e(3,3)
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      eps_e_i = eps_e(eps)
      !
      ETensor=ElastTensor(E,nu)
      !
      !
      d_HFEDposNoSplit_d_eps_e = DOUBLECONTRACTIONTwoFour(eps_e_i,ETensor)
 
    END FUNCTION d_HFEDposNoSplit_d_eps_e

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_HFEDnegNoSplit_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi neg no split algorithm

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      DIMENSION d_HFEDnegNoSplit_d_eps_e(3,3)
      !
      d_HFEDnegNoSplit_d_eps_e = zero
 
    END FUNCTION d_HFEDnegNoSplit_d_eps_e

!------------------------------------------------------------------------------------
!
!
! 2. Ableitung
!
!
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_d_HFEDposNoSplit_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi neg no split algorithm

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      DIMENSION d_d_HFEDposNoSplit_d_eps_e_d_eps_e(3,3,3,3)
      !
      E          = parHFEDMatrixPhase(1)
      nu         = parHFEDMatrixPhase(2)
      !
      d_d_HFEDposNoSplit_d_eps_e_d_eps_e = ElastTensor(E,nu)
 
    END FUNCTION d_d_HFEDposNoSplit_d_eps_e_d_eps_e

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_d_HFEDnegNoSplit_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi neg no split algorithm

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      DIMENSION d_d_HFEDnegNoSplit_d_eps_e_d_eps_e(3,3,3,3)
      !
      d_d_HFEDnegNoSplit_d_eps_e_d_eps_e = zero
 
    END FUNCTION d_d_HFEDnegNoSplit_d_eps_e_d_eps_e

    
!------------------------------------------------------------------------------------
!
!
!
! Amor et al
!
!
!
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION bulkModulus(nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Amor et al (2009)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK) :: E, nu, lambda, mu
      !
      bulkModulus = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      lambda          = nu * E / ((one - two * nu)*(one + nu))
      mu              = half * one / (one + nu) * E
      !
      bulkModulus = lambda + two / three * mu
      
    END FUNCTION bulkModulus

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION HFEDpos_sph_Liu_Split(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Amor et al (2009)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu, K, mu
      !
      HFEDpos_sph_Liu_Split = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      HFEDpos_sph_Liu_Split = half * K * MacAulay(trace(eps)) ** two
      
    END FUNCTION HFEDpos_sph_Liu_Split
    
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION HFEDpos_dev_Liu_Split(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Amor et al (2009)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu, K, mu
      !
      HFEDpos_dev_Liu_Split = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      HFEDpos_dev_Liu_Split = mu * DOUBLECONTRACTIONTwoTwo(eps_D,eps_D)
      
    END FUNCTION HFEDpos_dev_Liu_Split

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION HFEDposAmorSplit(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Amor et al (2009)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu, K, mu
      !
      HFEDposAmorSplit = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      HFEDposAmorSplit = half * K * MacAulay(trace(eps)) ** two + mu * DOUBLECONTRACTIONTwoTwo(eps_D,eps_D)
      
    END FUNCTION HFEDposAmorSplit

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION HFEDnegAmorSplit(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi neg Amor et al (2009)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu, K, mu
      !
      HFEDnegAmorSplit = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      HFEDnegAmorSplit = half * K * MacAulay(-trace(eps)) ** two
      
    END FUNCTION HFEDnegAmorSplit

!------------------------------------------------------------------------------------
!
!
! 1. Ableitung
!
!
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_HFEDposAmorSplit_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Amor et al (2009)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3), Ident(3,3)
      REAL(kind=AbqRK) :: E, nu, K, mu
      DIMENSION d_HFEDposAmorSplit_d_eps_e(3,3)
      !
      d_HFEDposAmorSplit_d_eps_e = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      ! Von n abhängig machen wenn Funktion auf 2/3D aktualisiert werden
      Ident = Identity(3)
      !Ident(3,3) = zero ! Für 2D tests
      !
      !
      d_HFEDposAmorSplit_d_eps_e = K * MacAulay(trace(eps)) * Ident + two * mu * eps_D
      !
      !
      ! DEBUG:
!~       IF (Incrementnumber .GT. 524 .AND. Elementnumber .EQ. 33 .AND. Integrationpointnumber .EQ. 1) THEN
!~ 		WRITE(6,*) '-----------------------------------------------'
!~ 		WRITE(6,*) '----------Split Positive Amor------------------'
!~ 		WRITE(6,*) '-----------------------------------------------'
!~ 		WRITE(6,*) 'eps: ', eps
!~ 		WRITE(6,*) 'eps deviatoric: ', eps
!~ 		WRITE(6,*) 'K * MacAulay(trace(eps)) * Ident', K * MacAulay(trace(eps)) * Ident
!~ 		WRITE(6,*) 'two * mu * eps_D', two * mu * eps_D
!~ 		WRITE(6,*) '-----------------------------------------------'
!~ 	  END IF
      !
    END FUNCTION d_HFEDposAmorSplit_d_eps_e

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_HFEDnegAmorSplit_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Amor et al (2009)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3), Ident(3,3)
      REAL(kind=AbqRK) :: E, nu, K, mu
      DIMENSION d_HFEDnegAmorSplit_d_eps_e(3,3)
      !
      d_HFEDnegAmorSplit_d_eps_e = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      ! Von n abhängig machen wenn Funktion auf 2/3D aktualisiert werden
      Ident = Identity(3)
      !Ident(3,3) = zero ! Für 2D tests
      !
      !
      d_HFEDnegAmorSplit_d_eps_e = - K * MacAulay(-trace(eps)) * Ident

    END FUNCTION d_HFEDnegAmorSplit_d_eps_e

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_d_HFEDposAmorSplit_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Amor et al (2009)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      INTEGER(kind=AbqIK) :: i1,i2,i3
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3),P_sph_dummy(3,3,3,3),P_dev_dummy(3,3,3,3),Split_dummy(3,3,3,3),Split_dummy_voigt(6,6),P_sph_dummy_voigt(6,6),P_dev_dummy_voigt(6,6)
      REAL(kind=AbqRK) :: E, nu, K, mu
      DIMENSION d_d_HFEDposAmorSplit_d_eps_e_d_eps_e(3,3,3,3)
      !
      d_d_HFEDposAmorSplit_d_eps_e_d_eps_e = zero
      !
      P_dev_dummy = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      d_d_HFEDposAmorSplit_d_eps_e_d_eps_e = d_MacAulay(trace(eps)) * K * three * Projection_sph(1) + two * mu * Projection_dev(1)
      !
	  !
    END FUNCTION d_d_HFEDposAmorSplit_d_eps_e_d_eps_e

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_d_HFEDnegAmorSplit_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Amor et al (2009)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      INTEGER(kind=AbqIK) :: i1,i2,i3 ! Debug print
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3), Split_dummy(3,3,3,3), Split_dummy_voigt(6,6) ! dummy voigt
      REAL(kind=AbqRK) :: E, nu, K, mu
      DIMENSION d_d_HFEDnegAmorSplit_d_eps_e_d_eps_e(3,3,3,3)
      !
      d_d_HFEDnegAmorSplit_d_eps_e_d_eps_e = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      d_d_HFEDnegAmorSplit_d_eps_e_d_eps_e = d_MacAulay(-trace(eps)) * K * three * Projection_sph(1)
      !
      !

    END FUNCTION d_d_HFEDnegAmorSplit_d_eps_e_d_eps_e
    
!------------------------------------------------------------------------------------
!
!
!
! Liu et al MIXED MODE (2023) DOI: 10.1016/j.ijmecsci.2023.108368
!
! Mixed-Mode Formulation for Tension, Amor Split for Compression
!
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION Gceff(Hsph, Hdev, nMixedModepar,parMixedModeMatrixPhase)
    ! Psi pos Liu et al MIXED MODE (2023)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nMixedModepar
      REAL(kind=AbqRK), INTENT(IN) :: parMixedModeMatrixPhase(nMixedModepar)
      REAL(kind=AbqRK), INTENT(IN) :: Hsph, Hdev
      REAL(kind=AbqRK) :: kappa, m, GcI, GcII
      !
      m = 2.6 ! Benzegaggh - Kenane mit empirischem Exponenten m
      ! 2.6 ist der von Liu gewählte Wert aus dem Referenzpaper (Quelle 51 des Liu Papers)
      !
      !
      kappa = 1e-7 ! numerical Parameter to avoid /0
	  !
      Gceff = zero
      !
      GcI               = parMixedModeMatrixPhase(1)
      GcII              = parMixedModeMatrixPhase(2)
      !
      Gceff = GcI + (GcII-GcI) * (Hdev/(Hsph+Hdev + kappa)) ** m
      
    END FUNCTION Gceff
    
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION CorrectionFactor(Hsph, Hdev, nMixedModepar,parMixedModeMatrixPhase)
    ! Psi pos Liu et al MIXED MODE (2023)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nMixedModepar
      REAL(kind=AbqRK), INTENT(IN) :: parMixedModeMatrixPhase(nMixedModepar)
      REAL(kind=AbqRK), INTENT(IN) :: Hsph, Hdev
      REAL(kind=AbqRK) :: H, kappa, GcI, GcII, Gc
      !
	  !
      CorrectionFactor = zero
      !
      kappa = 1e-7 ! numerical Parameter to avoid /0
      !
      GcI               = parMixedModeMatrixPhase(1)
      GcII              = parMixedModeMatrixPhase(2)
      !
      Gc = Gceff(Hsph, Hdev, nMixedModepar,parMixedModeMatrixPhase)
      !
      H = (Hsph+Hdev)
      !
      !
      CorrectionFactor = H*(GcI * GcII)/(Gc*(Hdev*GcI+Hsph*GcII) + kappa)
      
    END FUNCTION CorrectionFactor
    
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION Heff(Hsph, Hdev, nMixedModepar,parMixedModeMatrixPhase)
    ! Psi pos Liu et al MIXED MODE (2023)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nMixedModepar
      REAL(kind=AbqRK), INTENT(IN) :: parMixedModeMatrixPhase(nMixedModepar)
      REAL(kind=AbqRK), INTENT(IN) :: Hsph, Hdev
      !
      REAL(kind=AbqRK) :: Gc, kappa, GcI, GcII, alpha
      !
	  !
      Heff = zero
      !
      GcI               = parMixedModeMatrixPhase(1)
      GcII              = parMixedModeMatrixPhase(2)
      !
      alpha = CorrectionFactor(Hsph, Hdev, nMixedModepar,parMixedModeMatrixPhase)
      !
      Gc = Gceff(Hsph, Hdev, nMixedModepar,parMixedModeMatrixPhase)
      !
      !
      IF (GcI .GT. zero .AND. GcII .GT. zero) THEN
		Heff = (Hsph/GcI + Hdev/GcII) * alpha * Gc
	  ELSE
	  
	  END IF
      
    END FUNCTION Heff
    
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION HFEDposLiuMixedMode(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Liu et al MIXED MODE (2023)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu, K, mu
      !
      HFEDposLiuMixedMode = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      HFEDposLiuMixedMode = half * K * MacAulay(trace(eps)) ** two + mu * DOUBLECONTRACTIONTwoTwo(eps_D,eps_D)
      
    END FUNCTION HFEDposLiuMixedMode

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION HFEDnegLiuMixedMode(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi neg Liu et al MIXED MODE (2023)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu, K, mu
      !
      HFEDnegLiuMixedMode = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      HFEDnegLiuMixedMode = half * K * MacAulay(-trace(eps)) ** two
      
    END FUNCTION HFEDnegLiuMixedMode

!------------------------------------------------------------------------------------
!
!
! 1. Ableitung
!
!
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_HFEDposLiuMixedMode_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Liu et al MIXED MODE (2023)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3), Ident(3,3)
      REAL(kind=AbqRK) :: E, nu, K, mu
      DIMENSION d_HFEDposLiuMixedMode_d_eps_e(3,3)
      !
      d_HFEDposLiuMixedMode_d_eps_e = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      ! Von n abhängig machen wenn Funktion auf 2/3D aktualisiert werden
      Ident = Identity(3)
      !Ident(3,3) = zero ! Für 2D tests
      !
      !
      d_HFEDposLiuMixedMode_d_eps_e = K * MacAulay(trace(eps)) * Ident + two * mu * eps_D
      !
      !
      ! DEBUG:
!~       IF (Incrementnumber .GT. 524 .AND. Elementnumber .EQ. 33 .AND. Integrationpointnumber .EQ. 1) THEN
!~ 		WRITE(6,*) '-----------------------------------------------'
!~ 		WRITE(6,*) '----------Split Positive Amor------------------'
!~ 		WRITE(6,*) '-----------------------------------------------'
!~ 		WRITE(6,*) 'eps: ', eps
!~ 		WRITE(6,*) 'eps deviatoric: ', eps
!~ 		WRITE(6,*) 'K * MacAulay(trace(eps)) * Ident', K * MacAulay(trace(eps)) * Ident
!~ 		WRITE(6,*) 'two * mu * eps_D', two * mu * eps_D
!~ 		WRITE(6,*) '-----------------------------------------------'
!~ 	  END IF
      !
    END FUNCTION d_HFEDposLiuMixedMode_d_eps_e

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_HFEDnegLiuMixedMode_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Liu et al MIXED MODE (2023)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3), Ident(3,3)
      REAL(kind=AbqRK) :: E, nu, K, mu
      DIMENSION d_HFEDnegLiuMixedMode_d_eps_e(3,3)
      !
      d_HFEDnegLiuMixedMode_d_eps_e = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      ! Von n abhängig machen wenn Funktion auf 2/3D aktualisiert werden
      Ident = Identity(3)
      !Ident(3,3) = zero ! Für 2D tests
      !
      !
      d_HFEDnegLiuMixedMode_d_eps_e = - K * MacAulay(-trace(eps)) * Ident

    END FUNCTION d_HFEDnegLiuMixedMode_d_eps_e

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_d_HFEDposLiuMixedMode_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Liu et al MIXED MODE (2023)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      INTEGER(kind=AbqIK) :: i1,i2,i3
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3),P_sph_dummy(3,3,3,3),P_dev_dummy(3,3,3,3),Split_dummy(3,3,3,3),Split_dummy_voigt(6,6),P_sph_dummy_voigt(6,6),P_dev_dummy_voigt(6,6)
      REAL(kind=AbqRK) :: E, nu, K, mu
      DIMENSION d_d_HFEDposLiuMixedMode_d_eps_e_d_eps_e(3,3,3,3)
      !
      d_d_HFEDposLiuMixedMode_d_eps_e_d_eps_e = zero
      !
      P_dev_dummy = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      d_d_HFEDposLiuMixedMode_d_eps_e_d_eps_e = d_MacAulay(trace(eps)) * K * three * Projection_sph(1) + two * mu * Projection_dev(1)
      !
	  !
    END FUNCTION d_d_HFEDposLiuMixedMode_d_eps_e_d_eps_e

!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION d_d_HFEDnegLiuMixedMode_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos Liu et al MIXED MODE (2023)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE MathModul
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      INTEGER(kind=AbqIK) :: i1,i2,i3 ! Debug print
      REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3), Split_dummy(3,3,3,3), Split_dummy_voigt(6,6) ! dummy voigt
      REAL(kind=AbqRK) :: E, nu, K, mu
      DIMENSION d_d_HFEDnegLiuMixedMode_d_eps_e_d_eps_e(3,3,3,3)
      !
      d_d_HFEDnegLiuMixedMode_d_eps_e_d_eps_e = zero
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      mu              = half * one / (one + nu) * E
      !
      K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
      !
      !
      eps_D = deviator(eps)
      !
      d_d_HFEDnegLiuMixedMode_d_eps_e_d_eps_e = d_MacAulay(-trace(eps)) * K * three * Projection_sph(1)
      !
      !

    END FUNCTION d_d_HFEDnegLiuMixedMode_d_eps_e_d_eps_e
    
!------------------------------------------------------------------------------
    
END MODULE SplitEnergyModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


