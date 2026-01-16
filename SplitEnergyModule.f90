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
    
!------------------------------------------------------------------------------
    
END MODULE SplitEnergyModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


