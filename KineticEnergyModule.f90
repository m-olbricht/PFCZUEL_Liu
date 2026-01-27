!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Kinetic Energy Module
! 
! Martin Olbricht, TU Bergakademie Freiberg, 21.01.2026
!
! Functions of the kinetic Energy (not included in Kt or Fint, only for Energy)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE KineticEnergyModule

  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS
  
!------------------------------------------------------------------------------------
!
!
!
!       Kinetic Energy
!
!
!
!------------------------------------------------------------------------------------

    REAL(kind=AbqRK) FUNCTION kinED(density, vel, D)
    
    ! kinetic Energy Density
    
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      
      INTEGER(kind=AbqIK), INTENT(IN) :: D
      
      REAL(kind=AbqRK), INTENT(IN) :: density, vel(D)
      !
      kinED = zero
      !
	  !
      kinED =   half*density*SINGLECONTRACTIONOneOne(vel,vel)
      
    END FUNCTION kinED

!------------------------------------------------------------------------------------

END MODULE KineticEnergyModule
