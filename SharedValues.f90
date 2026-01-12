!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! tensor modules
!
! Martin Olbricht, TU Bergakademie Freiberg, 13.05.2025
!
! 13.05.2025: first implementation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE SharedValues

  USE ABQINTERFACE_PF
  

  IMPLICIT NONE

  ! Beispielvariablen (können jederzeit erweitert werden)
  PUBLIC :: Incrementnumber,Elementnumber,Integrationpointnumber
  
  INTEGER(kind=AbqIK) :: Incrementnumber
  INTEGER(kind=AbqIK) :: Elementnumber
  INTEGER(kind=AbqIK) :: Integrationpointnumber

END MODULE SharedValues
