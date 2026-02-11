#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#
# create user subroutine object files and shared libraries
#
# Stephan Roth, TU Bergakademie Freiberg, 25.04.2024
#
#####################################################################################################
#
# prefix
prefix="CZUELlib_PFF_CZ_240425"
#
#####################################################################################################
#
# all subroutines
#
# array of all files, consider file order!
declare -a fileArray
# general and common used modules
fileArray[1]="ABQinterface.f90"
fileArray[2]="UEL_lib_PF_240228.f90"
fileArray[3]="BasicModulesUELlib_PF_210615.f90"
fileArray[4]="TensorModule_200730.f90"
fileArray[5]="DegFunctionModule.f90"
fileArray[6]='SharedValues.f90'
# phase-field
fileArray[7]='KineticEnergyModule.f90'
fileArray[8]="PhaseField_Parameters_220316.f90"
fileArray[9]="FreeEnergyModule_PFF_240103.f90"
fileArray[10]="SplitEnergyModule.f90"
fileArray[11]="BulkEnergyModule_PFF_240409.f90"
fileArray[12]="CrackSurfaceEnergyModule_PFF.f90"
fileArray[13]="ViscousDissipationEnergyModule_PFF_240402.f90"
fileArray[14]="PhaseFieldC_modules_PFF_240514_combined.f90"
fileArray[15]="PF_UELlibmodules_PFF_240103.f90"
fileArray[16]="PFUEL_lib_PFF_240103.f90"
# cohesive zone
fileArray[17]="DegFunctionInterfaceModule.f90"
fileArray[18]="InterfaceEnergyModule.f90"
fileArray[19]="CZ_Parameters_240103.f90"
fileArray[20]="SplitIntEnergyModule.f90"
fileArray[21]="CohesiveEnergyModule_normed_PFF_240126.f90" # Choose: normed / unnormed
# fileArray[21]="CohesiveEnergyEnergyModuleMixedMode.f90"
fileArray[22]="ContactEnergyModule_PFF_240126.f90"
fileArray[23]="PenaltyEnergyModule_PFF_240109.f90"
fileArray[24]="ViscousDissipationEnergyModule_CZ.f90" # Loeschen fuer originales 25.11.2025
fileArray[25]="AliasModuleCZ.f90"
fileArray[26]="UMAT_TSL_UELlib_240126.f90"
fileArray[27]="CZUELlibmodules_240126.f90"
fileArray[28]="CZUEL_lib_240126.f90"
# main file
fileArray[29]="USER_MAIN_240103.f90"
#
# number of all files
numFiles=${#fileArray[*]} 
#
# combine all files in userfile-B.f
rm userfile-B.f
for (( filenum = 1; filenum <= $numFiles; filenum++ ))
  do
  cat ${fileArray[$filenum]} >> userfile-B.f
done
#
compilerOptions="-c -fPIC -qopenmp -extend-source -nostandard-realloc-lhs -WB -I%I"
ifort $compilerOptions userfile-B.f
mv userfile-B.o "$prefix".o
#
rm userfile-B.f *.mod
#
#####################################################################################################
