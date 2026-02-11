MODULE AliasModuleCZ

!  Alias Modul
!  
!  Speichert alle Umbenennungen, damit sie von verschiedenen Subroutinen aufgerufen werden können

   ! ==========================================================
   ! Penalty Energy
   ! ==========================================================
   USE PenaltyEnergyModule, ONLY: &
        PenED  => PenED_Exponential, &
        d_PenED_d_damageJump => d_PenED_d_damageJump_Exponential, &
        d_PenED_d_damageJump_d_damageJump => &
        d_PenED_d_damageJump_d_damageJump_Exponential


   ! ==========================================================
   ! Interface Energy (AT2 aktiv)
   ! Zum Wechseln auf AT1 einfach hier umstellen
   ! ==========================================================
   USE InterfaceEnergyModule, ONLY: &
        IED  => IED_AT2, &
        d_IED_d_phase => d_IED_AT2_d_phase, &
        d_IED_d_phase_d_phase => d_IED_AT2_d_phase_d_phase, &
        d_IED_d_grad_phase => d_IED_AT2_d_grad_phase, &
        d_IED_d_grad_phase_d_grad_phase => &
        d_IED_AT2_d_grad_phase_d_grad_phase, &
        d_IED_d_phase_d_grad_phase => &
        d_IED_AT2_d_phase_d_grad_phase

!~         IED  => IED_AT1, &
!~         d_IED_d_phase => d_IED_AT1_d_phase, &
!~         d_IED_d_phase_d_phase => d_IED_AT1_d_phase_d_phase, &
!~         d_IED_d_grad_phase => d_IED_AT1_d_grad_phase, &
!~         d_IED_d_grad_phase_d_grad_phase => &
!~         d_IED_AT1_d_grad_phase_d_grad_phase, &
!~         d_IED_d_phase_d_grad_phase => &
!~         d_IED_AT1_d_phase_d_grad_phase

   ! ==========================================================
   ! Viscous Dissipation CZ
   ! ==========================================================
   USE ViscousDissipationCZModule, ONLY: &
        VDED  => VDED_Liu, &
        d_VDED_d_phase => d_VDED_Liu_d_phase, &
        d_VDED_d_phase_d_phase => d_d_VDED_Liu_d_phase_d_phase


   ! ==========================================================
   ! Degradation Function
   ! ==========================================================
   USE DegradationFunctionInterfaceModule, ONLY: &
        degF  => quad_degF, &
        d_degF_d_phase => d_quad_degF_d_phase, &
        d_d_degF_d_phase_d_phase => &
        d_d_quad_degF_d_phase_d_phase
        
        
   IMPLICIT NONE
   
   PUBLIC :: PenED, d_PenED_d_damageJump, d_PenED_d_damageJump_d_damageJump
   PUBLIC :: IED, d_IED_d_phase, d_IED_d_phase_d_phase
   PUBLIC :: VDED, d_VDED_d_phase, d_VDED_d_phase_d_phase
   PUBLIC :: degF, d_degF_d_phase, d_d_degF_d_phase_d_phase
   

END MODULE AliasModuleCZ
