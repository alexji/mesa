      ! declarations for imex controls
         
         integer :: &
            problem_number, nz, initial_model_number, &
            newton_iter_max, nsteps_per_call, max_calls
         logical :: &
            include_P_in_momentum_flux, stage1_only, thermal_diffusion_only
         real(dp) :: &
            initial_time, initial_dt, time_end, &
            R_inner, M_inner, L_inner, Lsurf_factor, &
            gamma, Cv, crad, TC0, TCa, TCb, &
            CFL, CFL_front, max_timestep_factor, limit_dE_div_E, &
            tol_resid_norm, tol_resid_max
         
         namelist /imex_controls/ &
            problem_number, nz, initial_model_number, &
            newton_iter_max, nsteps_per_call, max_calls, &
            include_P_in_momentum_flux, stage1_only, &
            initial_time, initial_dt, time_end, &
            R_inner, M_inner, L_inner, Lsurf_factor, &
            gamma, Cv, crad, TC0, TCa, TCb, &
            CFL, CFL_front, max_timestep_factor, limit_dE_div_E, &
            tol_resid_norm, tol_resid_max
            
      ! end of declarations
