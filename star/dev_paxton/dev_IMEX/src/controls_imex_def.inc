      ! declarations for imex controls
         
         integer :: &
            problem_number, nz, initial_model_number, &
            newton_iter_max, nsteps_per_call, max_calls, &
            model_number_for_nsteps_per_call2, nsteps_per_call2
         logical :: &
            include_P_in_momentum_flux, stage1_only, use_JFNK, &
            use_HLLC, use_T_prev_for_L, thermal_diffusion_only, &
            low_energy_density_regime, time_centering, &
            piecewise_constant_reconstruction
         real(dp) :: &
            initial_time, initial_dt, time_end, &
            R_inner, M_inner, L_inner, Lsurf_factor, &
            gamma, cgas, crad, TC0, TCa, TCb, &
            CFL, CFL_front, max_timestep_factor, limit_dE_div_E, &
            tol_resid_norm, tol_resid_max
         
         namelist /imex_controls/ &
            problem_number, nz, initial_model_number, &
            newton_iter_max, nsteps_per_call, max_calls, &
            model_number_for_nsteps_per_call2, nsteps_per_call2, &
            include_P_in_momentum_flux, stage1_only, use_JFNK, &
            use_HLLC, use_T_prev_for_L, thermal_diffusion_only, &
            low_energy_density_regime, time_centering, &
            piecewise_constant_reconstruction, &
            initial_time, initial_dt, time_end, &
            R_inner, M_inner, L_inner, Lsurf_factor, &
            gamma, cgas, crad, TC0, TCa, TCb, &
            CFL, CFL_front, max_timestep_factor, limit_dE_div_E, &
            tol_resid_norm, tol_resid_max
            
      ! end of declarations
