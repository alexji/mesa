 
   module imex_input
      use const_def, only: dp
      include 'controls_imex_def.inc' ! declarations only
      real(dp), dimension(:), allocatable :: &
         r_init, rho_init, mom_init, etot_init ! (nz)         
   end module imex_input 
 
   module imex_work
      use const_def, only: dp
      integer :: model_number, stage, newton_iter
      real(dp) :: resid_norm, tol_resid_norm, &
         time, dt, dt_next, gam, delta, Rgas
      real(dp), dimension(:), allocatable :: &
         sub, diag, sup, bp, vp, xp, &
         r, area, Vol, dr, dr_bar, dVol, & 
         Eeos_x, rhs, deltaT, P_face, &
         T_start, T_1, d_ETOT_dt_I1, d_ETOT_dt_I2
      real(dp), dimension(:,:), allocatable :: &
         prim_start, cons_start, prim_1, cons_1, &
         grad_cell, flux_face, d_cons_dt_X0, d_cons_dt_X1
   end module imex_work

   module imex_output
      use const_def, only: dp
      real(dp), dimension(:,:), allocatable :: prim, cons
      real(dp), dimension(:), allocatable :: T, L, v_face, csound
   end module imex_output
 
 
   module imex
      use imex_input

      ! bates, knoll, rider, lowrie, mousseau,
      ! On Consistent Time-Integration Methods for Radiation Hydrodynamics in
      ! the Equilibrium Diffusion Limit: Low-Energy-Density Regime
      ! Journal of Computational Physics 167, 99–130 (2001)

      use const_def, only: dp, qp, clight, pi
      use math_lib
      use auto_diff
      use utils_lib, only: is_bad
      
      implicit none

      integer, parameter :: &
         i_rho = 1, i_mom = 2, i_etot = 3, nvar = i_etot         
         
         
      contains
      
      
!!! top level control   
   
      
      subroutine start_imex(total_energy_initial, ierr)
         use imex_work
         use imex_output
         real(dp), intent(out) :: total_energy_initial
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         write(*,*) 'start_imex'
         select case(problem_number)
         case (0)
            call initialize_static_problem()
         case (1)
            call initialize_Barenblatt_problem()
         case (2)
            call initialize_smooth_problem()
         case default
            ierr = -1
            write(*,*) 'bad problem_number for imex', problem_number
            return
         end select
         
         Rgas = Cv*(gamma - 1d0)
         model_number = initial_model_number
         time = initial_time
         dt_next = initial_dt
         
         call alloc_work_arrays()
         call set_grid_vars()
         call set_init_prim_cons_T()
         total_energy_initial = sum(cons(i_etot,1:nz))
         write(*,1) 'total_energy_initial', total_energy_initial
         !stop 'start_imex'
                  
         contains
         
         subroutine set_init_prim_cons_T()
            integer :: k, j
            include 'formats'
            do k=1,nz
               T(k) = get_T(rho_init(k), mom_init(k), etot_init(k))
               prim(i_rho,k) = rho_init(k)
               prim(i_mom,k) = mom_init(k)
               prim(i_etot,k) = etot_init(k)
               do j=1,nvar
                  cons(j,k) = prim(j,k)*dVol(k)
               end do
               !write(*,2) 'r T rho etot', k, r(k), T(k), rho_init(k), etot_init(k)
            end do
         end subroutine set_init_prim_cons_T
         
         subroutine alloc_work_arrays()
            allocate(&
               r(nz), area(nz), Vol(nz), dr(nz), dr_bar(nz), dVol(nz), & 
               Eeos_x(nz), rhs(nz), sub(nz), diag(nz), sup(nz), deltaT(nz), &
               bp(nz), vp(nz), xp(nz), P_face(nz), &
               T_start(nz), T_1(nz), T(nz), L(nz), v_face(nz), csound(nz), &
               d_ETOT_dt_I1(nz), d_ETOT_dt_I2(nz), &
               prim_start(nvar,nz), cons_start(nvar,nz), &
               prim_1(nvar,nz), cons_1(nvar,nz), &
               prim(nvar,nz), cons(nvar,nz), &
               grad_cell(nvar,nz), flux_face(nvar,nz), &
               d_cons_dt_X0(nvar,nz), d_cons_dt_X1(nvar,nz))            
         end subroutine alloc_work_arrays
      
         subroutine set_grid_vars()
            integer :: k
            real(dp) :: rm1, r00, volm1, vol00
            include 'formats'
            rm1 = R_inner
            volm1 = 4d0*pi/3d0*pow3(rm1)
            do k=nz,1,-1
               r(k) = r_init(k)
               r00 = r(k)
               dr(k) = r00 - rm1
               area(k) = 4d0*pi*pow2(r00)
               vol00 = area(k)*r00/3d0
               Vol(k) = vol00
               dVol(k) = vol00 - volm1
               rm1 = r00
               volm1 = vol00
            end do
            do k=2,nz
               dr_bar(k) = 0.5d0*(dr(k-1) + dr(k))
            end do
            dr_bar(1) = 0.5d0*dr(1) ! NOT USED
         end subroutine set_grid_vars
         
      end subroutine start_imex
      
      
      subroutine finish_imex(ierr)
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         write(*,*) 'finish_imex'
      end subroutine finish_imex
      
      
      subroutine get_imex_data( &
            r_out, rho_out, E_out, P_out, v_out, T_out, L_out, cs_out, v_div_cs_out)
         use imex_work
         use imex_output
         real(dp), intent(out), dimension(:) :: &
            r_out, rho_out, E_out, P_out, v_out, T_out, L_out, cs_out, v_div_cs_out
         integer :: k
         real(dp) :: rho, Pgas, Prad, mom, v
         include 'formats'
         do k=1,nz
            rho = prim(i_rho,k)
            mom = prim(i_mom,k)
            v = mom/rho
            r_out(k) = r(k)
            rho_out(k) = rho
            E_out(k) = prim(i_etot,k)
            Pgas = (gamma-1d0)*rho*Cv*T(k)
            Prad = crad*pow4(T(k))/3
            P_out(k) = Pgas + Prad
            v_out(k) = v
            T_out(k) = T(k)
            L_out(k) = L(k)
            cs_out(k) = csound(k)
            v_div_cs_out(k) = v_out(k)/cs_out(k)
         end do
      end subroutine get_imex_data

      
      subroutine steps_imex( &
            max_steps_for_this_call, age, timestep, total_energy, &
            final_step, mod_number, num_zones, ierr)
         use imex_work
         use imex_output
         integer, intent(in) :: max_steps_for_this_call
         real(dp), intent(out) :: age, timestep, total_energy
         logical, intent(out) :: final_step
         integer, intent(out) :: mod_number, num_zones, ierr
         integer :: step
         include 'formats'
         ierr = 0
         final_step = .true.
         do step = 1, max_steps_for_this_call
            model_number = model_number + 1
            call save_start_values ! T_start, prim_start, cons_start
            final_step = do_step( &
               T_start, prim_start, cons_start, & ! input
               T_1, prim_1, cons_1, grad_cell, flux_face, & ! work
               d_cons_dt_X0, d_ETOT_dt_I1, d_cons_dt_X1, d_ETOT_dt_I2, & ! work
               T, prim, cons, L, v_face, & ! output
               ierr)
            if (ierr /= 0 .or. final_step) exit
         end do
         age = time
         timestep = dt
         num_zones = nz
         mod_number = model_number
         total_energy = sum(cons(i_etot,1:nz))
         
         contains
         
         subroutine save_start_values
            integer :: k, j
            do k=1,nz
               T_start(k) = T(k)
               do j=1,nvar
                  prim_start(j,k) = prim(j,k)
                  cons_start(j,k) = cons(j,k)
               end do
            end do
         end subroutine save_start_values
         
      end subroutine steps_imex
      
      
      logical function do_step( & 
         ! doesn't assume or enforce exact consistency between T and prim
            T_0, prim_0, cons_0, & ! input
            T_1, prim_1, cons_1, grad_cell, flux_face, & ! work
            d_cons_dt_X0, d_ETOT_dt_I1, d_cons_dt_X1, d_ETOT_dt_I2, & ! work
            T_2, prim_2, cons_2, L, v_face, & ! output
            ierr) result(final_step)
         use imex_work, only: &
            time, dt, dt_next, newton_iter, model_number, stage, gam, delta
         real(dp), intent(in), dimension(:) :: T_0 ! input
         real(dp), intent(in), dimension(:,:) ::  cons_0, prim_0 ! input
         real(dp), intent(out), dimension(:) :: & ! work
            T_1, d_ETOT_dt_I1,  d_ETOT_dt_I2
         real(dp), intent(out), dimension(:,:) :: & ! work
            prim_1, cons_1, grad_cell, flux_face, &
            d_cons_dt_X0, d_cons_dt_X1
         real(dp), intent(out), dimension(:) :: T_2, L, v_face ! output
         real(dp), intent(out), dimension(:,:) :: prim_2, cons_2 ! output
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         newton_iter = 0
         dt = dt_next
         write(*,2) 'do_step dt', model_number, dt, time
         if (dt == 0d0) stop 'do_step'
         final_step = (time + dt >= time_end)
         if (final_step) dt = time_end - time
         stage = 1
         if (stage1_only) then
            gam = 1d0
            delta = 0d0
            call stage1( &
               T_0, prim_0, cons_0, & ! input
               grad_cell, flux_face, & ! work
               d_cons_dt_X0, d_ETOT_dt_I1, & ! output
               T_2, prim_2, cons_2, L, v_face, & ! output
               ierr)
            if (ierr /= 0) then
               write(*,2) 'stage1 failed', model_number
               return
            end if
         else 
            ! 2016_Wang_Shu_Zang; 1997_ascher_ruuth_spiteri
            gam = 1d0 - sqrt(2d0)/2d0
            delta = 1d0 - 1d0/(2d0*gam)
            call stage1( &
               T_0, prim_0, cons_0, & ! input
               grad_cell, flux_face, & ! work
               d_cons_dt_X0, d_ETOT_dt_I1, & ! output
               T_1, prim_1, cons_1, L, v_face, & ! output
               ierr)
            if (ierr /= 0) then
               write(*,2) 'stage1 failed', model_number
               return
            end if
            stage = 2
            call stage2( &
               cons_0, T_1, prim_1, cons_1, d_cons_dt_X0, d_ETOT_dt_I1, & ! input
               grad_cell, flux_face, d_cons_dt_X1, d_ETOT_dt_I2, & ! work
               T_2, prim_2, cons_2, L, v_face, & ! output
               ierr)
            if (ierr /= 0) then
               write(*,2) 'stage2 failed', model_number
               return
            end if
         end if
         time = time + dt
         stage = 0
         dt_next = pick_next_timestep(T_2, prim_2, T_0, prim_0)
         !write(*,1) 'dt_next', dt_next
      end function do_step
      
      
      subroutine stage1( &
            T_0, prim_0, cons_0, & ! input
            grad_cell, flux_face, & ! work
            d_cons_dt_X0, d_ETOT_dt_I1, & ! output
            T_1, prim_1, cons_1, L_1, v_face_1, & ! output
            ierr)
         use imex_work, only: dt, gam, dVol
         real(dp), intent(in), dimension(:) :: T_0
         real(dp), intent(in), dimension(:,:) ::  cons_0, prim_0
         real(dp), intent(out), dimension(:) :: &
            d_ETOT_dt_I1, T_1, L_1, v_face_1
         real(dp), intent(out), dimension(:,:) :: &
            grad_cell, flux_face, d_cons_dt_X0, prim_1, cons_1
         integer, intent(out) :: ierr
         integer :: k, j
         include 'formats'
         ierr = 0
         call calc_explicit( &
            T_0, prim_0, cons_0, grad_cell, flux_face, d_cons_dt_X0, v_face_1, ierr)
         if (ierr /= 0) return
         do k=1,nz ! update explicit intermediate result
            do j=1,nvar
               cons_1(j,k) = cons_0(j,k) + dt*gam*d_cons_dt_X0(j,k)
               prim_1(j,k) = cons_1(j,k)/dVol(k)
            end do
            !write(*,2) 'v rho mom etot', k, prim_1(i_mom,k)/prim_1(i_rho,k), prim_1(:,k)
         end do
         call calc_implicit(T_0, prim_1, cons_1, d_ETOT_dt_I1, T_1, L_1, ierr)
         if (ierr /= 0) return
         do k=1,nz ! update etot
            cons_1(i_etot,k) = cons_1(i_etot,k) + dt*gam*d_ETOT_dt_I1(k)
            prim_1(i_etot,k) = cons_1(i_etot,k)/dVol(k)
         end do
      end subroutine stage1
      
      
      subroutine stage2( &
            cons_0, T_1, prim_1, cons_1, d_cons_dt_X0, d_ETOT_dt_I1, & ! input
            grad_cell, flux_face, d_cons_dt_X1, d_ETOT_dt_I2, & ! work
            T_2, prim_2, cons_2, L_2, v_face_2, & ! output
            ierr)
         use imex_work, only: dt, gam, delta, dVol
         real(dp), intent(in), dimension(:) :: &
            T_1, d_ETOT_dt_I1
         real(dp), intent(in), dimension(:,:) :: &
            cons_0, prim_1, cons_1, d_cons_dt_X0
         real(dp), intent(out), dimension(:) :: &
            d_ETOT_dt_I2, T_2, L_2, v_face_2
         real(dp), intent(out), dimension(:,:) :: &
            grad_cell, flux_face, d_cons_dt_X1, prim_2, cons_2
         integer, intent(out) :: ierr
         integer :: k, j
         include 'formats'
         ierr = 0
         call calc_explicit( &
            T_1, prim_1, cons_1, grad_cell, flux_face, d_cons_dt_X1, v_face_2, ierr)
         if (ierr /= 0) return
         do k=1,nz ! update explicit intermediate result
            do j=1,nvar
               cons_2(j,k) = cons_0(j,k) + &
                  dt*(delta*d_cons_dt_X0(j,k) + (1d0-delta)*d_cons_dt_X1(j,k)) 
               if (j == i_etot) &
                  cons_2(j,k) = cons_2(j,k) + dt*(1d0-gam)*d_ETOT_dt_I1(k)
               prim_2(j,k) = cons_2(j,k)/dVol(k)
            end do
         end do
         call calc_implicit(T_1, prim_2, cons_2, d_ETOT_dt_I2, T_2, L_2, ierr)
         if (ierr /= 0) return
         do k=1,nz ! update etot
            cons_2(i_etot,k) = cons_2(i_etot,k) + dt*gam*d_ETOT_dt_I2(k)
            prim_2(i_etot,k) = cons_2(i_etot,k)/dVol(k)
         end do
      end subroutine stage2
            
      
      real(dp) function pick_next_timestep(T_2, prim_2, T_0, prim_0)
         use imex_work, only: dt
         real(dp), intent(in), dimension(:) :: T_2, T_0
         real(dp), intent(in), dimension(:,:) ::  prim_2, prim_0
         real(dp) :: dt_advection, dt_grid, dt_rel_dE, dt_front_v, dt_max_new
         include 'formats'
         dt_advection = get_min_dt_advection(T_2, prim_2)
         dt_grid = get_min_dt_grid()
         dt_rel_dE = get_dt_rel_dE(prim_0, prim_2)
         dt_front_v = get_dt_front_v(prim_0, prim_2)
         dt_max_new = dt*max_timestep_factor
         pick_next_timestep = &
            min(dt_advection, dt_grid, dt_rel_dE, dt_front_v, dt_max_new)
         !write(*,1) 'dt_advection', dt_advection
         !write(*,1) 'dt_grid', dt_grid
         !write(*,1) 'dt_rel_dE', dt_rel_dE
         !write(*,1) 'dt_front_v', dt_front_v
         !write(*,1) 'dt_max_new', dt_max_new
         !write(*,1) 'pick_next_timestep', pick_next_timestep
      end function pick_next_timestep
      
      
      real(dp) function get_dt_rel_dE(prim_start, prim_end) result(dt_rel_dE)
         use imex_work, only: dt
         real(dp), intent(in) :: prim_start(:,:), prim_end(:,:)
         real(dp) :: etot_min, max_dE_div_E
         etot_min = minval(prim_end(i_etot,1:nz))
         max_dE_div_E = &
            maxval(abs(prim_end(i_etot,1:nz) - prim_start(i_etot,1:nz))/ &
               (prim_end(i_etot,1:nz) + etot_min))
         if (max_dE_div_E <= 1d-5*limit_dE_div_E) then
            dt_rel_dE = dt*max_timestep_factor
         else
            dt_rel_dE = dt*sqrt(limit_dE_div_E/max_dE_div_E)
         end if
      end function get_dt_rel_dE
      
      
      real(dp) function get_dt_front_v(prim_start, prim_end) result(dt_front_v)
         use imex_work, only: dt
         real(dp), intent(in) :: prim_start(:,:), prim_end(:,:)
         real(dp) :: sum_dE_dt, sum_dE
         integer :: k
         sum_dE_dt = 0d0
         do k=1,nz
            sum_dE_dt = abs(prim_end(i_etot,k) - prim_start(i_etot,k))/dt
         end do
         sum_dE = 0d0
         do k=2,nz-1
            sum_dE = 0.5d0*abs(prim_end(i_etot,k-1) - prim_end(i_etot,k+1))
         end do
         if (sum_dE_dt <= 1d-5*sum_dE + 1d-20) then
            dt_front_v = dt*max_timestep_factor
         else
            dt_front_v = CFL_front*sum_dE/sum_dE_dt
         end if
      end function get_dt_front_v
      
      
      real(dp) function get_min_dt_advection(T, prim) result(dt_advection)
         use imex_work, only: dr
         use imex_output, only: csound
         real(dp), intent(in) :: T(:)
         real(dp), intent(in) :: prim(:,:)
         real(dp) :: dt_cell, vel, rho
         integer :: k
         dt_advection = 1d99
         do k=1,nz
            rho = prim(i_rho,k)
            csound(k) = get_csound(T(k), rho)
            vel = prim(i_mom,k)/rho
            dt_cell = dr(k)/abs(vel + csound(k))
            if (dt_cell < dt_advection) dt_advection = dt_cell
         end do
         dt_advection = CFL*dt_advection
      end function get_min_dt_advection
      
      
      real(dp) function get_min_dt_grid() result(dt_grid) ! cheng 2020   3.10
         use imex_work, only: dt, dr
         use imex_output, only: v_face
         real(dp) :: dt_cell
         integer :: k
         dt_grid = 1d99
         do k=2,nz
            dt_cell = min(dr(k),dr(k-1))/max(1d-20,abs(v_face(k)))
            if (dt_cell < dt_grid) dt_grid = dt_cell
         end do
         dt_grid = 0.45d0*dt_grid
      end function get_min_dt_grid
      
      
!!! implicit part      
      
      
      subroutine calc_implicit( &
            T_prev, prim, cons, & ! input
            d_ETOT_dt_I, T, L, & ! output
            ierr)
         use imex_work, only: newton_iter, stage, Eeos_x, rhs, dVol
         real(dp), intent(in) :: T_prev(:)
         real(dp), intent(in), dimension(:,:) :: prim, cons
         real(dp), intent(out), dimension(:) :: d_ETOT_dt_I, T, L ! output
         integer, intent(out) :: ierr       
         integer :: k
         real(dp) :: rho, v, etot, ekin, resid_norm, resid_max
         logical :: converged
         include 'formats'        
         ierr = 0 
         do k=1,nz
            rho = prim(i_rho,k)
            v = prim(i_mom,k)/rho
            etot = prim(i_etot,k)
            ekin = 0.5d0*rho*pow2(v)
            Eeos_X(k) = (etot - ekin)*dVol(k)
            T(k) = T_prev(k) ! initial guess
         end do         
         converged = .false.
         do newton_iter = 1, newton_iter_max         
            !$OMP PARALLEL DO PRIVATE(k)
            do k=1,nz
               call store_matrix_equation(k) ! sets d_ETOT_dt_I(k) and L(k)
            end do
            !$OMP END PARALLEL DO
            resid_norm = sum(abs(rhs(1:nz)))/nz
            resid_max = maxval(abs(rhs(1:nz)))
            !write(*,4) 'resid norm max', newton_iter, stage, model_number, resid_norm, resid_max
            if (resid_norm <= tol_resid_norm .and. resid_max <= tol_resid_max) then
               converged = .true.
               exit
            end if
            call solve_matrix_equation(ierr) ! updates T
         end do         
         newton_iter = 0
         if (.not. converged) then
            ierr = -1
            return
         end if
         
         contains
         
         subroutine store_matrix_equation(k)
            use imex_work, only: dt, gam, dVol, sub, diag, sup, rhs
            integer, intent(in) :: k
            type(auto_diff_real_4var_order1) :: &
               T_00, L_00, L_p1, Eeos_expected, Eeos_actual, residual
            include 'formats'
            T_00 = wrap_T_00(T,k)
            L_00 = get_L_face(T,prim,k)
            L(k) = L_00%val
            if (k == nz) then
               L_p1 = L_inner
            else
               L_p1 = shift_p1(get_L_face(T,prim,k+1))
            end if
            d_ETOT_dt_I(k) = L_p1%val - L_00%val
            Eeos_expected = Eeos_X(k) + dt*gam*(L_p1 - L_00)
            Eeos_actual = (crad*pow4(T_00) + rho*Cv*T_00)*dVol(k)
            residual = Eeos_actual - Eeos_expected
            rhs(k) = -residual%val
            if (k > 1) sub(k-1) = residual%d1val1
            diag(k) = residual%d1val2
            if (k < nz) sup(k) = residual%d1val3
         end subroutine store_matrix_equation
         
         subroutine solve_matrix_equation(ierr)
            use imex_work, only: newton_iter, stage, deltaT, model_number
            integer, intent(out) :: ierr
            integer :: k
            include 'formats'
            ierr = 0
            call solve_tridiag(deltaT, nz, ierr)
            if (ierr /= 0) then
               write(*,4) 'solve_tridiag failed', newton_iter, stage, model_number
               return
            end if
            do k=1,nz
               T(k) = T(k) + deltaT(k)
            end do            
         end subroutine solve_matrix_equation

      end subroutine calc_implicit
      
      
      function get_L_face(T,prim,k) result(L_face) ! luminosity (by radiative diffusion)
         use imex_work, only: area, dr_bar
         type(auto_diff_real_4var_order1) :: L_face
         real(dp), intent(in) :: T(:), prim(:,:)
         integer, intent(in) :: k
         type(auto_diff_real_4var_order1) :: T_m1, T_00, TC
         if (k > nz) then
            L_face = L_inner
            return
         end if
         T_00 = wrap_T_00(T,k)
         if (k == 1) then
            L_face = Lsurf_factor * area(k) * clight * crad * pow4(T_00)
         else
            T_m1 = wrap_T_m1(T,k)
            TC = get_TC_face(T_00,T_m1,k)
            L_face = -area(k)*TC*(T_m1 - T_00)/dr_bar(k)
         end if
         
         contains
         
         function get_TC_face(T_00,T_m1,k) result(TC_face) ! thermal conductivity at face
            use imex_work, only: dr
            integer, intent(in) :: k ! k > 1 and k <= nz
            type(auto_diff_real_4var_order1) :: T_00, T_m1, TC_face
            type(auto_diff_real_4var_order1) :: TC_00, TC_m1
            real(dp) :: rho_face
            TC_00 = TC0
            TC_m1 = TC0
            if (TCa /= 0d0) then
               TC_00 = TC_00*pow(prim(i_rho,k),TCa)
               TC_m1 = TC_m1*pow(prim(i_rho,k-1),TCa)
            end if
            if (TCb /= 0d0) then
               TC_00 = TC_00*pow(T_00,TCb)
               TC_m1 = TC_m1*pow(T_m1,TCb)
            end if
            TC_face = (dr(k-1)*TC_00 + dr(k)*TC_m1)/(dr(k-1) + dr(k))
         end function get_TC_face
         
      end function get_L_face


      subroutine solve_tridiag(x, n, ierr)
         !      sub - sub-diagonal
         !      diag - the main diagonal
         !      sup - sup-diagonal
         !      rhs - right hand side
         !      x - the answer
         !      n - number of equations
         use imex_work, only: sub, diag, sup, rhs, bp, vp, xp
         integer, intent(in) :: n
         real(dp), dimension(:), intent(out) :: x ! output
         integer, intent(out) :: ierr
         real(dp) :: m
         integer i
         ierr = 0
         sub(n) = 0d0
         sup(n) = 0d0
         bp(1) = diag(1)
         vp(1) = rhs(1)
         do i = 2,n
            m = sub(i-1)/bp(i-1)
            bp(i) = diag(i) - m*sup(i-1)
            vp(i) = rhs(i) - m*vp(i-1)
         end do
         xp(n) = vp(n)/bp(n)
         x(n) = xp(n)
         do i = n-1, 1, -1
            xp(i) = (vp(i) - sup(i)*xp(i+1))/bp(i)
            x(i) = xp(i)
         end do
      end subroutine solve_tridiag

      
!!! explicit part      
      
      
      subroutine calc_explicit( &
            T, prim, cons, & ! input
            grad_cell, flux_face, & ! work
            d_cons_dt_X, v_face, & ! output
            ierr)
         real(dp), intent(in) :: T(:), prim(:,:), cons(:,:)
         real(dp), intent(out) :: grad_cell(:,:), flux_face(:,:)
         real(dp), intent(out) :: d_cons_dt_X(:,:), v_face(:)
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         call get_grad_cell(prim, grad_cell, ierr)
         if (ierr /= 0) return
         call get_flux_face(prim, grad_cell, flux_face, v_face, ierr)
         if (ierr /= 0) return
         call get_d_cons_dt_X(flux_face, prim, T, d_cons_dt_X, ierr)
      end subroutine calc_explicit
      
      
      subroutine get_grad_cell( &
            prim, & ! input
            grad_cell, & ! output
            ierr)
         real(dp), intent(in) :: prim(:,:)
         real(dp), intent(out) :: grad_cell(:,:)
         integer, intent(out) :: ierr
         integer :: k, op_err
         ierr = 0
         !$OMP PARALLEL DO PRIVATE(k,op_err)
         do k=1,nz
            op_err = 0
            call get1_grad_cell(k, op_err)
            if (op_err /= 0) ierr = op_err
         end do
         !$OMP END PARALLEL DO
         
         contains
         
         subroutine get1_grad_cell(k, ierr)
            use imex_work, only: dr_bar
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            real(dp) :: a, b
            integer :: j
            include 'formats'
            ierr = 0
            if (k == 1 .or. k == nz) then
               grad_cell(1:nvar,k) = 0d0 
               return
            end if
            do j=1,nvar ! simple minmod
               a = (prim(j,k-1) - prim(j,k))/dr_bar(k)
               b = (prim(j,k) - prim(j,k+1))/dr_bar(k+1)
               if (a*b <= 0d0) then
                  grad_cell(j,k) = 0d0
               else if (abs(b) < abs(a)) then
                  grad_cell(j,k) = b
               else
                  grad_cell(j,k) = a
               end if
            end do            
         end subroutine get1_grad_cell
                  
      end subroutine get_grad_cell
   
      
      subroutine get_flux_face( &
            prim, grad_cell, & ! input
            flux_face, v_face, & ! output
            ierr)
         use imex_work, only: area
         real(dp), intent(in) :: prim(:,:), grad_cell(:,:)
         real(dp), intent(out) :: flux_face(:,:), v_face(:)
         integer, intent(out) :: ierr
         integer :: j, k, op_err
         include 'formats'
         ierr = 0
         
         !$OMP PARALLEL DO PRIVATE(k,op_err)
         do k=2,nz
            op_err = 0
            if (.false.) then
               call get1_flux_face_HLLC(k, op_err)
            else
               call get1_flux_face_simple(k, op_err)
            end if
            if (op_err /= 0) ierr = op_err
         end do
         !$OMP END PARALLEL DO
         
         ! "outflow" outer BC     
         do j=1,nvar
            flux_face(j,1) = flux_face(j,2)*area(2)/area(1)
         end do
         
         contains
         
         subroutine get1_flux_face_simple(k, ierr) ! Local Lax Friedrichs flux
            ! 2010a_kadioglu_knoll, eqn 9
            use imex_work, only: dr, P_face
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            integer :: j
            real(dp) :: rhoL, rhoR, momL, momR, vL, vR, etotL, etotR, &
               TL, TR, csL, csR, alpha, PgasL, PgasR, PradL, PradR, PL, PR
            real(dp), dimension(nvar) :: primL, primR, fluxL, fluxR
            include 'formats'
            ierr = 0
            do j=1,nvar
               primL(j) = prim(j,k) + grad_cell(j,k)*dr(k)/2d0 ! on left side of face k
               primR(j) = prim(j,k-1) - grad_cell(j,k-1)*dr(k-1)/2d0 ! on right side of face k
            end do
            rhoL = primL(i_rho)
            rhoR = primR(i_rho)
            momL = primL(i_mom)
            momR = primR(i_mom)
            vL = momL/rhoL
            vR = momR/rhoR
            v_face(k) = 0.5d0*(vL + vR)  
            etotL = primL(i_etot)
            etotR = primR(i_etot)
            TL = get_T(rhoL, momL, etotL)
            TR = get_T(rhoR, momR, etotR)
            csL = get_csound(TL,rhoL)
            csR = get_csound(TR,rhoR)
            PgasL = (gamma-1)*rhoL*Cv*TL
            PgasR = (gamma-1)*rhoR*Cv*TR
            PradL = crad*pow4(TL)/3d0
            PradR = crad*pow4(TR)/3d0
            PL = PgasL + PradL
            PR = PgasR + PradR
            P_face(k) = 0.5d0*(PL + PR)
            fluxL(i_rho) = rhoL*vL
            fluxR(i_rho) = rhoR*vR
            fluxL(i_mom) = momL*vL
            fluxR(i_mom) = momR*vR
            if (include_P_in_momentum_flux) then
               fluxL(i_mom) = momL*vL + PL
               fluxR(i_mom) = momR*vR + PR
            end if
            fluxL(i_etot) = (etotL + PL)*vL
            fluxR(i_etot) = (etotR + PR)*vR
            alpha = 0d0 ! min(abs(vL-csL),abs(vL+csL),abs(vR-csR),abs(vR+csR))
            do j=1,nvar
               flux_face(j,k) = 0.5d0*(fluxR(j)+fluxL(j) - alpha*(primR(j)-primL(j)))
            end do          
         end subroutine get1_flux_face_simple
         
         
         subroutine get1_flux_face_HLLC(k, ierr) 
            ! for P in momentum flux and geometry momentum source term
            use imex_work, only: dr
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            integer :: j
            real(dp) :: TL, TR, csL, csR, PL, PR, rhoL, rhoR, momL, momR, vL, vR, &
               etotL, etotR, sqrt_rhoL, sqrt_rhoR, rho_v_L, rho_v_R, &
               Sl, Sr, Ss, d_bar, eta_2, PgasL, PgasR, PradL, PradR, prim_star_prefactor
            real(dp), dimension(nvar) :: primL, primR, fluxL, fluxR, prim_star
            include 'formats'
            ierr = 0
            do j=1,nvar
               primL(j) = prim(j,k) + grad_cell(j,k)*dr(k)/2d0 ! on left side of face k
               primR(j) = prim(j,k-1) - grad_cell(j,k-1)*dr(k-1)/2d0 ! on right side of face k
            end do
            rhoL = primL(i_rho); sqrt_rhoL = sqrt(rhoL)
            rhoR = primR(i_rho); sqrt_rhoR = sqrt(rhoR)
            momL = primL(i_mom); momR = primR(i_mom)
            vL = momL/rhoL; vR = momR/rhoR
            etotL = primL(i_etot); etotR = primR(i_etot)
            TL = get_T(rhoL, momL, etotL); TR = get_T(rhoR, momR, etotR)
            PgasL = (gamma-1)*rhoL*Cv*TL; PgasR = (gamma-1)*rhoR*Cv*TR
            PradL = crad*pow4(TL)/3d0; PradR = crad*pow4(TR)/3d0
            PL = PgasL + PradL; PR = PgasR + PradR
            csL = get_csound(TL,rhoL); csR = get_csound(TR,rhoR)
            fluxL(i_rho) = rhoL*vL; fluxR(i_rho) = rhoR*vR
            fluxL(i_mom) = momL*vL + PL; fluxR(i_mom) = momR*vR + PR
            fluxL(i_etot) = (etotL + PL)*vL; fluxR(i_etot) = (etotR + PR)*vR
            v_face(k) = (sqrt_rhoL*vL + sqrt_rhoR*vR)/(sqrt_rhoL + sqrt_rhoR) 
            if (.true.) then ! debugging
               do j=1,nvar
                  flux_face(j,k) = 0.5d0*(fluxR(j) + fluxL(j))
               end do
               return
            else if (.false.) then ! Sl and Sr using 1988_einfeldt method
               eta_2 = 0.5*sqrt_rhoL*sqrt_rhoR/pow2(sqrt_rhoL + sqrt_rhoR)        ! Toro eqn 10.54
               d_bar = sqrt((sqrt_rhoL*pow2(csL) + sqrt_rhoR*pow2(csR))/ &
                            (sqrt_rhoL + sqrt_rhoR) + eta_2*pow2(vR - vL))       ! Toro eqn 10.53
               Sl = v_face(k) - d_bar
               Sr = v_face(k) + d_bar        ! Toro eqn 10.52    
            else ! simple wave speeds       ! Toro eqn 10.47
               Sl = min(vL - csL, vR - csR)
               Sr = max(vL + csL, vR + csR)
            end if  
            if (Sl >= 0d0) then
               do j=1,nvar
                  flux_face(j,k) = fluxL(j)
               end do
            else if (Sr <= 0d0) then
               do j=1,nvar
                  flux_face(j,k) = fluxR(j)
               end do
            else 
               rho_v_L = rhoL*(Sl - vL)
               rho_v_R = rhoR*(Sr - vR)
               Ss = (PR - PL + vL*rho_v_L - vR*rho_v_R)/(rho_v_L - rho_v_R)  !  Toro 10.37
               if (Ss > 0d0) then
                  call get_prim_star(Ss, Sl, vL, rhoL, etotL, PL, prim_star)
                  do j=1,nvar
                     flux_face(j,k) = fluxL(j) + Sl*(prim_star(j) - primL(j))      ! Toro 10.38
                  end do
               else 
                  call get_prim_star(Ss, Sr, vR, rhoR, etotR, PR, prim_star)
                  do j=1,nvar
                     flux_face(j,k) = fluxR(j) + Sr*(prim_star(j) - primR(j))      ! Toro 10.38
                  end do
               end if
            end if            
         end subroutine get1_flux_face_HLLC
      
         subroutine get_prim_star(Ss, Sk, vK, rhoK, etotK, PK, prim_star)     ! Toro 10.37
            real(dp), intent(in) :: Ss, Sk, vK, rhoK, etotK, PK
            real(dp), intent(out) :: prim_star(nvar)
            real(dp) :: prim_star_prefactor
            prim_star_prefactor = rhoK*(Sk - vK)/(Sk - Ss)
            prim_star(i_rho) = prim_star_prefactor
            prim_star(i_mom) = prim_star_prefactor*Ss
            prim_star(i_etot) = prim_star_prefactor* &
               (etotK/rhoK + (Ss - vK)*(Ss + PK/(rhoK*(Sk - vK))))
         end subroutine get_prim_star
         
      end subroutine get_flux_face
      
      
      real(dp) function get_csound(T, rho) ! 2020_cheng_shu_song, eqn 2.7
         real(dp), intent(in) :: T, rho
         real(dp) :: Pgas, Prad, z, Gamma1
         Pgas = (gamma-1d0)*Cv*rho*T
         Prad = crad*pow4(T)/3
         z = Prad/Pgas
         Gamma1 = &
            (gamma/(gamma-1d0) + z*(20d0 * 16d0*z))/((1d0/(gamma-1d0) + 12d0*z)*(1d0 + z))
         get_csound = sqrt(Gamma1*(Pgas + Prad)/rho)
      end function get_csound


      real(dp) function get_etot(rho, mom, T) result(etot)
         real(dp), intent(in) :: rho, mom, T
         real(dp) :: v, Egas, Erad, Ekin
         v = mom/rho
         Egas = rho*Cv*T
         Erad = crad*pow4(T)
         Ekin = 0.5d0*mom*v
         etot = Egas + Erad + Ekin
      end function get_etot
      
      
      real(dp) function get_T(rho, mom, etot) result(T) 
         ! solve for T: etot == rho*Cv*T + 0.5*mom*v + crad*T^4   2020_cheng_shu_song
         real(dp), intent(in) :: rho, mom, etot
         real(dp) :: rhoCv, v, c1, c2, s3, s4, s1, s2, s, sqrt_2s, b1, b1_3, b2, b3, b4, b5
         rhoCv = rho*Cv
         v = mom/rho
         if (crad > 1d-6) then
            c1 = rhoCv/crad
            c2 = -1d0/crad*(etot - 0.5d0*mom*v)
            s3 = sqrt(pow4(c1)/256d0 - pow3(c2)/27d0)
            s4 = pow2(c1)/16d0
            s1 = s3 + s4
            s2 = s3 - s4
            s = pow(s1,1d0/3d0) - pow(s2,1d0/3d0)  ! Cheng 2.13
            sqrt_2s = sqrt(2d0*s)
            T = 0.5d0*(-sqrt_2s + sqrt(-2d0*s + 2d0*c1/sqrt_2s))  ! Cheng 2.12
         else ! cheng 2.14
            b1 = (etot - 0.5d0*mom*v)/rhoCv
            b1_3 = pow3(b1)
            b2 = -b1*b1_3/rhoCv
            b3 = -2d0*b1_3*b2/rhoCv
            b4 = -(6d0*pow2(b1*b2) + 2d0*b1_3*b2)/rhoCv
            b5 = (2*b1_3*b2 + 2d0*pow2(b1)*b2*b3)/rhoCv
            T = b1 + crad*(b2 + crad*(b3 + crad*(b4 + crad*b5)))
         end if
         if (is_bad(T) .or. T <= 0d0) T = (etot - 0.5d0*mom*v)/rhoCv ! cheng 2.16
      end function get_T
      
      
      subroutine get_d_cons_dt_X( &
            flux_face, prim, T, & ! input
            d_cons_dt_X, & ! output
            ierr)
         real(dp), intent(in) :: flux_face(:,:), prim(:,:), T(:)
         real(dp), intent(out) :: d_cons_dt_X(:,:)
         integer, intent(out) :: ierr
         integer :: k, op_err
         ierr = 0
!         !$OMP PARALLEL DO PRIVATE(k,op_err)
         do k=1,nz
            op_err = 0
            call get1_d_cons_dt_X(k, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!         !$OMP END PARALLEL DO
         
         contains
         
         subroutine get1_d_cons_dt_X(k, ierr)
            use imex_work, only: area, dVol, dr, P_face
            use imex_output, only: L
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            integer :: j
            real(dp) :: Pgas, Prad, Ptot
            include 'formats'
            ierr = 0
            Pgas = (gamma-1)*prim(i_rho,k)*Cv*T(k)
            Prad = crad*pow4(T(k))/3d0
            Ptot = Pgas + Prad
            if (k == nz) then
               do j=1,nvar
                  d_cons_dt_X(j,k) = -area(k)*flux_face(j,k)
               end do
               if (include_P_in_momentum_flux) then
                  d_cons_dt_X(i_mom,k) = d_cons_dt_X(i_mom,k) + area(k)*Ptot
               else ! -dP/dr as source term
               end if
            else if (k == 1) then
               do j=1,nvar
                  d_cons_dt_X(j,k) = area(k+1)*flux_face(j,k+1)
               end do
               if (include_P_in_momentum_flux) then
                  d_cons_dt_X(i_mom,k) = d_cons_dt_X(i_mom,k) - area(k+1)*Ptot
               else ! -dP/dr as source term
               end if
            else
               do j=1,nvar
                  d_cons_dt_X(j,k) = area(k+1)*flux_face(j,k+1) - area(k)*flux_face(j,k)
               end do
               if (include_P_in_momentum_flux) then
                  d_cons_dt_X(i_mom,k) = d_cons_dt_X(i_mom,k) + &
                     (area(k) - area(k+1))*Ptot ! geometry source term
               else ! -dP/dr as source term
                  d_cons_dt_X(i_mom,k) = d_cons_dt_X(i_mom,k) + &
                     dVol(k)*(P_face(k+1) - P_face(k))/dr(k)
               end if
            end if            
         end subroutine get1_d_cons_dt_X
         
      end subroutine get_d_cons_dt_X
      
      
!!! auto_diff


      function wrap_T_m1(T, k) result(T_m1)
         real(dp), intent(in) :: T(:)
         type(auto_diff_real_4var_order1) :: T_m1
         integer, intent(in) :: k
         T_m1 = 0d0 
         if (k > 1) then
            T_m1 % val = T(k-1)
            T_m1 % d1val1 = 1d0
         end if
      end function wrap_T_m1


      function wrap_T_00(T, k) result(T_00)
         real(dp), intent(in) :: T(:)
         type(auto_diff_real_4var_order1) :: T_00
         integer, intent(in) :: k
         T_00 = 0d0 
         T_00 % val = T(k)
         T_00 % d1val2 = 1d0
      end function wrap_T_00


      function wrap_T_p1(T, k) result(T_p1)
         real(dp), intent(in) :: T(:)
         type(auto_diff_real_4var_order1) :: T_p1
         integer, intent(in) :: k
         T_p1 = 0d0 
         if (k < nz) then
            T_p1 % val = T(k+1)
            T_p1 % d1val3 = 1d0
         end if
      end function wrap_T_p1

      
      type(auto_diff_real_4var_order1) function shift_p1(val_00) result(val_p1)
         type(auto_diff_real_4var_order1), intent(in) :: val_00
         val_p1%val = val_00%val
         val_p1%d1val1 = 0d0
         val_p1%d1val2 = val_00%d1val1
         val_p1%d1val3 = val_00%d1val2
      end function shift_p1


      type(auto_diff_real_4var_order1) function shift_m1(val_00) result(val_m1)
         type(auto_diff_real_4var_order1), intent(in) :: val_00
         val_m1%val = val_00%val
         val_m1%d1val1 = val_00%d1val2
         val_m1%d1val2 = val_00%d1val3
         val_m1%d1val3 = 0d0
      end function shift_m1
      
      
!!! problem setup

      
      subroutine initialize_static_problem()
         integer :: k
         real(dp) :: R_max, deltaR, r00, rp1, T_init
         include 'formats'         
         allocate(r_init(nz), rho_init(nz), mom_init(nz), etot_init(nz))         
         ! grid points equally spaced in r
         R_max = 1d0
         deltaR = (R_max - R_inner)/nz
         rp1 = R_inner
         T_init = 1d0
         do k=nz,1,-1
            r00 = rp1 + deltaR
            r_init(k) = r00
            rho_init(k) = 1d0
            mom_init(k) = 0d0
            etot_init(k) = get_etot(rho_init(k), mom_init(k), T_init)
            rp1 = r00
         end do         
      end subroutine initialize_static_problem








      
      subroutine initialize_smooth_problem() ! 2001 bates et al "Smooth problem test"
         integer :: k
         real(dp) :: R_max, deltaR, E0, c0, rp1, r00, rmid, dV, E00, Ep1
         include 'formats'                  
         allocate(r_init(nz), rho_init(nz), mom_init(nz), etot_init(nz))         
         ! grid points equally spaced in r
         R_max = 1d0
         deltaR = (R_max - R_inner)/nz
         E0 = 100
         c0 = 0.25d0
         rp1 = R_inner
         do k=nz,1,-1
            r00 = rp1 + deltaR
            r_init(k) = r00
            rmid = 0.5d0*(r00 + rp1)
            rho_init(k) = 1d0/rmid
            dV = 4d0*pi/3d0*(pow3(r00) - pow3(rp1))
            E00 = E0*exp(-pow2(r00/c0))/pow3(c0*sqrt(pi))
            Ep1 = E0*exp(-pow2(rp1/c0))/pow3(c0*sqrt(pi))
            etot_init(k) = & ! Bates eqn 40
               (E0*(erf(r00/c0) - erf(rp1/c0)) - 2d0*pi*pow2(c0)*(r00*E00 - rp1*Ep1))/dV
            mom_init(k) = 0d0
            rp1 = r00
         end do         
      end subroutine initialize_smooth_problem








      
      subroutine initialize_Barenblatt_problem() ! 2001 bates et al "Barenblatt test"
         ! Nonlinear thermal conduction from a point source
         ! hydrodynamic motion disabled
         integer :: k
         real(dp) :: R_max, deltaR, rp1, r00, E0, dV1
         include 'formats'
         allocate(r_init(nz), rho_init(nz), mom_init(nz), etot_init(nz))
         ! grid points equally spaced in r
         R_max = 1d0
         deltaR = (R_max - R_inner)/nz
         rp1 = R_inner
         do k=nz,1,-1
            r00 = rp1 + deltaR
            r_init(k) = r00
            rho_init(k) = 1d0
            etot_init(k) = 0d0
            mom_init(k) = 0d0
            rp1 = r00
         end do
         E0 = 10d0
         dV1 = 4d0*pi/3d0*(pow3(r_init(1)) - pow3(R_inner))
         etot_init(nz) = E0/dV1      
      end subroutine initialize_Barenblatt_problem


   end module imex
   