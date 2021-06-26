 
   module imex_input
      use const_def, only: dp
      include 'controls_imex_def.inc' ! declarations only
      real(dp), dimension(:), allocatable :: &
         r_init, rho_init, mom_init, etot_init ! (nz)         
   end module imex_input 
 
 
   module imex_work
      use const_def, only: dp
      integer :: model_number, newton_iter, gmres_mr, &
         total_newton_iters, total_gmres_matvecs
      real(dp) :: resid_norm, tol_resid_norm, time, dt, dt_next, &
         Cv, dt_advection, dt_grid, dt_rel_dE, &
         dt_max_new, dt_front_v, total_KE
      real(dp), dimension(:), allocatable :: &
         sub, diag, sup, bp, vp, xp, &
         r, area, Vol, dr, dr_bar, dVol, & 
         rhs, deltaT, P_face, L_start, &
         T_start, d_ETOT_dt_implicit, &
         imex1, imex2, imex3, imex4, imex5, imex6, &
         gmr_r, gmr_cs, gmr_g, gmr_sn, gmr_y, gmr_p ! gmres
      real(dp), dimension(:,:), allocatable :: &
         prim_start, grad_cell, flux_face, &
         gmr_v, gmr_h ! gmres
   end module imex_work


   module imex_plot_data ! only set by call set_imex_plot_data
      use const_def, only: dp
      real(dp), dimension(:), allocatable :: &
         plot_r, plot_L, plot_rho, plot_T, plot_v, plot_csound, &
         plot_etot, plot_ekin, plot_eeos, plot_Peos
   end module imex_plot_data


   module imex_output
      use const_def, only: dp
      real(dp), dimension(:,:), allocatable :: prim
      real(dp), dimension(:), allocatable :: T, L, v_face, csound
   end module imex_output
 
 
   module imex
      use imex_input

      ! bates, knoll, rider, lowrie, mousseau
      ! On Consistent Time-Integration Methods for Radiation Hydrodynamics in
      ! the Equilibrium Diffusion Limit: Low-Energy-Density Regime
      ! Journal of Computational Physics 167, 99–130 (2001)
      
      ! kadioglu, knoll
      ! A fully second order implicit/explicit time integration technique
      ! for hydrodynamics plus nonlinear heat conduction problems
      ! Journal of Computational Physics 229 (2010) 3237–3249

      use const_def, only: dp, qp, clight, pi
      use math_lib
      use auto_diff
      use utils_lib, only: is_bad
      
      implicit none

      integer, parameter :: &
         i_rho = 1, i_mom = 2, i_etot = 3, nvar = i_etot         
         
         
      contains
   
      
      subroutine start_imex(total_energy_initial, ierr)
         use imex_work
         use imex_output
         real(dp), intent(out) :: total_energy_initial
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         write(*,*) 'start_imex'
         if (crad < 0d0) crad = 7.5657332502799993d-15 ! crad in mesa
         if (cgas < 0d0) cgas = 8.314462618d7 ! cgas in mesa
         Cv = cgas/(gamma - 1d0)
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
         model_number = initial_model_number
         time = initial_time
         dt_next = initial_dt
         total_newton_iters = 0
         total_gmres_matvecs = 0
         call alloc_work_arrays()
         call alloc_imex_plot_data_arrays()
         call set_grid_vars()
         call set_init_prim_and_T()
         total_energy_initial = dot_product(prim(i_etot,1:nz),dVol(1:nz))
                  
         contains
         
         subroutine set_init_prim_and_T()
            integer :: k, j
            real(dp) :: rho, mom, etot, &
               sum_EKIN_actual, sum_EKIN_for_ETOT, sum_EGAS, sum_ERAD, sum_ETOT
            include 'formats'
            do k=1,nz
               rho = rho_init(k)
               mom = mom_init(k)
               etot = etot_init(k)
               prim(i_rho,k) = rho
               prim(i_mom,k) = mom
               prim(i_etot,k) = etot
               T(k) = get_T(rho, mom, etot)
               L(k) = 0d0
               v_face(k) = 0d0
               imex1(k) = 0d0
               imex2(k) = 0d0
               imex3(k) = 0d0
               imex4(k) = 0d0
            end do
            call get_imex_total_energies(prim, T, &
               sum_EKIN_actual, sum_EKIN_for_ETOT, sum_EGAS, sum_ERAD, sum_ETOT)
         end subroutine set_init_prim_and_T
         
         subroutine alloc_work_arrays()
            integer :: neq, mr
            gmres_mr = 20
            mr = gmres_mr
            neq = nz ! only 1 var per zone for gmres
            allocate(&
               r(nz), area(nz), Vol(nz), dr(nz), dr_bar(nz), dVol(nz), & 
               rhs(nz), sub(nz), diag(nz), sup(nz), deltaT(nz), &
               bp(nz), vp(nz), xp(nz), P_face(nz), L_start(nz), &
               T_start(nz), T(nz), L(nz), v_face(nz), csound(nz), &
               d_ETOT_dt_implicit(nz), prim_start(nvar,nz), prim(nvar,nz), &
               imex1(nz), imex2(nz), imex3(nz), imex4(nz), imex5(nz), imex6(nz), &
               grad_cell(nvar,nz), flux_face(nvar,nz), &
               gmr_r(neq), gmr_cs(mr), gmr_sn(mr), gmr_g(mr+1), gmr_p(neq), &
               gmr_y(mr+1), gmr_v(neq,mr+1), gmr_h(mr+1,mr))
         end subroutine alloc_work_arrays

         subroutine alloc_imex_plot_data_arrays()
            use imex_plot_data
            allocate(&
               plot_r(nz), plot_L(nz), plot_rho(nz), plot_T(nz), plot_v(nz), &
               plot_csound(nz), plot_etot(nz), plot_ekin(nz), plot_eeos(nz), plot_Peos(nz))
         end subroutine alloc_imex_plot_data_arrays
      
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
      
      
      subroutine set_imex_plot_data()
         use imex_work
         use imex_output
         use imex_plot_data
         integer :: k
         real(dp) :: rho, Pgas, Prad, mom, v
         include 'formats'         
         do k=1,nz
            plot_r(k) = r(k)
            plot_L(k) = L(k)
            rho = prim(i_rho,k)
            plot_rho(k) = rho
            plot_T(k) = T(k)
            mom = prim(i_mom,k)
            v = mom/rho
            plot_v(k) =  v
            plot_csound(k) = csound(k)
            plot_etot(k) = prim(i_etot,k)
            plot_ekin(k) = 0.5d0*rho*pow2(v)
            plot_eeos(k) = plot_etot(k) - plot_ekin(k)
            Pgas = get_Pgas(rho,T(k))
            Prad = get_Prad(T(k))
            plot_Peos(k) = Pgas + Prad
         end do
      end subroutine set_imex_plot_data
      
      
      subroutine get_imex_total_energies( &
            prim, T, sum_EKIN_actual, sum_EKIN_for_ETOT, sum_EGAS, sum_ERAD, sum_ETOT)
         use imex_work
         real(dp), intent(in) :: prim(:,:), T(:)
         real(dp), intent(out) :: &
            sum_EKIN_actual, sum_EKIN_for_ETOT, sum_EGAS, sum_ERAD, sum_ETOT
         integer :: k
         real(dp) :: rho, Etot, Egas, Erad, Ekin, v
         sum_ETOT = 0d0
         sum_EKIN_for_ETOT = 0d0
         sum_EGAS = 0d0
         sum_ERAD = 0d0
         sum_EKIN_actual = 0d0
         do k=1,nz
            rho = prim(i_rho,k)
            Etot = prim(i_etot,k)
            Egas = rho*Cv*T(k)
            Erad = crad*pow4(T(k))
            Ekin = Etot - (Egas + Erad)
            v = prim(i_mom,k)/rho
            sum_EKIN_actual = sum_EKIN_actual + 0.5d0*rho*pow2(v)*dVol(k)
            sum_EKIN_for_ETOT = sum_EKIN_for_ETOT + Ekin*dVol(k)
            sum_EGAS = sum_EGAS + Egas*dVol(k)
            sum_ERAD = sum_ERAD + Erad*dVol(k)
            sum_ETOT = sum_ETOT + Etot*dVol(k)
         end do      
      end subroutine get_imex_total_energies

      
      subroutine steps_imex( &
            max_steps_for_this_call, age, timestep, &
            final_step, mod_number, num_zones, ierr)
         use imex_work
         use imex_output
         integer, intent(in) :: max_steps_for_this_call
         real(dp), intent(out) :: age, timestep
         logical, intent(out) :: final_step
         integer, intent(out) :: mod_number, num_zones, ierr
         integer :: step, k
         include 'formats'
         ierr = 0
         final_step = .true.
         do step = 1, max_steps_for_this_call
            model_number = model_number + 1
            call save_start_values ! T_start, prim_start
            final_step = do_step( &
               T_start, prim_start, & ! input
               grad_cell, flux_face, d_ETOT_dt_implicit, & ! work
               T, prim, L, v_face, & ! output
               ierr)
            if (ierr /= 0 .or. final_step) exit
         end do
         age = time
         timestep = dt
         num_zones = nz
         mod_number = model_number
         call report_energies()
         if (final_step) then
            select case(problem_number)
            case (0)
               call finish_static_problem()
            case (1)
               call finish_Barenblatt_problem()
            case (2)
               call finish_smooth_problem()
            end select
         end if
         
         contains
         
         subroutine save_start_values
            integer :: k, j
            do k=1,nz
               T_start(k) = T(k)
               do j=1,nvar
                  prim_start(j,k) = prim(j,k)
               end do
            end do
         end subroutine save_start_values
         
         subroutine report_energies()
            real(dp) :: &
               sum_EKIN_actual, sum_EKIN_for_ETOT, sum_EGAS, sum_ERAD, sum_ETOT
            include 'formats'
            call get_imex_total_energies(prim, T, &
               sum_EKIN_actual, sum_EKIN_for_ETOT, sum_EGAS, sum_ERAD, sum_ETOT)
         end subroutine report_energies
         
      end subroutine steps_imex
      
      
      logical function do_step( & 
            T_0, prim_0, & ! input
            grad_cell, flux_face, d_ETOT_dt_implicit, & ! work
            T, prim, L, v_face, & ! output
            ierr) result(final_step)
         use imex_work, only: &
            time, dt, dt_next, newton_iter, model_number
         real(dp), intent(in), dimension(:) :: T_0 ! input
         real(dp), intent(in), dimension(:,:) ::  prim_0 ! input
         real(dp), intent(out), dimension(:) :: d_ETOT_dt_implicit
         real(dp), intent(out), dimension(:,:) :: grad_cell, flux_face
         real(dp), intent(out), dimension(:) :: T, L, v_face ! output
         real(dp), intent(out), dimension(:,:) :: prim ! output
         integer, intent(out) :: ierr
         integer :: k, j
         include 'formats'
         ierr = 0
         newton_iter = 0
         dt = dt_next
         if (time_centering) call save_start_values()
         !write(*,2) 'do_step dt time', model_number, dt, time
         if (dt == 0d0) then
            write(*,2) 'dt == 0d0', model_number, dt, time
            stop 'do_step'
         end if
         final_step = (time + dt >= time_end)
         if (final_step) dt = time_end - time
         ! bates 2001 setup with explicit followed by implicit
         if (thermal_diffusion_only) then
            ! no explicit
         else
            call calc_explicit( &
               T_0, prim_0, grad_cell, flux_face, prim, v_face, ierr)
            if (ierr /= 0) then
               write(*,2) 'do_step calc_explicit failed', model_number
               return
            end if
         end if
         call calc_implicit( &
            T_0, prim_0, grad_cell, flux_face, d_ETOT_dt_implicit, &
            prim, v_face, T, L, ierr)
         if (ierr /= 0) then
            write(*,2) 'do_step calc_implicit failed', model_number
            return
         end if
         time = time + dt
         dt_next = pick_next_timestep(T, prim, T_0, prim_0)
         !write(*,1) 'dt_next', dt_next
      end function do_step
      
      
      subroutine save_start_values()
         use imex_output, only: L
         use imex_work, only: L_start
         integer :: k
         do k=1,nz
            L_start(k) = L(k)
         end do
      end subroutine save_start_values
            
      
      real(dp) function pick_next_timestep(T, prim, T_0, prim_0)
         use imex_work, only: dt, &
            dt_advection, dt_grid, dt_rel_dE, dt_max_new, dt_front_v
         real(dp), intent(in), dimension(:) :: T, T_0
         real(dp), intent(in), dimension(:,:) ::  prim, prim_0
         include 'formats'
         dt_advection = get_min_dt_advection(T, prim)
         dt_grid = get_min_dt_grid()
         dt_rel_dE = get_dt_rel_dE(prim_0, prim)
         dt_front_v = get_dt_front_v(prim_0, prim)
         dt_max_new = dt*max_timestep_factor
         pick_next_timestep = &
            min(dt_advection, dt_rel_dE, dt_grid, dt_max_new, dt_front_v)
      end function pick_next_timestep
      
      
      real(dp) function get_dt_rel_dE(prim_start, prim_end) result(dt_rel_dE)
         ! bates 2001, eqn 27 (with bug fixed)
         use imex_work, only: dt, model_number
         real(dp), intent(in) :: prim_start(:,:), prim_end(:,:)
         real(dp) :: etot_end, etot_start, dE_div_E, max_dE_div_E
         integer :: k
         include 'formats'
         max_dE_div_E = 0
         do k=1,nz
            etot_end = prim_end(i_etot,k)
            etot_start = prim_start(i_etot,k)
            dE_div_E = abs(etot_end - etot_start)/max(etot_end,1d-4)
            if (dE_div_E > max_dE_div_E) max_dE_div_E = dE_div_E
         end do
         if (max_dE_div_E <= 1d-20*limit_dE_div_E + 1d-20) then
            dt_rel_dE = dt*max_timestep_factor
         else
            dt_rel_dE = dt*sqrt(limit_dE_div_E/max_dE_div_E)
         end if
      end function get_dt_rel_dE
      
      
      real(dp) function get_dt_front_v(prim_start, prim_end) result(dt_front_v)
         use imex_work, only: dt, model_number
         real(dp), intent(in) :: prim_start(:,:), prim_end(:,:)
         real(dp) :: sum_dE_dt, sum_dE
         integer :: k
         include 'formats'
         sum_dE_dt = 0d0
         do k=1,nz
            sum_dE_dt = abs(prim_end(i_etot,k) - prim_start(i_etot,k))/dt
         end do
         sum_dE = 0d0
         do k=2,nz-1
            sum_dE = 0.5d0*abs(prim_end(i_etot,k-1) - prim_end(i_etot,k+1))
         end do
         sum_dE = abs(prim_end(i_etot,1) - prim_end(i_etot,2))
         sum_dE = abs(prim_end(i_etot,nz-1) - prim_end(i_etot,nz))
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
            csound(k) = get_csound(rho, T(k))
            vel = prim(i_mom,k)/rho
            dt_cell = dr(k)/(abs(vel) + csound(k))
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
      
      
      subroutine calc_implicit( &
            T_prev, prim_0, grad_cell, flux_face, d_ETOT_dt_implicit, &
            prim, v_face, T, L, ierr)
         use imex_work, only: &
            total_newton_iters, dt, newton_iter, rhs, dVol, model_number
         real(dp), intent(in) :: T_prev(:), prim_0(:,:)
         real(dp), intent(out), dimension(:,:) :: grad_cell, flux_face, prim
         real(dp), intent(out), dimension(:) :: d_ETOT_dt_implicit, v_face, T, L
         integer, intent(out) :: ierr       
         integer :: k, j
         real(dp) :: rho, v, etot, ekin, resid_norm, resid_max
         logical :: converged
         include 'formats'        
         ierr = 0 
         do k=1,nz
            T(k) = T_prev(k) ! initial guess
         end do         
         if (thermal_diffusion_only) then
            do k=1,nz
               do j=1,nvar
                  prim(j,k) = prim_0(j,k)
               end do
            end do
         ! else prim set by explicit
         end if
         converged = .false.
         do newton_iter = 1, newton_iter_max
            total_newton_iters = total_newton_iters + 1
            if (.not. thermal_diffusion_only) then
               call calc_explicit( & ! redo with modified T
                  T, prim_0, grad_cell, flux_face, prim, v_face, ierr)
               if (ierr /= 0) then
                  write(*,2) 'stage1 calc_explicit failed', model_number
                  return
               end if
            end if         
            call update_T(T, T_prev, prim, d_ETOT_dt_implicit, L, ierr)
            if (ierr /= 0) then
               write(*,2) 'stage1 update_T failed', model_number
               return
            end if            
            resid_norm = sum(abs(rhs(1:nz)))/nz
            resid_max = maxval(abs(rhs(1:nz)))
            if (is_bad(resid_norm) .or. is_bad(resid_max)) then
               write(*,3) 'resid norm max', newton_iter, model_number, resid_norm, resid_max
               do k=1,nz
                  if (is_bad(rhs(k))) then
                     write(*,2) 'rhs', k, rhs(k)
                  end if
               end do
            end if
            if (resid_norm <= tol_resid_norm .and. resid_max <= tol_resid_max) then
               converged = .true.
               exit
            end if
         end do         
         newton_iter = 0
         if (.not. converged) then
            write(*,2) 'calc_implicit failed to converge', model_number, resid_norm, resid_max
            ierr = -1
            return
         end if
         do k=1,nz ! update etot
            prim(i_etot,k) = prim(i_etot,k) + dt*d_ETOT_dt_implicit(k)/dVol(k)
         end do

      end subroutine calc_implicit
      
      
      subroutine update_T(T, T_prev, prim, d_ETOT_dt_implicit, L, ierr)
         use imex_work, only: dt, newton_iter, rhs, dVol, deltaT, model_number
         real(dp), intent(inout) :: T(:)
         real(dp), intent(in) :: T_prev(:), prim(:,:)
         real(dp), intent(out), dimension(:) :: d_ETOT_dt_implicit, L
         integer, intent(out) :: ierr       
         integer :: k
         include 'formats'        
         ierr = 0 
         
         !$OMP PARALLEL DO PRIVATE(k)
         do k=1,nz
            call store_matrix_equation(T,T_prev,k) ! sets d_ETOT_dt_implicit(k) and L(k)
         end do
         !$OMP END PARALLEL DO    
         
         if (use_JFNK) then
            deltaT(1:nz) = 0d0
            call solve_matrix_equation_with_GMRES(T, deltaT, ierr)
         else
            call solve_matrix_equation(ierr)
         end if
         
         contains
         
         subroutine store_matrix_equation(T,T_prev,k)
            use imex_work, only: dt, Cv, dVol, sub, diag, sup, rhs
            integer, intent(in) :: k
            real(dp), intent(in) :: T(:), T_prev(:)
            type(auto_diff_real_4var_order1) :: &
               T_00, L_00, L_p1, dL, ETOT_expected, ETOT_actual, residual
            real(dp) :: rho, v
            include 'formats'
            T_00 = wrap_T_00(T,k)
            if (use_T_prev_for_L) then
               L_00 = get_L_face(T_prev,prim,k)
            else
               L_00 = get_L_face(T,prim,k)
            end if
            L(k) = L_00%val
            if (k == nz) then
               L_p1 = L_inner
            else
               if (use_T_prev_for_L) then
                  L_p1 = get_L_face(T_prev,prim,k+1)
               else
                  L_p1 = get_L_face(T,prim,k+1)
               end if
               L_p1 = shift_p1(L_p1)
            end if
            dL = L_00 - L_p1
            if (use_T_prev_for_L) then
               dL%d1val1 = 0d0
               dL%d1val2 = 0d0
               dL%d1val3 = 0d0
            end if
            d_ETOT_dt_implicit(k) = -dL%val
            ETOT_expected = prim(i_etot,k)*dVol(k) - dt*dL
            if (is_bad(ETOT_expected%val)) then
               !$omp critical (store_matrix_equation_crit1)
               write(*,2) 'ETOT_expected%val', k, ETOT_expected%val
               write(*,2) 'etot(k)', k, prim(i_etot,k)
               write(*,2) 'L_00%val', k, L_00%val
               write(*,2) 'L_p1%val', k, L_p1%val
               write(*,2) 'dt', k, dt
               write(*,*) 'use_T_prev_for_L', use_T_prev_for_L
               write(*,2) 'T_prev', k, T_prev(k)
               if (k < nz) write(*,2) 'T_prev', k+1, T_prev(k+1)
               stop 'store_matrix_equation'
               !$omp end critical (store_matrix_equation_crit1)
            end if
            rho = prim(i_rho,k)
            v = prim(i_mom,k)/rho
            ETOT_actual = (rho*Cv*T_00 + 0.5d0*rho*pow2(v))*dVol(k)
            if (.not. low_energy_density_regime) &
               ETOT_actual = ETOT_actual + crad*pow4(T_00)*dVol(k)
            if (is_bad(ETOT_actual%val)) then
               !$omp critical (store_matrix_equation_crit2)
               write(*,2) 'ETOT_actual%val', k, ETOT_actual%val
               write(*,2) 'T_00%val', k, T_00%val
               write(*,2) 'L_p1%val', k, L_p1%val
               write(*,2) 'dVol(k)', k, dVol(k)
               stop 'store_matrix_equation'
               !$omp end critical (store_matrix_equation_crit2)
            end if
            residual = ETOT_actual - ETOT_expected
            rhs(k) = -residual%val
            if (k > 1) sub(k-1) = residual%d1val1
            diag(k) = residual%d1val2
            if (k < nz) sup(k) = residual%d1val3
         end subroutine store_matrix_equation
         
         subroutine solve_matrix_equation(ierr)
            use imex_work, only: newton_iter, deltaT, model_number
            integer, intent(out) :: ierr
            integer :: k
            include 'formats'
            ierr = 0
            call solve_tridiag(deltaT, nz, ierr)
            if (ierr /= 0) then
               write(*,3) 'solve_tridiag failed', newton_iter, model_number
               return
            end if
            do k=1,nz
               T(k) = T(k) + deltaT(k)
            end do            
            !write(*,2) 'T', nz, T(nz)
         end subroutine solve_matrix_equation

      end subroutine update_T
      
      
      function get_L_face(T,prim,k) result(L_face) ! luminosity (by radiative diffusion)
         use imex_work, only: area, dr_bar, L_start
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
            if (low_energy_density_regime) then
               L_face = 0d0
            else
               L_face = Lsurf_factor * area(k) * clight * crad * pow4(T_00)
            end if
         else
            T_m1 = wrap_T_m1(T,k)
            TC = get_TC_face(T_00,T_m1,k)
            L_face = -area(k)*TC*(T_m1 - T_00)/dr_bar(k)
         end if
         if (time_centering) L_face = 0.5d0*(L_face + L_start(k))
         
         contains
         
         function get_TC_face(T_00,T_m1,k) result(TC_face) ! thermal conductivity at face
            use imex_work, only: dr
            integer, intent(in) :: k ! k > 1 and k <= nz
            type(auto_diff_real_4var_order1) :: T_00, T_m1, TC_face, T_avg
            type(auto_diff_real_4var_order1) :: TC_00, TC_m1
            real(dp) :: rho_avg
            rho_avg = 0.5d0*(prim(i_rho,k) + prim(i_rho,k-1))
            T_avg = 0.5d0*(T_00 + T_m1)
            TC_00 = TC0
            TC_m1 = TC0
            if (TCa /= 0d0) then
               TC_00 = TC_00*pow(rho_avg,TCa)
               !TC_00 = TC_00*pow(prim(i_rho,k),TCa)
               !TC_m1 = TC_m1*pow(prim(i_rho,k-1),TCa)
            end if
            if (TCb /= 0d0) then
               TC_00 = TC_00*pow(T_avg,TCb)
               !TC_00 = TC_00*pow(T_00,TCb)
               !TC_m1 = TC_m1*pow(T_m1,TCb)
            end if
            TC_face = TC_00
            !TC_face = (dr(k-1)*TC_00 + dr(k)*TC_m1)/(dr(k-1) + dr(k))
         end function get_TC_face
         
      end function get_L_face
         
         
      subroutine solve_matrix_equation_with_GMRES(T, deltaT, ierr)
         ! tridiagonal equation info in rhs, sub, diag, sup
         use imex_work, only: newton_iter, model_number, rhs, diag, &
            gmres_mr, gmr_p, gmr_r, gmr_v, gmr_cs, gmr_g, gmr_h, gmr_sn, gmr_y
         real(dp), dimension(:), intent(inout) :: T, deltaT 
         integer, intent(out) :: ierr
         integer :: k, neq, itr_max
         real(dp) :: tol_abs, tol_rel
         include 'formats'
         ierr = 0
         itr_max = 20
         tol_abs = 1.0D-08
         tol_rel = 1.0D-08
         neq = nz ! only 1 variable per zone for this problem
         do k=1,nz
            gmr_p(k) = 1d0/diag(k)
         end do
         call mgmres( &
            neq, matvec, psolve, deltaT, rhs, itr_max, gmres_mr, tol_abs, tol_rel, &
            gmr_r, gmr_v, gmr_cs, gmr_g, gmr_h, gmr_sn, gmr_y)
         do k=1,nz
            T(k) = T(k) + deltaT(k)
         end do            
         
         contains
   
         subroutine psolve(x) ! set x = gmr_p*x
            real(dp), intent(inout) :: x(:) ! (neq)
            integer :: k
            !$omp simd
            do k=1,nz
               x(k) = gmr_p(k)*x(k)
            end do
         end subroutine psolve
   
         subroutine matvec(x, r) ! set r = Jacobian*x
            use imex_work, only: sub, diag, sup, rhs, bp, vp, xp, total_gmres_matvecs
            real(dp), intent(in) :: x(:) ! (neq)
            real(dp), intent(out) :: r(:) ! (neq)
            integer :: k
            total_gmres_matvecs = total_gmres_matvecs + 1
            !$omp simd
            do k=2,nz-1
               r(k) = sub(k-1)*x(k-1) + diag(k)*x(k) + sup(k)*x(k+1)
            end do
            k = 1
            r(k) = diag(k)*x(k) + sup(k)*x(k+1)
            k = nz
            r(k) = sub(k-1)*x(k-1) + diag(k)*x(k)
         end subroutine matvec
         
      end subroutine solve_matrix_equation_with_GMRES


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
            T, prim_0, & ! input
            grad_cell, flux_face, & ! work
            prim_1, v_face, & ! output
            ierr)
         real(dp), intent(in) :: T(:), prim_0(:,:)
         real(dp), intent(out) :: grad_cell(:,:), flux_face(:,:)
         real(dp), intent(out) :: prim_1(:,:), v_face(:)
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         if (thermal_diffusion_only) then
            v_face(:) = 0d0
            return
         end if
         call get_grad_cell(prim_0, grad_cell, ierr)
         if (ierr /= 0) return
         call get_flux_face(prim_0, grad_cell, T, flux_face, v_face, ierr)
         if (ierr /= 0) return
         call get_new_prim(flux_face, prim_0, T, prim_1, ierr)
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
            T, & ! debugging only
            flux_face, v_face, & ! output
            ierr)
         use imex_work, only: area
         real(dp), intent(in) :: prim(:,:), grad_cell(:,:), T(:)
         real(dp), intent(out) ::  flux_face(:,:), v_face(:)
         integer, intent(out) :: ierr
         integer :: j, k, op_err
         include 'formats'
         ierr = 0
         
         if (use_HLLC .and. .not. include_P_in_momentum_flux) &
            stop 'cannot use_HLLC .and. .not. include_P_in_momentum_flux'
         
         !$OMP PARALLEL DO PRIVATE(k,op_err)
         do k=2,nz
            op_err = 0
            if (use_HLLC) then
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
            csL = get_csound(rhoL, TL)
            csR = get_csound(rhoR, TR)
            PgasL = get_Pgas(rhoL, TL)
            PgasR = get_Pgas(rhoR, TR)
            PradL = get_Prad(TL)
            PradR = get_Prad(TR)
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
            alpha = min(abs(vL-csL),abs(vL+csL),abs(vR-csR),abs(vR+csR))
            do j=1,nvar
               flux_face(j,k) = 0.5d0*(fluxR(j)+fluxL(j) - alpha*(primR(j)-primL(j)))
            end do          
         end subroutine get1_flux_face_simple
         
         
         subroutine get1_flux_face_HLLC(k, ierr) 
            ! for P in momentum flux and geometry momentum source term
            use imex_work, only: Cv, dr, imex1, imex2, imex3
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
            rhoL = primL(i_rho)
            sqrt_rhoL = sqrt(rhoL)
            rhoR = primR(i_rho)
            sqrt_rhoR = sqrt(rhoR)
            momL = primL(i_mom)
            momR = primR(i_mom)
            vL = momL/rhoL
            vR = momR/rhoR
            etotL = primL(i_etot)
            etotR = primR(i_etot)
            TL = get_T(rhoL, momL, etotL)
            TR = get_T(rhoR, momR, etotR)
            PgasL = get_Pgas(rhoL, TL)
            PgasR = get_Pgas(rhoR, TR)
            PradL = get_Prad(TL)
            PradR = get_Prad(TR)
            PL = PgasL + PradL
            PR = PgasR + PradR
            csL = get_csound(rhoL, TL)
            csR = get_csound(rhoR, TR)
            fluxL(i_rho) = rhoL*vL
            fluxR(i_rho) = rhoR*vR
            fluxL(i_mom) = momL*vL + PL
            fluxR(i_mom) = momR*vR + PR
            fluxL(i_etot) = (etotL + PL)*vL
            fluxR(i_etot) = (etotR + PR)*vR
            v_face(k) = (sqrt_rhoL*vL + sqrt_rhoR*vR)/(sqrt_rhoL + sqrt_rhoR) 
            if (.true.) then ! debugging
               do j=1,nvar
                  flux_face(j,k) = 0.5d0*(fluxR(j) + fluxL(j))
               end do
            else 
               if (.false.) then ! Sl and Sr using 1988_einfeldt method
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
            end if            
            !imex1(k) = TL ! fluxL(i_rho) ! flux_face(i_rho,k)
            !imex2(k) = TR ! fluxL(i_mom) ! flux_face(i_mom,k)
            !imex3(k) = T(k) ! fluxL(i_etot) ! flux_face(i_etot,k)
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
      
      
      real(dp) function get_csound(rho, T) ! 2020_cheng_shu_song, eqn 2.7
         use imex_work, only: Cv
         real(dp), intent(in) :: rho, T
         real(dp) :: Pgas, Prad, z, Gamma1
         Pgas = get_Pgas(rho, T)
         Prad = get_Prad(T)
         z = Prad/Pgas
         Gamma1 = &
            (gamma/(gamma-1d0) + z*(20d0 * 16d0*z))/((1d0/(gamma-1d0) + 12d0*z)*(1d0 + z))
         get_csound = sqrt(Gamma1*(Pgas + Prad)/rho)
      end function get_csound


      real(dp) function get_egas(rho, T) result(Egas)
         use imex_work, only: Cv
         real(dp), intent(in) :: rho, T
         Egas = rho*Cv*T
      end function get_egas


      real(dp) function get_erad(T) result(Erad)
         use imex_work, only: Cv
         real(dp), intent(in) :: T
         if (low_energy_density_regime) then
            Erad = 0d0
         else
            Erad = crad*pow4(T)
         end if
      end function get_erad


      real(dp) function get_Pgas(rho, T) result(Pgas)
         use imex_work, only: Cv
         real(dp), intent(in) :: rho, T
         Pgas = (gamma - 1d0)*get_egas(rho, T)
      end function get_Pgas


      real(dp) function get_Prad(T) result(Prad)
         real(dp), intent(in) :: T
         Prad = get_erad(T)/3d0
      end function get_Prad


      real(dp) function get_etot(rho, mom, T) result(etot)
         use imex_work, only: Cv
         real(dp), intent(in) :: rho, mom, T
         real(dp) :: v, Egas, Erad, Ekin
         v = mom/rho
         Egas = get_egas(rho, T)
         Erad = get_erad(T)
         Ekin = 0.5d0*mom*v
         etot = Egas + Erad + Ekin
      end function get_etot
      
      
      real(dp) function get_T(rho, mom, etot) result(T) 
         use imex_work, only: Cv
         ! solve for T: etot == rho*Cv*T + 0.5*mom*v + crad*T^4   2020_cheng_shu_song
         real(dp), intent(in) :: rho, mom, etot
         real(dp) :: rhoCv, egas, v, c1, c2, s3, s4, s1, s2, s, sqrt_2s, b1, b1_3, b2, b3, b4, b5
         include 'formats'
         rhoCv = rho*Cv
         v = mom/rho
         if (low_energy_density_regime) then
            egas = etot - 0.5d0*mom*v ! ignore Erad
            T = egas/rhoCv 
            if (is_bad(T) .or. T < 0d0) then
               write(*,*) 'T rho rhoCv mom etot egas', T, rho, rhoCv, mom, etot, egas
               stop 'get_T'
            end if
            return
         end if
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
      
      
      subroutine get_new_prim( &
            flux_face, prim_0, T, & ! input
            prim_1, & ! output
            ierr)
         real(dp), intent(in) :: flux_face(:,:), prim_0(:,:), T(:)
         real(dp), intent(out) :: prim_1(:,:)
         integer, intent(out) :: ierr
         integer :: k, op_err
         ierr = 0
         !$OMP PARALLEL DO PRIVATE(k,op_err)
         do k=1,nz
            op_err = 0
            call get1_new_prim(k, op_err)
            if (op_err /= 0) ierr = op_err
         end do
         !$OMP END PARALLEL DO
         
         contains
         
         subroutine get1_new_prim(k, ierr)
            use imex_work, only: dt, area, dVol, dr, P_face, imex1, imex2, imex3
            use imex_output, only: L
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            integer :: j
            real(dp) :: rho, Pgas, Prad, Ptot, dV
            include 'formats'
            ierr = 0
            rho = prim_0(i_rho,k)
            Pgas = get_Pgas(rho,T(k))
            Prad = get_Prad(T(k))
            Ptot = Pgas + Prad    
            dV = dVol(k)
            if (k == nz) then
               do j=1,nvar
                  prim_1(j,k) = prim_0(j,k) - dt*area(k)*flux_face(j,k)/dV
               end do
               if (include_P_in_momentum_flux) then ! geometry source term
                  prim_1(i_mom,k) = prim_1(i_mom,k) + dt*area(k)*Ptot/dV
               else ! -dP/dr as source term ?
               end if
            else if (k == 1) then
               do j=1,nvar
                  prim_1(j,k) = prim_0(j,k) + dt*area(k+1)*flux_face(j,k+1)/dV
               end do
               if (include_P_in_momentum_flux) then ! geometry source term
                  prim_1(i_mom,k) = prim_1(i_mom,k) - dt*area(k+1)*Ptot/dV
               else ! -dP/dr as source term ?
               end if
            else
               do j=1,nvar
                  prim_1(j,k) = prim_0(j,k) + &
                     dt*(area(k+1)*flux_face(j,k+1) - area(k)*flux_face(j,k))/dV
               end do
               if (include_P_in_momentum_flux) then ! geometry source term
                  prim_1(i_mom,k) = prim_1(i_mom,k) + &
                     dt*(area(k) - area(k+1))*Ptot/dV
               else ! -dP/dr as source term
                  prim_1(i_mom,k) = prim_1(i_mom,k) - dt*(P_face(k) - P_face(k+1))/dr(k)
               end if
            end if     
               
         end subroutine get1_new_prim
         
      end subroutine get_new_prim
      
      
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
         real(dp) :: R_max, deltaR, r00, rp1, T_init, mom_seed
         include 'formats'         
         allocate(r_init(nz), rho_init(nz), mom_init(nz), etot_init(nz))         
         ! grid points equally spaced in r
         R_max = 1d0
         deltaR = (R_max - R_inner)/nz
         rp1 = R_inner
         T_init = 1d0
         mom_seed = 2d-9 ! bigger than this goes unstable for 2nd order
         do k=nz,1,-1
            r00 = rp1 + deltaR
            r_init(k) = r00
            rho_init(k) = 1d0
            mom_init(k) = mom_seed*2d0*(rand() - 0.5d0)
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
         write(*,*) 'initialize_smooth_problem'
      end subroutine initialize_smooth_problem


      
      subroutine initialize_Barenblatt_problem() ! 2001 bates et al "Barenblatt test"
         ! Nonlinear thermal conduction from a point source
         ! hydrodynamic motion disabled
         integer :: k
         real(dp) :: R_max, deltaR, rp1, r00, E0, dV0, rho_ambient, T_ambient, Etot_ambient
         include 'formats'
         allocate(r_init(nz), rho_init(nz), mom_init(nz), etot_init(nz))
         ! grid points equally spaced in r
         R_max = 1d0
         deltaR = (R_max - R_inner)/nz
         rp1 = R_inner
         rho_ambient = 1d0
         T_ambient = 1d-4
         Etot_ambient = get_egas(rho_ambient, T_ambient)
         do k=nz,1,-1
            r00 = rp1 + deltaR
            r_init(k) = r00
            rho_init(k) = rho_ambient
            etot_init(k) = Etot_ambient
            mom_init(k) = 0d0
            rp1 = r00
         end do
         E0 = 10d0
         dV0 = 4d0*pi/3d0*(pow3(r_init(nz)) - pow3(R_inner))
         etot_init(nz) = E0/dV0 
         write(*,*) 'initialize_Barenblatt_problem'
      end subroutine initialize_Barenblatt_problem


      subroutine finish_static_problem()
         write(*,*) 'finish_static_problem'
      end subroutine finish_static_problem


      subroutine finish_Barenblatt_problem()
         use imex_work
         use imex_output
         integer :: k
         include 'formats'
         write(*,*) 'finish_Barenblatt_problem'
         do k=1,nz
            if (L(k) > 0d0) then
               write(*,2) 'r for outermost L /= 0', k, r(k)
               exit
            end if
         end do
         write(*,2) 'T(nz)', nz, T(nz)
         write(*,2) 'model number, dt', model_number, dt
      end subroutine finish_Barenblatt_problem


      subroutine finish_smooth_problem()
         write(*,*) 'finish_smooth_problem'
      end subroutine finish_smooth_problem



      !*****************************************************************************
      !
      !  MGMRES applies restarted GMRES.   derived from MGMRES_ST
      !
      !  Discussion of MGMRES_ST:
      !
      !    The linear system A*X=B is solved iteratively.
      !
      !    The matrix A is assumed to be stored in sparse triplet form.  Only
      !    the nonzero entries of A are stored.  For instance, the K-th nonzero
      !    entry in the matrix is stored by:
      !
      !      A(K) = value of entry,
      !      IA(K) = row of entry,
      !      JA(K) = column of entry.
      !
      !    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
      !    corrections to the code on 31 May 2007.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    13 July 2007
      !
      !  Author:
      !
      !    Original C version by Lili Ju.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
      !    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
      !    Charles Romine, Henk van der Vorst,
      !    Templates for the Solution of Linear Systems:
      !    Building Blocks for Iterative Methods,
      !    SIAM, 1994.
      !    ISBN: 0898714710,
      !    LC: QA297.8.T45.
      !
      !    Tim Kelley,
      !    Iterative Methods for Linear and Nonlinear Equations,
      !    SIAM, 2004,
      !    ISBN: 0898713528,
      !    LC: QA297.8.K45.
      !
      !    Yousef Saad,
      !    Iterative Methods for Sparse Linear Systems,
      !    Second Edition,
      !    SIAM, 2003,
      !    ISBN: 0898715342,
      !    LC: QA188.S17.
      !
      !  Parameters:
      !
      !    Input, integer :: N, the order of the linear system.
      !
      !    Input, integer :: NZ_NUM, the number of nonzero matrix values.
      !
      !    Input, integer :: IA(NZ_NUM), JA(NZ_NUM), the row and column
      !    indices of the matrix values.
      !
      !    Input, real(dp) :: A(NZ_NUM), the matrix values.
      !
      !    Input/output, real(dp) :: X(N); on input, an approximation to
      !    the solution.  On output, an improved approximation.
      !
      !    Input, real(dp) :: RHS(N), the right hand side of the linear system.
      !
      !    Input, integer :: ITR_MAX, the maximum number of (outer)
      !    iterations to take.
      !
      !    Input, integer :: MR, the maximum number of (inner) iterations
      !    to take.  0 < MR <= N.
      !
      !    Input, real(dp) :: TOL_ABS, an absolute tolerance applied to the
      !    current residual.
      !
      !    Input, real(dp) :: TOL_REL, a relative tolerance comparing the
      !    current residual to the initial residual.
      !
      
      subroutine mgmres ( &
            n, matvec, psolve, x, rhs, itr_max, mr, tol_abs, tol_rel, &
            r, v, c, g, h, s, y )
         integer, intent(in) :: n
         interface
            subroutine matvec(x, r) ! set r = A*x
               use const_def, only: dp
               real(dp), intent(in) :: x(:)
               real(dp), intent(out) :: r(:)
            end subroutine matvec
            subroutine psolve(x) ! set x = Precond*x
               use const_def, only: dp
               real(dp), intent(inout) :: x(:)
            end subroutine psolve
         end interface
         real(dp), intent(inout) :: x(:) ! (n)   initial guess on input, result on output
         real(dp), intent(in) :: rhs(:) ! (n)
         real(dp), intent(out), dimension(:) :: r, c, g, s, y
         real(dp), intent(out), dimension(:,:) :: v, h
         integer, intent(in) :: itr_max, mr
         real(dp), intent(in) :: tol_abs, tol_rel
         real(dp) :: av, mu, rho, rho_tol, htmp
         integer :: i, itr, itr_used, j, k, k_copy
         real(dp), parameter :: delta = 1.0D-03
         logical, parameter :: verbose = .false.
         itr_used = 0
         if ( n < mr ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'MGMRES_ST - Fatal error!'
            write ( *, '(a)' ) '  N < MR.'
            write ( *, '(a,i8)' ) '  N = ', n
            write ( *, '(a,i8)' ) '  MR = ', mr
            stop
         end if
         do itr = 1, itr_max ! loop back to here for restarts
            call matvec ( x, r )
            !$omp simd
            do j=1,n
               r(j) = rhs(j) - r(j)
            end do
            call psolve ( r ) ! apply pcond to residual
            rho = sqrt ( dot_product ( r(1:n), r(1:n) ) )
            if ( verbose ) &
               write ( *, '(a,i8,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
            if (is_bad(rho)) stop 'bad residual'
            if ( itr == 1 ) rho_tol = rho * tol_rel
            if ( rho <= rho_tol .and. rho <= tol_abs ) exit
            !$omp simd
            do j=1,n
               v(j,1) = r(j) / rho
            end do
            g(1) = rho
            g(2:mr+1) = 0.0D+00
            h(1:mr+1,1:mr) = 0.0D+00
            k_copy = 1 ! to avoid compiler warnings
            do k = 1, mr
               k_copy = k
               call matvec ( v(1:n,k), v(1:n,k+1) )
               call psolve ( v(1:n,k+1) ) ! apply pcond to result of matvec
               av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
               do j = 1, k
                  h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
                  !$omp simd
                  do i=1,n
                     v(i,k+1) = v(i,k+1) - h(j,k) * v(i,j)
                  end do
               end do
               h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
               if ( av + delta * h(k+1,k) == av ) then
                  do j = 1, k
                     htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
                     h(j,k) = h(j,k) + htmp
                     !$omp simd
                     do i=1,n
                        v(i,k+1) = v(i,k+1) - htmp * v(i,j)
                     end do
                  end do
                  h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
               end if
               if ( h(k+1,k) /= 0.0D+00 ) then
                  !$omp simd
                  do i=1,n
                     v(i,k+1) = v(i,k+1) / h(k+1,k)
                  end do
               end if
               if ( 1 < k ) then
                  !$omp simd
                  do i=1,k+1
                     y(i) = h(i,k)
                  end do
                  do j = 1, k - 1
                     call mult_givens ( c(j), s(j), j, y(1:k+1) )
                  end do
                  !$omp simd
                  do i=1,k+1
                     h(i,k) = y(i)
                  end do
               end if
               mu = sqrt ( pow2(h(k,k)) + pow2(h(k+1,k)) )
               c(k) = h(k,k) / mu
               s(k) = -h(k+1,k) / mu
               h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
               h(k+1,k) = 0.0D+00
               call mult_givens ( c(k), s(k), k, g(1:k+1) )
               rho = abs ( g(k+1) )
               itr_used = itr_used + 1
               if ( verbose ) then
                  write ( *, '(a,i8,a,g14.6)' ) '  K =   ', k, '  Residual = ', rho
               end if
               if ( rho <= rho_tol .and. rho <= tol_abs ) exit
            end do ! k loop
            k = k_copy - 1
            y(k+1) = g(k+1) / h(k+1,k+1)
            do i = k, 1, -1
               y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
            end do
            do i = 1, n
               x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
            end do
            if ( rho <= rho_tol .and. rho <= tol_abs ) exit
            ! else restart
         end do
         if ( verbose ) then
            write ( *, '(a)'       ) ' '
            write ( *, '(a)'       ) 'MGMRES:'
            write ( *, '(a,i8)'    ) '  Iterations = ', itr_used
            write ( *, '(a,g14.6)' ) '  Final residual = ', rho
         end if
         
         contains
         
         subroutine mult_givens ( c, s, k, g )
            real(dp), intent(in) :: c, s ! cos and sin
            integer, intent(in) :: k
            real(dp), intent(inout) :: g(:) ! (1:k+1)
            real(dp) :: g1, g2
            g1 = c * g(k) - s * g(k+1)
            g2 = s * g(k) + c * g(k+1)
            g(k)   = g1
            g(k+1) = g2
         end subroutine mult_givens
         
      end subroutine mgmres
      
      

      !*****************************************************************************
      !
      !! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
      !
      !  Discussion:
      !
      !    In order to make it easier to compare this code with the Original C,
      !    the vector indexing is 0-based.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    08 August 2006
      !
      !  Author:
      !
      !    Original C version by Lili Ju.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Reference:
      !
      !    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
      !    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
      !    Charles Romine, Henk van der Vorst,
      !    Templates for the Solution of Linear Systems:
      !    Building Blocks for Iterative Methods,
      !    SIAM, 1994.
      !    ISBN: 0898714710,
      !    LC: QA297.8.T45.
      !
      !    Tim Kelley,
      !    Iterative Methods for Linear and Nonlinear Equations,
      !    SIAM, 2004,
      !    ISBN: 0898713528,
      !    LC: QA297.8.K45.
      !
      !    Yousef Saad,
      !    Iterative Methods for Sparse Linear Systems,
      !    Second Edition,
      !    SIAM, 2003,
      !    ISBN: 0898715342,
      !    LC: QA188.S17.
      !
      !  Parameters:
      !
      !    Input, real(dp) :: C, S, the cosine and sine of a Givens
      !    rotation.
      !
      !    Input, integer :: K, indicates the location of the first
      !    vector entry.
      !
      !    Input/output, real(dp) :: G(1:K+1), the vector to be modified.
      !    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
      !

      
      
      subroutine test_MGMRES ( )

      !*****************************************************************************
      !
      !! TEST02 tests MGMRES_ST on a 9 by 9 matrix.
      !
      !  Discussion:
      !
      !    A = 
      !      2  0  0 -1  0  0  0  0  0
      !      0  2 -1  0  0  0  0  0  0
      !      0 -1  2  0  0  0  0  0  0
      !     -1  0  0  2 -1  0  0  0  0
      !      0  0  0 -1  2 -1  0  0  0
      !      0  0  0  0 -1  2 -1  0  0
      !      0  0  0  0  0 -1  2 -1  0
      !      0  0  0  0  0  0 -1  2 -1
      !      0  0  0  0  0  0  0 -1  2
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license. 
      !
      !  Modified:
      !
      !    13 July 2007
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer, parameter :: n = 9
        integer, parameter :: mr = n - 1
        integer, parameter :: nz_num = 23

        real(dp), dimension(nz_num) :: a = (/ &
           2.0D+00, -1.0D+00, &
           2.0D+00, -1.0D+00, &
          -1.0D+00,  2.0D+00, &
          -1.0D+00,  2.0D+00, -1.0D+00, &
          -1.0D+00,  2.0D+00, -1.0D+00, &
          -1.0D+00,  2.0D+00, -1.0D+00, &
          -1.0D+00,  2.0D+00, -1.0D+00, &
          -1.0D+00,  2.0D+00, -1.0D+00, &
          -1.0D+00,  2.0D+00 /)
        integer :: i
        integer, dimension(nz_num) :: ia = (/ &
          1, 1, &
          2, 2, &
          3, 3, &
          4, 4, 4, &
          5, 5, 5, &
          6, 6, 6, &
          7, 7, 7, &
          8, 8, 8, &
          9, 9 /)
        integer :: itr_max
        integer, dimension(nz_num) :: ja = (/ &
          1, 4, &
          2, 3, &
          2, 3, &
          1, 4, 5, &
          4, 5, 6, &
          5, 6, 7, &
          6, 7, 8, &
          7, 8, 9, &
          8, 9 /)
        real(dp), dimension(n) :: rhs = (/ &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00 /)
        integer :: seed = 123456789
        integer :: test
        real(dp) :: tol_abs
        real(dp) :: tol_rel
        real(dp) :: x_error
        real(dp) :: x_estimate(n)
        real(dp) :: r(n)
        real(dp) :: v(1:n,1:mr+1)
        real(dp) :: c(1:mr)
        real(dp) :: g(1:mr+1)
        real(dp) :: h(1:mr+1,1:mr)
        real(dp) :: s(1:mr)
        real(dp) :: y(1:mr+1)
        real(dp), dimension(n) :: x_exact = (/ &
          3.5D+00, &
          1.0D+00, &
          1.0D+00, &
          6.0D+00, &
          7.5D+00, &
          8.0D+00, &
          7.5D+00, &
          6.0D+00, &
          3.5D+00 /)

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02'
        write ( *, '(a)' ) '  Test MGMRES_ST on a matrix that is not quite '
        write ( *, '(a,i8)' ) '  the -1,2,-1 matrix, of order N = ', n

        do test = 1, 2

          if ( test == 1 ) then

            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) '  First try, use zero initial vector:'
            write ( *, '(a)' ) ' '

            x_estimate(1:n) = 0.0D+00

          else

            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) '  Second try, use random initial vector:'
            write ( *, '(a)' ) ' '

            call r8vec_uniform_01 ( n, seed, x_estimate )

          end if

          x_error = sqrt ( sum ( pow2( x_exact(1:n) - x_estimate(1:n) ) ) )

          write ( *, '(a,g14.6)' ) '  Before solving, X_ERROR = ', x_error

          itr_max = 20
          tol_abs = 1.0D-08
          tol_rel = 1.0D-08

          call mgmres ( &
             n, ax, psolve, x_estimate, rhs, itr_max, mr, tol_abs, tol_rel, &
             r, v, c, g, h, s, y )

          x_error = sqrt ( sum ( pow2( x_exact(1:n) - x_estimate(1:n) ) ) )

          write ( *, '(a,g14.6)' ) '  After solving, X_ERROR = ', x_error

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Final solution estimate:'
          write ( *, '(a)' ) ' '
          do i = 1, n
            write ( *, '(2x,i8,2x,g14.6)' ) i, x_estimate(i)
          end do

        end do
        
        
        contains
        
        
        subroutine psolve (x)
           real(dp), intent(inout) :: x(:)
        end subroutine psolve
        
        
        subroutine ax ( x, w )

        !*****************************************************************************
        !
        !! AX_ST computes A*x for a matrix stored in sparset triplet form.
        !
        !  Discussion:
        !
        !    The matrix A is assumed to be sparse.  To save on storage, only
        !    the nonzero entries of A are stored.  For instance, the K-th nonzero
        !    entry in the matrix is stored by:
        !
        !      A(K) = value of entry,
        !      IA(K) = row of entry,
        !      JA(K) = column of entry.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    08 August 2006
        !
        !  Author:
        !
        !    Original C version by Lili Ju.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Reference:
        !
        !    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
        !    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
        !    Charles Romine, Henk van der Vorst,
        !    Templates for the Solution of Linear Systems:
        !    Building Blocks for Iterative Methods,
        !    SIAM, 1994.
        !    ISBN: 0898714710,
        !    LC: QA297.8.T45.
        !
        !    Tim Kelley,
        !    Iterative Methods for Linear and Nonlinear Equations,
        !    SIAM, 2004,
        !    ISBN: 0898713528,
        !    LC: QA297.8.K45.
        !
        !    Yousef Saad,
        !    Iterative Methods for Sparse Linear Systems,
        !    Second Edition,
        !    SIAM, 2003,
        !    ISBN: 0898715342,
        !    LC: QA188.S17.
        !
        !  Parameters:
        !
        !    Input, integer :: N, the order of the system.
        !
        !    Input, integer :: NZ_NUM, the number of nonzeros.
        !
        !    Input, integer :: IA(NZ_NUM), JA(NZ_NUM), the row and column
        !    indices of the matrix values.
        !
        !    Input, real(dp) :: A(NZ_NUM), the matrix values.
        !
        !    Input, real(dp) :: X(N), the vector to be multiplied by A.
        !
        !    Output, real(dp) :: W(N), the value of A*X.
        !

          real(dp), intent(in) :: x(:)
          real(dp), intent(out) :: w(:)

          integer :: i
          integer :: j
          integer :: k

          w(1:n) = 0.0D+00

          do k = 1, nz_num
            i = ia(k)
            j = ja(k)
            w(i) = w(i) + a(k) * x(j)
          end do

        end subroutine ax


        subroutine r8vec_uniform_01 ( n, seed, r )

        !*****************************************************************************
        !
        !! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
        !
        !  Discussion:
        !
        !    An R8VEC is a vector of real ( kind = 8 ) values.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    05 July 2006
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Reference:
        !
        !    Paul Bratley, Bennett Fox, Linus Schrage,
        !    A Guide to Simulation,
        !    Second Edition,
        !    Springer, 1987,
        !    ISBN: 0387964673,
        !    LC: QA76.9.C65.B73.
        !
        !    Bennett Fox,
        !    Algorithm 647:
        !    Implementation and Relative Efficiency of Quasirandom
        !    Sequence Generators,
        !    ACM Transactions on Mathematical Software,
        !    Volume 12, Number 4, December 1986, pages 362-376.
        !
        !    Pierre L'Ecuyer,
        !    Random Number Generation,
        !    in Handbook of Simulation,
        !    edited by Jerry Banks,
        !    Wiley, 1998,
        !    ISBN: 0471134031,
        !    LC: T57.62.H37.
        !
        !    Peter Lewis, Allen Goodman, James Miller,
        !    A Pseudo-Random Number Generator for the System/360,
        !    IBM Systems Journal,
        !    Volume 8, Number 2, 1969, pages 136-143.
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) N, the number of entries in the vector.
        !
        !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
        !    should NOT be 0.  On output, SEED has been updated.
        !
        !    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
        !

          integer:: n

          integer:: i
          integer:: k
          integer:: seed
          real(dp) :: r(:)

          if ( seed == 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
            write ( *, '(a)' ) '  Input value of SEED = 0.'
            stop
          end if

          do i = 1, n

            k = seed / 127773

            seed = 16807 * ( seed - k * 127773 ) - k * 2836

            if ( seed < 0 ) then
              seed = seed + 2147483647
            end if

            r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

          end do

        end subroutine r8vec_uniform_01
        

      end subroutine test_MGMRES
      
      
   end module imex
   