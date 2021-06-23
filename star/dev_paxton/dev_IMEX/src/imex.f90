 
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
         time, dt, dt_next, gam, delta, Cv
      real(dp), dimension(:), allocatable :: &
         sub, diag, sup, bp, vp, xp, &
         r, area, Vol, dr, dr_bar, dVol, & 
         Eeos_x, rhs, deltaT, P_face, L_start, &
         T_start, T_1, d_ETOT_dt_I1, d_ETOT_dt_I2, &
         imex1, imex2, imex3, imex4, imex5, imex6
      real(dp), dimension(:,:), allocatable :: &
         prim_start, cons_start, prim_1, cons_1, &
         grad_cell, flux_face, d_cons_dt_X0, d_cons_dt_X1
      real(dp) :: dt_advection, dt_grid, dt_max_new, dt_front_v, total_KE
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
         if (crad < 0d0) crad = 7.5657332502799993d-015 ! crad in mesa
         if (cgas < 0d0) cgas = 8.314462618d7 ! cgas in mesa
         Cv = cgas/(gamma - 1d0)
         write(*,1) 'cgas', cgas
         write(*,1) 'gamma', gamma
         write(*,1) 'Cv', Cv
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
         call alloc_work_arrays()
         call set_grid_vars()
         call set_init_prim_cons_T()
         total_energy_initial = sum(cons(i_etot,1:nz))
         !stop 'start_imex'
                  
         contains
         
         subroutine set_init_prim_cons_T()
            integer :: k, j
            real(dp) :: &
               sum_EKIN_actual, sum_EKIN_for_ETOT, sum_EGAS, sum_ERAD, sum_ETOT
            include 'formats'
            do k=1,nz
               !write(*,2) 'set_init_prim_cons_T rho_init', k, rho_init(k)
               T(k) = get_T(rho_init(k), mom_init(k), etot_init(k))
               prim(i_rho,k) = rho_init(k)
               prim(i_mom,k) = mom_init(k)
               prim(i_etot,k) = etot_init(k)
               do j=1,nvar
                  cons(j,k) = prim(j,k)*dVol(k)
               end do
               L(k) = 0d0
               v_face(k) = 0d0
               imex1(k) = rho_init(k)
               imex2(k) = mom_init(k)
               imex3(k) = etot_init(k)
               imex4(k) = T(k)
               !write(*,2) 'r T rho etot', k, r(k), T(k), rho_init(k), etot_init(k)
            end do
            call get_imex_total_energies(prim, T, &
               sum_EKIN_actual, sum_EKIN_for_ETOT, sum_EGAS, sum_ERAD, sum_ETOT)
            !write(*,1) 'initial sum_EKIN_actual', sum_EKIN_actual
            !write(*,1) 'initial sum_EKIN_for_ETOT', sum_EKIN_for_ETOT
            !write(*,1) 'initial sum_EGAS', sum_EGAS
            !write(*,1) 'initial sum_ERAD', sum_ERAD
            !write(*,1) 'initial sum_ETOT', sum_ETOT  
            !write(*,2) 'prim(i_etot,nz)', nz, prim(i_etot,nz)
            !write(*,2) 'cons(i_etot,nz)', nz, cons(i_etot,nz)
            !write(*,2) 'T', nz, T(nz)
            !write(*,2) 'T', nz-1, T(nz-1)
            !write(*,*)          
         end subroutine set_init_prim_cons_T
         
         subroutine alloc_work_arrays()
            allocate(&
               r(nz), area(nz), Vol(nz), dr(nz), dr_bar(nz), dVol(nz), & 
               Eeos_x(nz), rhs(nz), sub(nz), diag(nz), sup(nz), deltaT(nz), &
               bp(nz), vp(nz), xp(nz), P_face(nz), L_start(nz), &
               T_start(nz), T_1(nz), T(nz), L(nz), v_face(nz), csound(nz), &
               d_ETOT_dt_I1(nz), d_ETOT_dt_I2(nz), &
               prim_start(nvar,nz), cons_start(nvar,nz), &
               prim_1(nvar,nz), cons_1(nvar,nz), &
               prim(nvar,nz), cons(nvar,nz), &
               imex1(nz), imex2(nz), imex3(nz), imex4(nz), imex5(nz), imex6(nz), &
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
            Pgas = get_Pgas(rho,T(k))
            Prad = get_Prad(T(k))
            P_out(k) = Pgas + Prad
            v_out(k) = v
            T_out(k) = T(k)
            L_out(k) = L(k)
            cs_out(k) = csound(k)
            v_div_cs_out(k) = v_out(k)/cs_out(k)
         end do
      end subroutine get_imex_data
      
      
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
                  cons_start(j,k) = cons(j,k)
               end do
            end do
         end subroutine save_start_values
         
         subroutine report_energies()
            real(dp) :: &
               sum_EKIN_actual, sum_EKIN_for_ETOT, sum_EGAS, sum_ERAD, sum_ETOT
            include 'formats'
            call get_imex_total_energies(prim, T, &
               sum_EKIN_actual, sum_EKIN_for_ETOT, sum_EGAS, sum_ERAD, sum_ETOT)
            !write(*,1) 'sum_EKIN_actual', sum_EKIN_actual
            !write(*,1) 'sum_EKIN_for_ETOT', sum_EKIN_for_ETOT
            !write(*,1) 'sum_EGAS', sum_EGAS
            !write(*,1) 'sum_ERAD', sum_ERAD
            !write(*,1) 'sum_ETOT', sum_ETOT  
            !write(*,*)          
         end subroutine report_energies
         
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
         if (time_centering) call save_start_values()
         !write(*,2) 'do_step dt time', model_number, dt, time
         if (dt == 0d0) then
            write(*,2) 'dt == 0d0', model_number, dt, time
            stop 'do_step'
         end if
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
            ! 1997_ascher_ruuth_spiteri L-stable, 2 stage, 2nd order (section 2.6)
            ! also try L-stable, 3 stage implicit, 4 stage explicit, 3rd order DIRK (section 2.7)
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
      
      
      subroutine save_start_values()
         use imex_output, only: L
         use imex_work, only: L_start
         integer :: k
         do k=1,nz
            L_start(k) = L(k)
         end do
      end subroutine save_start_values
      
      
      subroutine stage1( &
            T_0, prim_0, cons_0, & ! input
            grad_cell, flux_face, & ! work
            d_cons_dt_X0, d_ETOT_dt_I1, & ! output
            T_1, prim_1, cons_1, L_1, v_face_1, & ! output
            ierr)
         use imex_work, only: dt, gam, dVol, model_number
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
         if (ierr /= 0) then
            write(*,2) 'stage1 calc_explicit failed', model_number
            return
         end if
         do k=1,nz ! update explicit intermediate result
            do j=1,nvar
               cons_1(j,k) = cons_0(j,k) + dt*gam*d_cons_dt_X0(j,k)
               prim_1(j,k) = cons_1(j,k)/dVol(k)
            end do
         end do
         call calc_implicit(T_0, prim_1, cons_1, d_ETOT_dt_I1, T_1, L_1, ierr)
         if (ierr /= 0) then
            write(*,2) 'stage1 calc_implicit failed', model_number
            return
         end if
         k = nz
         !write(*,2) 'ETOT dETOT T', nz, cons_1(i_etot,k), dt*gam*d_ETOT_dt_I1(k), T_1(k)
         !stop 'stage1'
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
         use imex_work, only: dt, gam, delta, dVol, model_number
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
         if (ierr /= 0) then
            write(*,2) 'stage2 calc_explicit failed', model_number
            return
         end if
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
         if (ierr /= 0) then
            write(*,2) 'stage2 calc_implicit failed', model_number
            return
         end if
         do k=1,nz ! update etot
            cons_2(i_etot,k) = cons_2(i_etot,k) + dt*gam*d_ETOT_dt_I2(k)
            prim_2(i_etot,k) = cons_2(i_etot,k)/dVol(k)
         end do
      end subroutine stage2
            
      
      real(dp) function pick_next_timestep(T_2, prim_2, T_0, prim_0)
         use imex_work, only: dt, &
            dt_advection, dt_grid, dt_max_new, dt_front_v
         real(dp), intent(in), dimension(:) :: T_2, T_0
         real(dp), intent(in), dimension(:,:) ::  prim_2, prim_0
         include 'formats'
         dt_advection = get_min_dt_advection(T_2, prim_2)
         dt_grid = get_min_dt_grid()
         dt_front_v = get_dt_front_v(prim_0, prim_2)
         dt_max_new = dt*max_timestep_factor
         pick_next_timestep = min(dt_advection, dt_grid, dt_max_new, dt_front_v)
      end function pick_next_timestep
      
      
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
         !write(*,2) 'dt_front_v sum_dE_dt sum_dE', model_number, dt_front_v, sum_dE_dt, sum_dE
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
      
      
!!! implicit part      
      
      
      subroutine calc_implicit( &
            T_prev, prim, cons, & ! input
            d_ETOT_dt_I, T, L, & ! output
            ierr)
         use imex_work, only: newton_iter, stage, Eeos_x, rhs, dVol, model_number
         real(dp), intent(in) :: T_prev(:)
         real(dp), intent(in), dimension(:,:) :: prim, cons
         real(dp), intent(out), dimension(:) :: d_ETOT_dt_I, T, L ! output
         integer, intent(out) :: ierr       
         integer :: k
         real(dp) :: rho, v, etot, ekin, resid_norm, resid_max
         logical :: converged
         include 'formats'        
         ierr = 0 
         if (TC0 == 0d0) then ! no thermal diffusion
            do k=1,nz
               d_ETOT_dt_I(k) = 0d0
            end do
            return
         end if
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
               call store_matrix_equation(T,T_prev,k) ! sets d_ETOT_dt_I(k) and L(k)
            end do
            !$OMP END PARALLEL DO
            resid_norm = sum(abs(rhs(1:nz)))/nz
            resid_max = maxval(abs(rhs(1:nz)))
            if (is_bad(resid_norm) .or. is_bad(resid_max)) then
               write(*,4) 'resid norm max', newton_iter, stage, model_number, resid_norm, resid_max
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
            call solve_matrix_equation(ierr) ! updates T
         end do         
         newton_iter = 0
         if (.not. converged) then
            write(*,2) 'calc_implicit failed to converge', model_number, resid_norm, resid_max
            ierr = -1
            return
         end if
         
         contains
         
         subroutine store_matrix_equation(T,T_prev,k)
            use imex_work, only: dt, Cv, gam, dVol, sub, diag, sup, rhs
            integer, intent(in) :: k
            real(dp), intent(in) :: T(:), T_prev(:)
            type(auto_diff_real_4var_order1) :: &
               T_00, L_00, L_p1, dL, Eeos_expected, Eeos_actual, residual
            real(dp) :: rho
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
            d_ETOT_dt_I(k) = -dL%val
            Eeos_expected = Eeos_X(k) - dt*gam*dL
            if (is_bad(Eeos_expected%val)) then
               !$omp critical (store_matrix_equation_crit1)
               write(*,2) 'Eeos_expected%val', k, Eeos_expected%val
               write(*,2) 'Eeos_X(k)', k, Eeos_X(k)
               write(*,2) 'L_00%val', k, L_00%val
               write(*,2) 'L_p1%val', k, L_p1%val
               write(*,2) 'dt', k, dt
               write(*,2) 'gam', k, gam
               write(*,*) 'use_T_prev_for_L', use_T_prev_for_L
               write(*,2) 'T_prev', k, T_prev(k)
               if (k < nz) write(*,2) 'T_prev', k+1, T_prev(k+1)
               stop 'store_matrix_equation'
               !$omp end critical (store_matrix_equation_crit1)
            end if
            rho = prim(i_rho,k)
            Eeos_actual = rho*Cv*T_00*dVol(k)
            if (.not. low_energy_density_regime) &
               Eeos_actual = Eeos_actual + crad*pow4(T_00)*dVol(k)
            if (is_bad(Eeos_actual%val)) then
               !$omp critical (store_matrix_equation_crit2)
               write(*,2) 'Eeos_actual%val', k, Eeos_actual%val
               write(*,2) 'T_00%val', k, T_00%val
               write(*,2) 'L_p1%val', k, L_p1%val
               write(*,2) 'dVol(k)', k, dVol(k)
               write(*,2) 'gam', k, gam
               stop 'store_matrix_equation'
               !$omp end critical (store_matrix_equation_crit2)
            end if
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
            !write(*,2) 'T', nz, T(nz)
         end subroutine solve_matrix_equation

      end subroutine calc_implicit
      
      
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
         if (thermal_diffusion_only) then
            d_cons_dt_X(:,:) = 0d0
            v_face(:) = 0d0
            return
         end if
         call get_grad_cell(prim, grad_cell, ierr)
         if (ierr /= 0) return
         call get_flux_face(prim, grad_cell, T, flux_face, v_face, ierr)
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
      
      
      subroutine get_d_cons_dt_X( &
            flux_face, prim, T, & ! input
            d_cons_dt_X, & ! output
            ierr)
         real(dp), intent(in) :: flux_face(:,:), prim(:,:), T(:)
         real(dp), intent(out) :: d_cons_dt_X(:,:)
         integer, intent(out) :: ierr
         integer :: k, op_err
         ierr = 0
         !$OMP PARALLEL DO PRIVATE(k,op_err)
         do k=1,nz
            op_err = 0
            call get1_d_cons_dt_X(k, op_err)
            if (op_err /= 0) ierr = op_err
         end do
         !$OMP END PARALLEL DO
         
         contains
         
         subroutine get1_d_cons_dt_X(k, ierr)
            use imex_work, only: area, dVol, dr, P_face, imex1, imex2, imex3
            use imex_output, only: L
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            integer :: j
            real(dp) :: Pgas, Prad, Ptot
            include 'formats'
            ierr = 0
            Pgas = get_Pgas(prim(i_rho,k),T(k))
            Prad = get_Prad(T(k))
            Ptot = Pgas + Prad            
            if (k == nz) then
               do j=1,nvar
                  d_cons_dt_X(j,k) = -area(k)*flux_face(j,k)
               end do
               if (include_P_in_momentum_flux) then
                  !imex1(k) = d_cons_dt_X(i_mom,k)
                  !imex2(k) = area(k)*Ptot
                  d_cons_dt_X(i_mom,k) = d_cons_dt_X(i_mom,k) + area(k)*Ptot
                  !imex3(k) = d_cons_dt_X(i_mom,k)
               else ! -dP/dr as source term
               end if
            else if (k == 1) then
               do j=1,nvar
                  d_cons_dt_X(j,k) = area(k+1)*flux_face(j,k+1)
               end do
               if (include_P_in_momentum_flux) then
                  !imex1(k) = d_cons_dt_X(i_mom,k)
                  !imex2(k) = - area(k+1)*Ptot
                  d_cons_dt_X(i_mom,k) = d_cons_dt_X(i_mom,k) - area(k+1)*Ptot
                  !imex3(k) = d_cons_dt_X(i_mom,k)
               else ! -dP/dr as source term
               end if
            else
               do j=1,nvar
                  d_cons_dt_X(j,k) = area(k+1)*flux_face(j,k+1) - area(k)*flux_face(j,k)
               end do
               if (include_P_in_momentum_flux) then
                  !imex1(k) = d_cons_dt_X(i_mom,k)
                  !imex2(k) = (area(k) - area(k+1))*Ptot
                  d_cons_dt_X(i_mom,k) = d_cons_dt_X(i_mom,k) + &
                     (area(k) - area(k+1))*Ptot ! geometry source term
                  !imex3(k) = d_cons_dt_X(i_mom,k)
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


   end module imex
   