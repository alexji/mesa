 
   module imex

      !use star_lib
      !use star_def
      use const_def
      use math_lib
      use auto_diff
      use utils_lib, only: is_bad
      
      implicit none

      integer, parameter :: &
         i_rho = 1, i_mom = 2, i_etot = 3, nvar = i_etot
      
      ! input
      integer :: nz ! number of zones
      real(dp), dimension(:), allocatable :: r ! (nz)
      real(dp), dimension(:), allocatable :: rho, mom, T ! (nz)
      
      ! boundaries
      real(dp) :: R_inner, M_inner, L_inner ! values at inner boundary
      real(dp) :: Lsurf_factor
         ! L_surf = Lsurf_factor * area(1) * clight * crad * T(1)^4

      integer :: model_number ! number of steps so far
      real(dp) :: time ! current time
      real(dp) :: time_end ! final time for entire run.  stop when reach this.
      real(dp) :: dt_next ! suggested next timestep
      
      ! eos and kap
      real(dp) :: gamma ! for gamma law ideal gas
      real(dp) :: Rgas ! gas constant 
      real(dp) :: Cv ! = Rgas/(gamma - 1d0)
      real(dp) :: D ! diffusion coeffient for radiative energy
      
      ! solver
      real(dp) :: CFL ! Courant-Friedrichs_Lewy control parameter for timestep
      real(dp) :: atol_resid_norm, rtol_resid_norm
      integer :: newton_iter_max
      
      ! work
      integer :: stage, newton_iter
      real(dp) :: resid_norm, tol_resid_norm, dt, L_surf
      real(dp) :: gam, delta
      real(dp), dimension(:), allocatable :: &
         area, Vol, dr, dr_bar, dVol, & 
         rhs, sub, diag, sup, deltaT, bp, vp, xp, &
         T_start, T_1, L, v_bar, d_ETOT_dt_I1, d_ETOT_dt_I2
      real(dp), dimension(:,:), allocatable :: &
         prim_start, cons_start, prim_1, cons_1, prim, cons, &
         grad_cell, flux_face, d_cons_dt_X0, d_cons_dt_X1
         
      contains
      
      
!!! top level control      

      
      subroutine start_imex(init_problem, ierr)
         interface
            subroutine init_problem()
            end subroutine init_problem
         end interface
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         write(*,*) 'start_imex'
         call init_problem() ! set all input variables listed above
         call alloc_work_arrays()
         call set_grid_vars(r)
         ! 2016_Wang_Shu_Zang; 1997_ascher_ruuth_spiteri
         gam = 1d0 - sqrt(2d0)/2d0
         delta = 1d0 - 1d0/(2d0*gam)
         
         contains
         
         subroutine alloc_work_arrays()
            allocate(&
               area(nz), Vol(nz), dr(nz), dr_bar(nz), dVol(nz), & 
               rhs(nz), sub(nz), diag(nz), sup(nz), deltaT(nz), &
               bp(nz), vp(nz), xp(nz), &
               T_start(nz), T_1(nz), L(nz), v_bar(nz), &
               d_ETOT_dt_I1(nz), d_ETOT_dt_I2(nz), &
               prim_start(nvar,nz), cons_start(nvar,nz), &
               prim_1(nvar,nz), cons_1(nvar,nz), &
               prim(nvar,nz), cons(nvar,nz), &
               grad_cell(nvar,nz), flux_face(nvar,nz), &
               d_cons_dt_X0(nvar,nz), d_cons_dt_X1(nvar,nz))
         end subroutine alloc_work_arrays
      
         subroutine set_grid_vars(r)
            real(dp), intent(in) :: r(:)
            ! use r to set area, Vol, dr, dr_bar, dVol
            integer :: k
            real(dp) :: rm1, r00, volm1, vol00
            include 'formats'
            rm1 = R_inner
            volm1 = 4d0*pi/3d0*pow3(rm1)
            do k=nz,1,-1
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
            dr_bar(1) = 0.5d0*dr(k) ! NOT USED
         end subroutine set_grid_vars
         
      end subroutine start_imex
      
      
      subroutine finish_imex(ierr)
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         write(*,*) 'finish_imex'
      end subroutine finish_imex

      
      subroutine steps_imex( &
            max_steps_for_this_call, age, timestep, &
            final_step, mod_number, num_zones, ierr)
         integer, intent(in) :: max_steps_for_this_call
         real(dp), intent(out) :: age, timestep
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
               T, prim, cons, L, v_bar, & ! output
               ierr)
            if (ierr /= 0 .or. final_step) exit
         end do
         age = time
         timestep = dt
         num_zones = nz
         mod_number = model_number
         
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
            T_2, prim_2, cons_2, L, v_bar, & ! output
            ierr) result(final_step)
         real(dp), intent(in), dimension(:) :: T_0 ! input
         real(dp), intent(in), dimension(:,:) ::  cons_0, prim_0 ! input
         real(dp), intent(out), dimension(:) :: & ! work
            T_1, d_ETOT_dt_I1,  d_ETOT_dt_I2
         real(dp), intent(out), dimension(:,:) :: & ! work
            prim_1, cons_1, grad_cell, flux_face, &
            d_cons_dt_X0, d_cons_dt_X1
         real(dp), intent(out), dimension(:) :: T_2, L, v_bar ! output
         real(dp), intent(out), dimension(:,:) :: prim_2, cons_2 ! output
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         newton_iter = 0
         dt = dt_next
         final_step = (time + dt >= time_end)
         if (final_step) dt = time_end - time
         stage = 1
         call stage1( &
            T_0, prim_0, cons_0, & ! input
            grad_cell, flux_face, & ! work
            d_cons_dt_X0, d_ETOT_dt_I1, & ! output
            T_1, prim_1, cons_1, L, v_bar, & ! output
            ierr)
         if (ierr /= 0) then
            write(*,2) 'stage1 failed', model_number
            return
         end if
         stage = 2
         call stage2( &
            cons_0, T_1, prim_1, cons_1, d_cons_dt_X0, d_ETOT_dt_I1, & ! input
            grad_cell, flux_face, d_cons_dt_X1, d_ETOT_dt_I2, & ! work
            T_2, prim_2, cons_2, L, v_bar, & ! output
            ierr)
         if (ierr /= 0) then
            write(*,2) 'stage2 failed', model_number
            return
         end if
         time = time + dt
         stage = 0
         dt_next = pick_next_timestep(T_2, prim_2, T_0, prim_0)
      end function do_step
      
      
      subroutine stage1( &
            T_0, prim_0, cons_0, & ! input
            grad_cell, flux_face, & ! work
            d_cons_dt_X0, d_ETOT_dt_I1, & ! output
            T_1, prim_1, cons_1, L_1, v_bar_1, & ! output
            ierr)
         real(dp), intent(in), dimension(:) :: T_0
         real(dp), intent(in), dimension(:,:) ::  cons_0, prim_0
         real(dp), intent(out), dimension(:) :: &
            d_ETOT_dt_I1, T_1, L_1, v_bar_1
         real(dp), intent(out), dimension(:,:) :: &
            grad_cell, flux_face, d_cons_dt_X0, prim_1, cons_1
         integer, intent(out) :: ierr
         integer :: k, j
         include 'formats'
         ierr = 0
         call calc_explicit( &
            T_0, prim_0, cons_0, grad_cell, flux_face, d_cons_dt_X0, v_bar_1, ierr)
         if (ierr /= 0) return
         do k=1,nz ! update explicit intermediate result
            do j=1,nvar
               cons_1(j,k) = cons_0(j,k) + dt*gam*d_cons_dt_X0(j,k)
               prim_1(j,k) = cons_1(j,k)/dVol(k)
            end do
         end do
         call calc_implicit(T_0, prim_1, cons_1, &
            sub, diag, sup, rhs, deltaT, &
            d_ETOT_dt_I1, T_1, L_1, ierr)
         if (ierr /= 0) return
         do k=1,nz ! update etot
            cons_1(i_etot,k) = cons_1(i_etot,k) + dt*gam*d_ETOT_dt_I1(k)
            prim_1(i_etot,k) = cons_1(i_etot,k)/dVol(k)
         end do
      end subroutine stage1
      
      
      subroutine stage2( &
            cons_0, T_1, prim_1, cons_1, d_cons_dt_X0, d_ETOT_dt_I1, & ! input
            grad_cell, flux_face, d_cons_dt_X1, d_ETOT_dt_I2, & ! work
            T_2, prim_2, cons_2, L_2, v_bar_2, & ! output
            ierr)
         real(dp), intent(in), dimension(:) :: &
            T_1, d_ETOT_dt_I1
         real(dp), intent(in), dimension(:,:) :: &
            cons_0, prim_1, cons_1, d_cons_dt_X0
         real(dp), intent(out), dimension(:) :: &
            d_ETOT_dt_I2, T_2, L_2, v_bar_2
         real(dp), intent(out), dimension(:,:) :: &
            grad_cell, flux_face, d_cons_dt_X1, prim_2, cons_2
         integer, intent(out) :: ierr
         integer :: k, j
         include 'formats'
         ierr = 0
         call calc_explicit( &
            T_1, prim_1, cons_1, grad_cell, flux_face, d_cons_dt_X1, v_bar_2, ierr)
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
         call calc_implicit(T_1, prim_2, cons_2, &
            sub, diag, sup, rhs, deltaT, &
            d_ETOT_dt_I2, T_2, L_2, ierr)
         if (ierr /= 0) return
         do k=1,nz ! update etot
            cons_2(i_etot,k) = cons_2(i_etot,k) + dt*gam*d_ETOT_dt_I2(k)
            prim_2(i_etot,k) = cons_2(i_etot,k)/dVol(k)
         end do
      end subroutine stage2
            
      
      real(dp) function pick_next_timestep(T_2, prim_2, T_0, prim_0)
         real(dp), intent(in), dimension(:) :: T_2, T_0
         real(dp), intent(in), dimension(:,:) ::  prim_2, prim_0
         real(dp) :: dt_advection, dt_Etot
         include 'formats'
         dt_advection = get_min_dt_advection(T_2, prim_2)
         dt_Etot = get_min_dt_Etot(prim_0, prim_2)
         pick_next_timestep = min(dt_advection, dt_Etot)
      end function pick_next_timestep
      
      
      real(dp) function get_min_dt_advection(T, prim) result(dt_advection)
         real(dp), intent(in) :: T(:)
         real(dp), intent(in) :: prim(:,:)
         real(dp) :: dt_cell, csound, vel
         integer :: k
         dt_advection = 1d99
         do k=1,nz
            csound = sqrt(Rgas*T(k))
            vel = prim(i_mom,k)/prim(i_rho,k)
            dt_cell = dr(k)/abs(vel + csound)
            if (dt_cell < dt_advection) dt_advection = dt_cell
         end do
         dt_advection = CFL*dt_advection
      end function get_min_dt_advection
      
      
      real(dp) function get_min_dt_Etot(prim_start, prim_final) result(dt_Etot)
         ! W.J. Rider, D.A. Knoll,
         ! Time step size selection for radiation diffusion calculations,
         ! J. Comput. Phys. 152-2 (1999) 790–795.
         real(dp), intent(in), dimension(:,:) :: prim_start, prim_final
         real(dp) :: sum_abs_dEtot_time, sum_abs_dEtot_grid
         integer :: k
         sum_abs_dEtot_time = 0d0
         do k=1,nz
            sum_abs_dEtot_time = sum_abs_dEtot_time + &
               abs(prim_final(i_etot,k) - prim_start(i_etot,k))
         end do
         sum_abs_dEtot_grid = 0d0
         do k=2,nz-1
            sum_abs_dEtot_grid = sum_abs_dEtot_grid + &
               abs(prim_final(i_etot,k-1) - prim_final(i_etot,k+1))
         end do
         sum_abs_dEtot_grid = 0.5d0*sum_abs_dEtot_grid
         dt_Etot = dt*sum_abs_dEtot_grid/sum_abs_dEtot_time
      end function get_min_dt_Etot
      
      
!!! implicit part      
      
      
      subroutine calc_implicit( &
            T_prev, prim, cons, & ! input
            sub, diag, sup, rhs, deltaT, & ! work
            d_ETOT_dt_I, T, L, & ! output
            ierr)
         real(dp), intent(in) :: T_prev(:)
         real(dp), intent(in), dimension(:,:) :: prim, cons
         real(dp), intent(out), dimension(:) :: & ! work
            sub, diag, sup, rhs, deltaT
         real(dp), intent(out), dimension(:) :: d_ETOT_dt_I, T, L ! output
         integer, intent(out) :: ierr       
         integer :: k
         real(dp) :: resid_norm, tol_resid_norm
         logical :: converged
         include 'formats'        
         ierr = 0 
         do k=1,nz ! initial guess
            T(k) = T_prev(k)
         end do         
         converged = .false.
         do newton_iter = 1, newton_iter_max         
            !$OMP PARALLEL DO PRIVATE(k)
            do k=1,nz
               call store_matrix_equation(k) ! sets d_ETOT_dt_I(k) and L(k)
            end do
            !$OMP END PARALLEL DO
            resid_norm = sqrt(sum(pow2(rhs(1:nz))))
            if (newton_iter == 1) then
               tol_resid_norm = min(atol_resid_norm, rtol_resid_norm*resid_norm)
            else if (resid_norm <= tol_resid_norm) then
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
            integer, intent(in) :: k
            type(auto_diff_real_4var_order1) :: &
               T_00, T4_m1, T4_00, T4_p1, L_00, L_p1, &
               Eeos_expected, Eeos_actual, residual
            real(dp) :: Etot_X, rho, v, etot, ekin, Eeos_X
            include 'formats'
            T_00 = wrap_T_00(T,k)
            T4_00 = pow4(T_00)
            if (k > 1) then
               T4_m1 = pow4(wrap_T_m1(T,k))
               L_00 = -area(k)*D*(T4_m1 - T4_00)/dr_bar(k)
            else ! outer boundary
               L_surf = Lsurf_factor * area(k) * clight * crad * T4_00%val
               L_00 = L_surf
            end if
            L(k) = L_00%val
            if (k < nz) then
               T4_p1 = pow4(wrap_T_p1(T,k))
               L_p1 = -area(k+1)*D*(T4_00 - T4_p1)/dr_bar(k+1)
            else ! inner boundary
               L_p1 = L_inner
            end if
            d_ETOT_dt_I(k) = L_p1%val - L_00%val
            rho = prim(i_rho,k)
            v = prim(i_mom,k)/rho
            etot = prim(i_etot,k)
            ekin = 0.5d0*rho*pow2(v)
            Eeos_X = (etot - ekin)*dVol(k)
            Eeos_expected = Eeos_X + dt*gam*(L_p1 - L_00)
            Eeos_actual = (crad*T4_00 + rho*Cv*T_00)*dVol(k)
            residual = Eeos_actual - Eeos_expected
            rhs(k) = -residual%val
            if (k > 1) sub(k-1) = residual%d1val1
            diag(k) = residual%d1val2
            if (k < nz) sup(k) = residual%d1val3
         end subroutine store_matrix_equation
         
         subroutine solve_matrix_equation(ierr)
            integer, intent(out) :: ierr
            integer :: k
            include 'formats'
            ierr = 0
            sub(nz) = 0d0
            sup(nz) = 0d0
            call solve_tridiag(sub, diag, sup, rhs, deltaT, nz, &
               bp, vp, xp, ierr)
            if (ierr /= 0) then
               write(*,4) 'solve_tridiag failed', newton_iter, stage, model_number
               return
            end if
            do k=1,nz
               T(k) = T(k) + deltaT(k)
            end do            
         end subroutine solve_matrix_equation

      end subroutine calc_implicit


      subroutine solve_tridiag(sub, diag, sup, rhs, x, n, bp, vp, xp, ierr)
         !      sub - sub-diagonal
         !      diag - the main diagonal
         !      sup - sup-diagonal
         !      rhs - right hand side
         !      x - the answer
         !      n - number of equations
         integer, intent(in) :: n
         real(dp), dimension(:), intent(in) :: sup, diag, sub, rhs ! input
         real(dp), dimension(:), intent(out) :: bp, vp, xp ! work
         real(dp), dimension(:), intent(out) :: x ! output
         integer, intent(out) :: ierr
         real(dp) :: m
         integer i
         ierr = 0
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
            d_cons_dt_X, v_bar, & ! output
            ierr)
         real(dp), intent(in) :: T(:), prim(:,:), cons(:,:)
         real(dp), intent(out) :: grad_cell(:,:), flux_face(:,:)
         real(dp), intent(out) :: d_cons_dt_X(:,:), v_bar(:)
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         call get_grad_cell(prim, grad_cell, ierr)
         if (ierr /= 0) return
         call get_flux_face(prim, grad_cell, flux_face, v_bar, ierr)
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
         
         subroutine get1_grad_cell( &
               k, ierr)
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
            flux_face, v_bar, & ! output
            ierr)
         real(dp), intent(in) :: prim(:,:), grad_cell(:,:)
         real(dp), intent(out) :: flux_face(:,:), v_bar(:)
         integer, intent(out) :: ierr
         integer :: j, k, op_err
         include 'formats'
         ierr = 0
         
         !$OMP PARALLEL DO PRIVATE(k,op_err)
         do k=2,nz
            op_err = 0
            call get1_flux_face(k, op_err)
            if (op_err /= 0) ierr = op_err
         end do
         !$OMP END PARALLEL DO
         
         ! ??????????????????
         ! perhaps want area*flux constant across cell 1 for "free flow" outer BC
         do j=1,nvar
            flux_face(j,1) = flux_face(j,2)*area(2)/area(1)
         end do
         
         contains
         
         subroutine get1_flux_face(k, ierr)
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
            PgasL = (gamma-1)*rhoL*Cv*TL
            PgasR = (gamma-1)*rhoR*Cv*TR
            PradL = crad*pow4(TL)/3d0
            PradR = crad*pow4(TR)/3d0
            PL = PgasL + PradL
            PR = PgasR + PradR
            csL = sqrt(Gamma1(PradL/PgasL)*PL/rhoL)
            csR = sqrt(Gamma1(PradR/PgasR)*PR/rhoR)
            fluxL(i_rho) = rhoL*vL
            fluxL(i_mom) = momL*vL + PL
            fluxL(i_etot) = (etotL + PL)*vL
            fluxR(i_rho) = rhoR*vR
            fluxR(i_mom) = momR*vR + PR
            fluxR(i_etot) = (etotR + PR)*vR
            v_bar(k) = (sqrt_rhoL*vL + sqrt_rhoR*vR)/(sqrt_rhoL + sqrt_rhoR)     ! eqn 10.50

            ! Sl and Sr using 1988_einfeldt method
            eta_2 = 0.5*sqrt_rhoL*sqrt_rhoR/pow2(sqrt_rhoL + sqrt_rhoR)        ! eqn 10.54
            d_bar = sqrt((sqrt_rhoL*pow2(csL) + sqrt_rhoR*pow2(csR))/ &
                         (sqrt_rhoL + sqrt_rhoR) + eta_2*pow2(vR - vL))       ! eqn 10.53
            Sl = v_bar(k) - d_bar        ! eqn 10.52
            Sr = v_bar(k) + d_bar        ! eqn 10.52
         
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
               Ss = (PR - PL + vL*rho_v_L - vR*rho_v_R)/(rho_v_L - rho_v_R)  !  10.37
               if (Ss > 0d0) then
                  call get_prim_star(Ss, Sl, vL, rhoL, etotL, PL, prim_star)
                  do j=1,nvar
                     flux_face(j,k) = fluxL(j) + Sl*(prim_star(j) + primL(j))      ! 10.38
                  end do
               else 
                  call get_prim_star(Ss, Sr, vR, rhoR, etotR, PR, prim_star)
                  do j=1,nvar
                     flux_face(j,k) = fluxR(j) + Sr*(prim_star(j) + primR(j))      ! 10.38
                  end do
               end if
            end if
            
         end subroutine get1_flux_face
      
         subroutine get_prim_star(Ss, Sk, vK, rhoK, etotK, PK, prim_star) !     10.37
            real(dp), intent(in) :: Ss, Sk, vK, rhoK, etotK, PK
            real(dp), intent(out) :: prim_star(nvar)
            real(dp) :: prim_star_prefactor
            prim_star_prefactor = rhoK*(Sk - vK)/(Sk - Ss)
            prim_star(i_rho) = prim_star_prefactor
            prim_star(i_mom) = prim_star_prefactor*Ss
            prim_star(i_etot) = prim_star_prefactor* &
               (etotK/rhoK + (Ss - vK)*(Ss + PK/(rhoK*(Sk - vK))))
         end subroutine get_prim_star
      
         real(dp) function Gamma1(z) ! 2020_cheng_shu_song, eqn 2.7
            real(dp) :: z
            Gamma1 = &
               (gamma/(gamma-1d0) + z*(20d0 * 16d0*z))/((1d0/(gamma-1d0) + 12d0*z)*(1d0 + z))
         end function Gamma1
         
      end subroutine get_flux_face
      
      
      real(dp) function get_T(rho, mom, etot) result(T) ! solve for T: etot == rho*Cv*T + 0.5*mom*v + crad*T^4
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
         !$OMP PARALLEL DO PRIVATE(k,op_err)
         do k=1,nz
            op_err = 0
            call get1_d_cons_dt_X(k, op_err)
            if (op_err /= 0) ierr = op_err
         end do
         !$OMP END PARALLEL DO
         
         contains
         
         subroutine get1_d_cons_dt_X(k, ierr)
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            integer :: j
            real(dp) :: Pgas, Prad, Ptot
            include 'formats'
            ierr = 0
            if (k == nz) then
               do j=1,nvar
                  d_cons_dt_X(j,k) = -area(k)*flux_face(j,k)
               end do
               d_cons_dt_X(i_etot,k) = d_cons_dt_X(i_etot,k) + L_inner
            else
               do j=1,nvar
                  d_cons_dt_X(j,k) = area(k+1)*flux_face(j,k+1) - area(k)*flux_face(j,k)
               end do
               Pgas = (gamma-1)*prim(i_rho,k)*Cv*T(k)
               Prad = crad*pow4(T(k))/3
               Ptot = Pgas + Prad
               d_cons_dt_X(i_mom,k) = &
                  d_cons_dt_X(i_mom,k) + (area(k) - area(k+1))*Ptot ! geometry source term
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
      
      
!!! problem setup


      subroutine initialize_problem_A ! Bates 2001
         ! A. Nonlinear Thermal Conduction from a Point Source (the Barenblatt Problem)
         integer :: k
         real(dp) :: E0, rm1, r00, dVol_cntr, R_max, deltaR
         include 'formats'
         nz = 100
         allocate(r(nz), rho(nz), mom(nz), T(nz))
         R_inner = 0d0
         R_max = 1d0
         gamma = 5d0/3d0
         Rgas = 1d0
         Cv = Rgas/(gamma - 1d0)
         ! grid points equally spaced in r
         deltaR = (R_max - R_inner)/nz
         rm1 = R_inner
         do k=nz,1,-1
            r00 = rm1 + deltaR
            r(k) = r00
         end do
         do k=1,nz
            T(k) = 1d-4
            rho(k) = 1d0
            mom(k) = 0d0
         end do
         E0 = 10
         dVol_cntr = 4d0*pi/3d0*(pow3(r(nz)) - pow3(R_inner))
         T(nz) = E0/(rho(nz)*Cv*dVol_cntr)
         dt_next = 1d-6
      end subroutine initialize_problem_A
      

   end module imex
      
