! ***********************************************************************
!
!   Copyright (C) 2013-2019  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************


      module star_solver

      use star_private_def
      use const_def, only: dp
      use num_def
      use mtx_def
      use mtx_lib, only: block_multiply_xa
      use solver_support

      implicit none

      private
      public :: solver, get_solver_work_sizes


      contains


      subroutine solver( &
            s, nvar, skip_global_corr_coeff_limit, &
            gold_tolerances_level, tol_max_correction, tol_correction_norm, &
            work, lwork, iwork, liwork, &
            convergence_failure, ierr)
         use alloc, only: non_crit_get_quad_array, non_crit_return_quad_array
         use utils_lib, only: realloc_if_needed_1, quad_realloc_if_needed_1, fill_with_NaNs

         type (star_info), pointer :: s
         ! the primary variables
         integer, intent(in) :: nvar ! number of variables per zone
         logical, intent(in) :: skip_global_corr_coeff_limit

         ! work arrays. required sizes provided by the routine solver_work_sizes.
         ! for standard use, set work and iwork to 0 before calling.
         ! NOTE: these arrays contain some optional parameter settings and outputs.
         ! see num_def for details.
         integer, intent(in) :: lwork, liwork
         real(dp), intent(inout), target :: work(:) ! (lwork)
         integer, intent(inout), target :: iwork(:) ! (liwork)

         ! convergence criteria
         integer, intent(in) :: gold_tolerances_level ! 0, 1, or 2
         real(dp), intent(in) :: tol_max_correction, tol_correction_norm
            ! a trial solution is considered to have converged if
            ! max_correction <= tol_max_correction and
            !
            ! either
            !          (correction_norm <= tol_correction_norm)
            !    .and. (residual_norm <= tol_residual_norm)
            ! or
            !          (correction_norm*residual_norm <= tol_corr_resid_product)
            !    .and. (abs(slope) <= tol_abs_slope_min)
            !
            ! where "slope" is slope of the line for line search in the solver,
            ! and is analogous to the slope of df/ddx in a 1D solver root finder.

         ! output
         logical, intent(out) :: convergence_failure
         integer, intent(out) :: ierr ! 0 means okay.

         integer :: ldAF, neqns

         include 'formats'

         ierr = 0

         neqns = nvar*s% nz
         ldAF = 3*nvar

         call realloc_if_needed_1(s% AF1,ldAF*neqns,(ldAF+2)*200,ierr)
         if (ierr /= 0) return

         if (s% fill_arrays_with_NaNs) call fill_with_NaNs(s% AF1)

         call do_solver( &
            s, nvar, s% AF1, ldAF, neqns, skip_global_corr_coeff_limit, &
            gold_tolerances_level, tol_max_correction, tol_correction_norm, &
            work, lwork, iwork, liwork, &
            convergence_failure, ierr)

      end subroutine solver


      subroutine do_solver( &
            s, nvar, AF1, ldAF, neq, skip_global_corr_coeff_limit, &
            gold_tolerances_level, tol_max_correction, tol_correction_norm, &
            work, lwork, iwork, liwork, &
            convergence_failure, ierr)

         type (star_info), pointer :: s

         integer, intent(in) :: nvar, ldAF, neq
         logical, intent(in) :: skip_global_corr_coeff_limit

         real(dp), pointer, dimension(:) :: AF1 ! =(ldAF, neq)

         ! controls
         integer, intent(in) :: gold_tolerances_level
         real(dp), intent(in) :: tol_max_correction, tol_correction_norm

         ! work arrays
         integer, intent(in) :: lwork, liwork
         real(dp), intent(inout), target :: work(:) ! (lwork)
         integer, intent(inout), target :: iwork(:) ! (liwork)

         ! output
         logical, intent(out) :: convergence_failure
         integer, intent(out) :: ierr

         ! info saved in work arrays

         real(dp), dimension(:,:), pointer :: &
            dxsave, ddxsave, B, grad_f, soln, &
            row_scale_factors, col_scale_factors, &
            rhs, ddx, xder
         real(dp), dimension(:), pointer :: &
            dxsave1, ddxsave1, B1, grad_f1, soln1, &
            row_scale_factors1, col_scale_factors1, &
            rhs1, ddx1, xder1, &
            save_blks, save_ublk1, save_dblk1, save_lblk1
         integer, dimension(:), pointer :: ipiv1
         character (len=s%nz) :: equed1

         real(dp), dimension(:,:), pointer :: A, Acopy
         real(dp), dimension(:), pointer :: A1, Acopy1
         real(dp), dimension(:), pointer :: lblk1, dblk1, ublk1
         real(dp), dimension(:), pointer :: lblkF1, dblkF1, ublkF1

         ! locals
         real(dp)  ::  &
            coeff, f, slope, residual_norm, max_residual, &
            corr_norm_min, resid_norm_min, correction_factor, temp_correction_factor, &
            correction_norm, corr_norm_initial, max_correction, slope_extra, &
            tol_residual_norm, tol_max_residual, &
            tol_residual_norm2, tol_max_residual2, &
            tol_residual_norm3, tol_max_residual3, &
            tol_abs_slope_min, tol_corr_resid_product, &
            min_corr_coeff, max_corr_min, max_resid_min, max_abs_correction
         integer :: nz, iter, max_tries, zone, tiny_corr_cnt, i, j, k, &
            force_iter_value, iter_for_resid_tol2, iter_for_resid_tol3, &
            max_corr_k, max_corr_j, max_resid_k, max_resid_j
         integer(8) :: test_time1, time0, time1, clock_rate
         character (len=strlen) :: err_msg
         logical :: first_try, dbg_msg, passed_tol_tests, &
            doing_extra, okay, disabled_resid_tests, pass_resid_tests, &
            pass_corr_tests_without_coeff, pass_corr_tests_with_coeff

         integer, parameter :: num_tol_msgs = 15
         character (len=32) :: tol_msg(num_tol_msgs)
         character (len=64) :: message

         real(dp), pointer, dimension(:) :: equ1
         real(dp), pointer, dimension(:,:) :: equ ! (nvar,nz)
         real(dp), pointer, dimension(:,:) :: AF ! (ldAF,neq)
         real(dp), pointer, dimension(:,:,:) :: ublk, dblk, lblk ! (nvar,nvar,nz)
         real(dp), dimension(:,:,:), pointer :: lblkF, dblkF, ublkF ! (nvar,nvar,nz)

         include 'formats'
         
         nz = s% nz

         AF(1:ldAF,1:neq) => AF1(1:ldAF*neq)

         tol_msg(1) = 'avg corr'
         tol_msg(2) = 'max corr '
         tol_msg(3) = 'avg+max corr'
         tol_msg(4) = 'avg resid'
         tol_msg(5) = 'avg corr+resid'
         tol_msg(6) = 'max corr, avg resid'
         tol_msg(7) = 'avg+max corr, avg resid'
         tol_msg(8) = 'max resid'
         tol_msg(9) = 'avg corr, max resid'
         tol_msg(10) = 'max corr+resid'
         tol_msg(11) = 'avg+max corr, max resid'
         tol_msg(12) = 'avg+max resid'
         tol_msg(13) = 'avg corr, avg+max resid'
         tol_msg(14) = 'max corr, avg+max resid'
         tol_msg(15) = 'avg+max corr+resid'

         ierr = 0
         iter = 0
         s% solver_iter = iter

         call set_param_defaults
         dbg_msg = s% report_solver_progress
         
         if (gold_tolerances_level == 2) then
            tol_residual_norm = s% gold2_tol_residual_norm1
            tol_max_residual = s% gold2_tol_max_residual1
            tol_residual_norm2 = s% gold2_tol_residual_norm2
            tol_max_residual2 = s% gold2_tol_max_residual2
            tol_residual_norm3 = s% gold2_tol_residual_norm3
            tol_max_residual3 = s% gold2_tol_max_residual3
         else if (gold_tolerances_level == 1) then
            tol_residual_norm = s% gold_tol_residual_norm1
            tol_max_residual = s% gold_tol_max_residual1
            tol_residual_norm2 = s% gold_tol_residual_norm2
            tol_max_residual2 = s% gold_tol_max_residual2
            tol_residual_norm3 = s% gold_tol_residual_norm3
            tol_max_residual3 = s% gold_tol_max_residual3
         else
            tol_residual_norm = s% tol_residual_norm1
            tol_max_residual = s% tol_max_residual1
            tol_residual_norm2 = s% tol_residual_norm2
            tol_max_residual2 = s% tol_max_residual2
            tol_residual_norm3 = s% tol_residual_norm3
            tol_max_residual3 = s% tol_max_residual3
         end if

         tol_abs_slope_min = -1 ! unused
         tol_corr_resid_product = -1 ! unused
         if (skip_global_corr_coeff_limit) then
            min_corr_coeff = 1
         else
            min_corr_coeff = s% corr_coeff_limit
         end if
         
         if (gold_tolerances_level == 2) then
            iter_for_resid_tol2 = s% gold2_iter_for_resid_tol2
            iter_for_resid_tol3 = s% gold2_iter_for_resid_tol3
         else if (gold_tolerances_level == 1) then
            iter_for_resid_tol2 = s% gold_iter_for_resid_tol2
            iter_for_resid_tol3 = s% gold_iter_for_resid_tol3
         else
            iter_for_resid_tol2 = s% iter_for_resid_tol2
            iter_for_resid_tol3 = s% iter_for_resid_tol3
         end if

         call pointers(ierr)
         if (ierr /= 0) return

         doing_extra = .false.
         passed_tol_tests = .false. ! goes true when pass the tests
         convergence_failure = .false. ! goes true when time to give up
         coeff = 1.d0

         residual_norm=0
         max_residual=0
         corr_norm_min=1d99
         max_corr_min=1d99
         max_resid_min=1d99
         resid_norm_min=1d99
         correction_factor=0
         f=0d0
         slope=0d0

         call set_xscale_info(s, nvar, ierr)
         if (ierr /= 0) then
            if (dbg_msg) &
               write(*, *) 'solver failure: set_xscale_info returned ierr', ierr
            convergence_failure = .true.
            return
         end if
         
         call do_equations(ierr)                 
         if (ierr /= 0) then
            if (dbg_msg) &
               write(*, *) 'solver failure: eval_equations returned ierr', ierr
            convergence_failure = .true.
            return
         end if
         
         call sizequ(s, nvar, residual_norm, max_residual, max_resid_k, max_resid_j, ierr)
         if (ierr /= 0) then
            if (dbg_msg) &
               write(*, *) 'solver failure: sizequ returned ierr', ierr
            convergence_failure = .true.
            return
         end if

         first_try = .true.
         iter = 1
         s% solver_iter = iter
         if (s% doing_first_model_of_run) then
            max_tries = s% max_tries1
         else if (s% retry_cnt > 20) then
            max_tries = s% max_tries_after_20_retries
         else if (s% retry_cnt > 10) then
            max_tries = s% max_tries_after_10_retries
         else if (s% retry_cnt > 5) then
            max_tries = s% max_tries_after_5_retries
         else if (s% retry_cnt > 0) then
            max_tries = s% max_tries_for_retry
         else
            max_tries = s% solver_max_tries_before_reject
         end if
         tiny_corr_cnt = 0

         s% num_solver_iterations = 0
         
      iter_loop: do while (.not. passed_tol_tests)

            if (dbg_msg .and. first_try) write(*, *)
            
            max_resid_j = -1
            max_corr_j = -1

            if (iter >= iter_for_resid_tol2) then
               if (iter < iter_for_resid_tol3) then
                  tol_residual_norm = tol_residual_norm2
                  tol_max_residual = tol_max_residual2
                  if (dbg_msg .and. iter == iter_for_resid_tol2) &
                     write(*,1) 'tol2 residual tolerances: norm, max', &
                        tol_residual_norm, tol_max_residual
               else
                  tol_residual_norm = tol_residual_norm3
                  tol_max_residual = tol_max_residual3
                  if (dbg_msg .and. iter == iter_for_resid_tol3) &
                     write(*,1) 'tol3 residual tolerances: norm, max', &
                        tol_residual_norm, tol_max_residual
               end if
            else if (dbg_msg .and. iter == 1) then
               write(*,2) 'solver_call_number', s% solver_call_number
               write(*,2) 'gold tolerances level', gold_tolerances_level
               write(*,1) 'correction tolerances: norm, max', &
                  tol_correction_norm, tol_max_correction
               write(*,1) 'tol1 residual tolerances: norm, max', &
                  tol_residual_norm, tol_max_residual
            end if
            
            call solver_test_partials(nvar, xder, size(A,dim=1), A1, ierr)
            if (ierr /= 0) then
               call write_msg('solver_test_partials returned ierr /= 0')
               convergence_failure = .true.
               exit iter_loop
            end if
            
            s% num_solver_iterations = s% num_solver_iterations + 1
            if (s% model_number == 1 .and. &
                s% num_solver_iterations > 60 .and. &
                mod(s% num_solver_iterations,10) == 0) &
                  write(*,*) 'first model is slow to converge: num tries', &
                     s% num_solver_iterations
         
            if (.not. solve_equ()) then ! either singular or horribly ill-conditioned
               write(err_msg, '(a, i5, 3x, a)') 'info', ierr, 'bad_matrix'
               call oops(err_msg)
               exit iter_loop
            end if

            call inspectB(s, nvar, soln, ierr)
            if (ierr /= 0) then
               call oops('inspectB returned ierr')
               exit iter_loop
            end if

            ! compute size of scaled correction
            call sizeB(s, nvar, soln, &
                  max_correction, correction_norm, max_corr_k, max_corr_j, ierr)
            if (ierr /= 0) then
               call oops('correction rejected by sizeB')
               exit iter_loop
            end if

            correction_norm = abs(correction_norm)
            max_abs_correction = abs(max_correction)
            corr_norm_min = min(correction_norm, corr_norm_min)
            max_corr_min = min(max_abs_correction, max_corr_min)

            if (is_bad_num(correction_norm) .or. is_bad_num(max_abs_correction)) then
               ! bad news -- bogus correction
               call oops('bad result from sizeB -- correction info either NaN or Inf')
               if (s% stop_for_bad_nums) then
                  write(*,1) 'correction_norm', correction_norm
                  write(*,1) 'max_correction', max_correction
                  stop 'solver'
               end if
               exit iter_loop
            end if

            if (.not. s% ignore_too_large_correction) then
               if ((correction_norm > s% corr_param_factor*s% scale_correction_norm) .and. &
                     .not. s% doing_first_model_of_run) then
                  call oops('avg corr too large')
                  exit iter_loop
               endif
            end if

            ! shrink the correction if it is too large
            correction_factor = 1d0
            temp_correction_factor = 1d0

            if (correction_norm*correction_factor > s% scale_correction_norm) then
               correction_factor = min(correction_factor,s% scale_correction_norm/correction_norm)
            end if
            
            if (max_abs_correction*correction_factor > s% scale_max_correction) then
               temp_correction_factor = s% scale_max_correction/max_abs_correction
            end if

            if (iter > s% solver_itermin_until_reduce_min_corr_coeff) then
               if (min_corr_coeff == 1d0 .and. &
                  s% solver_reduced_min_corr_coeff < 1d0) then
                     min_corr_coeff = s% solver_reduced_min_corr_coeff
               end if
            end if

            correction_factor = max(min_corr_coeff, correction_factor)
            if (.not. s% ignore_min_corr_coeff_for_scale_max_correction) then
               temp_correction_factor = max(min_corr_coeff, temp_correction_factor)
            end if
            correction_factor = min(correction_factor, temp_correction_factor)

            ! fix B if out of definition domain
            call Bdomain(s, nvar, soln, correction_factor, ierr)
            if (ierr /= 0) then ! correction cannot be fixed
               call oops('correction rejected by Bdomain')
               exit iter_loop
            end if

            if (min_corr_coeff < 1d0) then
               ! compute gradient of f = equ<dot>jacobian
               ! NOTE: NOT jacobian<dot>equ
               call block_multiply_xa(nvar, nz, lblk1, dblk1, ublk1, equ1, grad_f1)

               slope = eval_slope(nvar, nz, grad_f, soln)
               if (is_bad_num(slope) .or. slope > 0d0) then ! a very bad sign
                  if (is_bad_num(slope) .and. s% stop_for_bad_nums) then
                     write(*,1) 'slope', slope
                     stop 'solver'
                  end if
                  slope = 0d0
                  min_corr_coeff = 1d0
               end if

            else

               slope = 0d0

            end if
            
            f = 0d0
            call adjust_correction( &
               min_corr_coeff, correction_factor, grad_f1, f, slope, coeff, err_msg, ierr)
            if (ierr /= 0) then
               call oops(err_msg)
               exit iter_loop
            end if
            s% solver_adjust_iter = 0

            ! coeff is factor by which adjust_correction rescaled the correction vector
            if (coeff > s% tiny_corr_factor*min_corr_coeff .or. min_corr_coeff >= 1d0) then
               tiny_corr_cnt = 0
            else
               tiny_corr_cnt = tiny_corr_cnt + 1
            end if

            ! check the residuals for the equations

            call sizequ(s, nvar, residual_norm, max_residual, max_resid_k, max_resid_j, ierr)
            if (ierr /= 0) then
               call oops('sizequ returned ierr')
               exit iter_loop
            end if

            if (is_bad_num(residual_norm)) then
               call oops('residual_norm is a a bad number (NaN or Infinity)')
               if (s% stop_for_bad_nums) then
                  write(*,1) 'residual_norm', residual_norm
                  stop 'solver'
               end if
               exit iter_loop
            end if
            
            if (is_bad_num(max_residual)) then
               call oops('max_residual is a a bad number (NaN or Infinity)')
               if (s% stop_for_bad_nums) then
                  write(*,1) 'max_residual', max_residual
                  stop 'solver'
               end if
               exit iter_loop
            end if

            residual_norm = abs(residual_norm)
            max_residual = abs(max_residual)
            s% residual_norm = residual_norm
            s% max_residual = max_residual
            resid_norm_min = min(residual_norm, resid_norm_min)
            max_resid_min = min(max_residual, max_resid_min)
            
            disabled_resid_tests = &
               tol_max_residual > 1d2 .and. tol_residual_norm > 1d2
            pass_resid_tests = &
               .not. disabled_resid_tests .and. &
               max_residual <= tol_max_residual .and. &
               residual_norm <= tol_residual_norm
            pass_corr_tests_without_coeff = &
               max_abs_correction <= tol_max_correction .and. &
               correction_norm <= tol_correction_norm
            pass_corr_tests_with_coeff = &
               max_abs_correction <= tol_max_correction*coeff .and. &
               correction_norm <= tol_correction_norm*coeff

            passed_tol_tests = &
               (pass_resid_tests .and. pass_corr_tests_with_coeff) .or. &
               (disabled_resid_tests .and. pass_corr_tests_without_coeff)
            
            if (.not. passed_tol_tests) then

               if (iter >= max_tries) then
                  call get_message
                  message = trim(message) // ' -- give up'
                  if (len_trim(s% retry_message) == 0) &
                     s% retry_message = trim(message) // ' in solver'
                  if (dbg_msg) call write_msg(message)
                  if (.not. pass_resid_tests .and. .not. disabled_resid_tests) then
                     if (residual_norm > tol_residual_norm) &
                        write(*,2) 'residual_norm > tol_residual_norm', &
                           s% model_number, residual_norm, tol_residual_norm
                     if (max_residual > tol_max_residual) &
                        write(*,2) 'max_residual > tol_max_residual', &
                           s% model_number, max_residual, tol_max_residual
                  end if
                  if (disabled_resid_tests) then ! no coeff for corrections
                     if (correction_norm > tol_correction_norm) &
                        write(*,2) 'correction_norm > tol_correction_norm', &
                           s% model_number, correction_norm, tol_correction_norm, coeff
                     if (max_abs_correction > tol_max_correction) &
                        write(*,2) 'max_abs_correction > tol_max_correction', &
                           s% model_number, max_abs_correction, tol_max_correction, coeff
                  else ! include coeff for corrections
                     if (correction_norm > tol_correction_norm*coeff) &
                        write(*,2) 'correction_norm > tol_correction_norm*coeff', &
                           s% model_number, correction_norm, tol_correction_norm*coeff, coeff
                     if (max_abs_correction > tol_max_correction*coeff) &
                        write(*,2) 'max_abs_correction > tol_max_correction*coeff', &
                           s% model_number, max_abs_correction, tol_max_correction*coeff, coeff
                  end if
                  convergence_failure = .true.
                  exit iter_loop
               else if (.not. first_try .and. .not. s% doing_first_model_of_run) then
                  if (correction_norm > s% corr_norm_jump_limit*corr_norm_min) then
                     call oops('avg correction jumped')
                     exit iter_loop
                  else if (residual_norm > s% resid_norm_jump_limit*resid_norm_min) then
                     call oops('avg residual jumped')
                     exit iter_loop
                  else if (max_abs_correction > s% max_corr_jump_limit*max_corr_min) then
                     call oops('max correction jumped')
                     exit iter_loop
                  else if (max_residual > s% max_resid_jump_limit*max_resid_min) then
                     call oops('max residual jumped')
                     exit iter_loop
                  else if (tiny_corr_cnt >= s% tiny_corr_coeff_limit &
                        .and. min_corr_coeff < 1) then
                     call oops('tiny corrections')
                     exit iter_loop
                  end if
               else if (.not. s% doing_first_model_of_run) then
                  if (coeff < min(min_corr_coeff,correction_factor)) then
                     call oops('coeff too small')
                     exit iter_loop
                  end if
               end if
            end if

            if (dbg_msg) then
               if (.not. passed_tol_tests) then
                  call get_message
               end if
               if (.not. passed_tol_tests) then
                  call write_msg(message)
               else if (iter < s% solver_itermin) then
                  call write_msg('iter < itermin')
               else
                  call write_msg('okay!')
               end if
            end if

            if (passed_tol_tests .and. (iter+1 < max_tries)) then
               ! about to declare victory... but may want to do another iteration
               force_iter_value = force_another_iteration(s, iter, s% solver_itermin)
               if (force_iter_value > 0) then
                  passed_tol_tests = .false. ! force another
                  tiny_corr_cnt = 0 ! reset the counter
                  corr_norm_min = 1d99
                  resid_norm_min = 1d99
                  max_corr_min = 1d99
                  max_resid_min = 1d99
               else if (force_iter_value < 0) then ! failure
                  call oops('force iter')
                  exit iter_loop
               end if
            end if

            if (s% use_other_solver_monitor .and. &
                  associated(s% other_solver_monitor)) then
               call s% other_solver_monitor( &
                  s% id, iter, passed_tol_tests, &
                  correction_norm, max_correction, &
                  residual_norm, max_residual, ierr)
               if (ierr /= 0) then
                  call oops('other_solver_monitor')
                  exit iter_loop
               end if
            end if

            iter=iter+1
            s% solver_iter = iter
            first_try = .false.

         end do iter_loop
            
         if (max_residual > s% warning_limit_for_max_residual .and. .not. convergence_failure) &
            write(*,2) 'WARNING: max_residual > warning_limit_for_max_residual', &
               s% model_number, max_residual, s% warning_limit_for_max_residual


         contains
         
         
         subroutine solver_test_partials(nvar, xder, ldA, A1, ierr)
            ! create jacobian by using numerical differences for partial derivatives
            integer, intent(in) :: nvar
            real(dp), pointer, dimension(:,:) :: xder ! (nvar, nz)
            integer, intent(in) :: ldA ! leading dimension of A
            real(dp), pointer, dimension(:) :: A1
            integer, intent(out) :: ierr
            
            integer :: j, k, i_var, i_var_sink, i_equ, k_off, cnt_00, cnt_m1, cnt_p1, k_lo, k_hi
            real(dp), dimension(:,:), pointer :: save_equ, save_dx
            real(dp) :: dvar, dequ, dxtra, &
               dx_0, dvardx, dvardx_0, xdum, err
            logical :: testing_partial

            include 'formats'

            ierr = 0
            testing_partial = & ! check inlist parameters
               s% solver_test_partials_dx_0 > 0d0 .and. &
               s% solver_test_partials_k > 0 .and. &
               s% solver_call_number == s% solver_test_partials_call_number .and. &
               s% solver_test_partials_iter_number == iter
            if (.not. testing_partial) return

            call do_equations(ierr)
            if (ierr /= 0) return

            allocate(save_dx(nvar,nz), save_equ(nvar,nz))

            do k=1,nz
               do j=1,nvar
                  save_dx(j,k) = s% solver_dx(j,k)
                  save_equ(j,k) = equ(j,k)
               end do
            end do
            
            s% doing_check_partials = .true. ! let set_vars_for_solver know
            k_lo = s% solver_test_partials_k_low
            if (k_lo > 0 .and. k_lo <= s% nz) then
               k_hi = s% solver_test_partials_k_high
               if (k_hi <= 0) then
                  k_hi = s% nz
               else
                  k_hi = min(k_hi,s% nz)
               end if
               do k = k_lo, k_hi
                  call test_cell_partials(s, k, save_dx, save_equ, ierr)
                  if (ierr /= 0) stop 'failed solver_test_partials'
               end do
            else
               k = s% solver_test_partials_k
               call test_cell_partials(s, k, save_dx, save_equ, ierr) 
               if (ierr /= 0) stop 'failed solver_test_partials'
            end if
            deallocate(save_dx, save_equ)
            stop 'done solver_test_partials'

         end subroutine solver_test_partials


         subroutine get_message
            include 'formats'
            i = 0
            if (correction_norm > tol_correction_norm*coeff) i = i+1
            if (max_abs_correction > tol_max_correction*coeff) i = i+2
            if (residual_norm > tol_residual_norm*coeff) i = i+4
            if (max_residual > tol_max_residual*coeff) i = i+8
            if (i == 0) then
               message = 'out of tries'
            else
               message = tol_msg(i)
            end if
         end subroutine get_message


         subroutine set_param_defaults
            if (s% corr_param_factor == 0) s% corr_param_factor = 10d0
            if (s% scale_max_correction == 0) s% scale_max_correction = 1d99
            if (s% corr_norm_jump_limit == 0) s% corr_norm_jump_limit = 1d99
            if (s% max_corr_jump_limit == 0) s% max_corr_jump_limit = 1d99
         end subroutine set_param_defaults


         subroutine oops(msg)
            character (len=*), intent(in) :: msg
            character (len=strlen) :: full_msg
            include 'formats'
            full_msg = trim(msg) // ' -- give up'
            if (len_trim(s% retry_message) == 0) s% retry_message = trim(full_msg) // ' in solver'
            call write_msg(full_msg)
            convergence_failure = .true.
         end subroutine oops
         
         
         subroutine do_equations(ierr)
            integer, intent(out) :: ierr
            call prepare_solver_matrix(nvar, xder, size(A,dim=1), A1, ierr)
            if (ierr /= 0) return
            call eval_equations(s, nvar, ierr)
            if (ierr /= 0) return
            call s% other_after_solver_setmatrix(s% id, ierr)
         end subroutine do_equations


         subroutine prepare_solver_matrix(nvar, xder, ldA, A1, ierr)
            integer, intent(in) :: nvar
            real(dp), pointer, dimension(:,:) :: xder ! (nvar, nz)
            integer, intent(in) :: ldA ! leading dimension of A
            real(dp), pointer, dimension(:) :: A1
            integer, intent(out) :: ierr
            real(dp), pointer, dimension(:,:) :: A ! (ldA, neqns)
            integer :: i, j, nz, neqns
            include 'formats'
            ierr = 0
            nz = s% nz
            neqns = nvar*nz
            A(1:ldA,1:neqns) => A1(1:ldA*neqns)         
            i = nvar*nvar*nz
            if (size(A1,dim=1) < 3*i) then
               write(*,*) 'prepare_solver_matrix: size(A1,dim=1) < 3*i', size(A1,dim=1), 3*i
               ierr = -1
               return
            end if
            s% ublk(1:nvar,1:nvar,1:nz) => A1(1:i)
            s% dblk(1:nvar,1:nvar,1:nz) => A1(i+1:2*i)
            s% lblk(1:nvar,1:nvar,1:nz) => A1(2*i+1:3*i)
         end subroutine prepare_solver_matrix


         subroutine adjust_correction( &
               min_corr_coeff_in, max_corr_coeff, grad_f, f, slope, coeff,  &
               err_msg, ierr)
            real(dp), intent(in) :: min_corr_coeff_in
            real(dp), intent(in) :: max_corr_coeff
            real(dp), intent(in) :: grad_f(:) ! (neq) ! gradient df/ddx at xold
            real(dp), intent(out) :: f ! 1/2 fvec^2. minimize this.
            real(dp), intent(in) :: slope
            real(dp), intent(out) :: coeff

            ! the new correction is coeff*xscale*soln
            ! with min_corr_coeff <= coeff <= max_corr_coeff
            ! if all goes well, the new x will give an improvement in f

            character (len=*), intent(out) :: err_msg
            integer, intent(out) :: ierr

            integer :: i, j, k, iter, k_max_corr, i_max_corr
            character (len=strlen) :: message
            logical :: first_time
            real(dp) :: a1, alam, alam2, alamin, a2, disc, f2, &
               rhs1, rhs2, temp, test, tmplam, max_corr, fold, min_corr_coeff
            real(dp) :: frac, f_target
            logical :: skip_eval_f, dbg_adjust

            real(dp), parameter :: alf = 1d-2 ! ensures sufficient decrease in f

            real(dp), parameter :: alam_factor = 0.2d0

            include 'formats'

            ierr = 0
            coeff = 0
            dbg_adjust = .false.

            skip_eval_f = (min_corr_coeff_in == 1)
            if (skip_eval_f) then
               f = 0
            else
               do k=1,nz
                  do i=1,nvar
                     dxsave(i,k) = s% solver_dx(i,k)
                     ddxsave(i,k) = ddx(i,k)
                  end do
               end do
               f = eval_f(nvar,nz,equ)
               if (is_bad_num(f)) then
                  ierr = -1
                  write(err_msg,*) 'adjust_correction failed in eval_f'
                  if (dbg_msg) write(*,*) &
                     'adjust_correction: eval_f(nvar,nz,equ)', eval_f(nvar,nz,equ)
                  if (s% stop_for_bad_nums) then
                     write(*,1) 'f', f
                     stop 'solver adjust_correction'
                  end if
                  return
               end if
            end if
            fold = f
            min_corr_coeff = min(min_corr_coeff_in, max_corr_coeff) ! make sure min <= max
            alam = max_corr_coeff
            first_time = .true.
            f2 = 0
            alam2 = 0
            if (dbg_adjust) then
               write(*,4) 'max_corr_coeff', k, s% solver_iter, &
                  s% model_number, max_corr_coeff
               write(*,4) 'slope', k, s% solver_iter, &
                  s% model_number, slope
               write(*,4) 'f', k, s% solver_iter, &
                  s% model_number, f
            end if

         search_loop: do iter = 1, 1000

               coeff = max(min_corr_coeff, alam)
               s% solver_adjust_iter = iter

               call apply_coeff(nvar, nz, dxsave, soln, coeff, skip_eval_f)
               
               call do_equations(ierr)               
               if (ierr /= 0) then
                  if (alam > min_corr_coeff .and. s% model_number == 1) then
                     ! try again with smaller correction vector.
                     ! need this to rescue create pre-main-sequence model in some nasty cases.
                     alam = max(alam/10, min_corr_coeff)
                     ierr = 0
                     cycle
                  end if
                  write(err_msg,'(a)') 'adjust_correction failed in eval_equations'
                  if (dbg_msg .or. dbg_adjust) &
                     write(*,2) 'adjust_correction: eval_equations returned ierr', &
                        ierr, min_corr_coeff, max_corr_coeff
                  exit search_loop
               end if

               if (min_corr_coeff == 1) return

               if (dbg_adjust) then
                  do k=1,nz
                     do i=1,nvar
                        write(*,5) trim(s% nameofequ(i)), k, iter, s% solver_iter, &
                           s% model_number, equ(i,k)
                     end do
                  end do
               end if

               f = eval_f(nvar,nz,equ)
               if (is_bad_num(f)) then
                  if (s% stop_for_bad_nums) then
                     write(*,1) 'f', f
                     stop 'solver adjust_correction eval_f'
                  end if
                  if (alam > min_corr_coeff) then
                     alam = max(alam/10, min_corr_coeff)
                     ierr = 0
                     cycle
                  end if
                  err_msg = 'equ norm is NaN or other bad num'
                  ierr = -1
                  exit search_loop
               end if

               f_target = max(fold/2, fold + alf*coeff*slope)
               if (f <= f_target) then
                  return ! sufficient decrease in f
               end if

               if (alam <= min_corr_coeff) then
                  return ! time to give up
               end if

               ! reduce alam and try again
               if (first_time) then
                  tmplam = -slope/(2*(f-fold-slope))
                  first_time = .false.
                  if (dbg_adjust) then
                     write(*,5) 'slope', k, iter, s% solver_iter, &
                        s% model_number, slope
                     write(*,5) 'f', k, iter, s% solver_iter, &
                        s% model_number, f
                     write(*,5) 'fold', k, iter, s% solver_iter, &
                        s% model_number, fold
                     write(*,5) '2*(f-fold-slope)', k, iter, s% solver_iter, &
                        s% model_number, 2*(f-fold-slope)
                  end if
               else ! have two prior f values to work with
                  rhs1 = f - fold - alam*slope
                  rhs2 = f2 - fold - alam2*slope
                  a1 = (rhs1/(alam*alam) - rhs2/(alam2*alam2))/(alam - alam2)
                  a2 = (-alam2*rhs1/(alam*alam) + alam*rhs2/(alam2*alam2))/(alam - alam2)
                  if (dbg_adjust) then
                     write(*,5) 'slope', k, iter, s% solver_iter, &
                        s% model_number, slope
                     write(*,5) 'f', k, iter, s% solver_iter, &
                        s% model_number, f
                     write(*,5) 'f2', k, iter, s% solver_iter, &
                        s% model_number, f2
                     write(*,5) 'fold', k, iter, s% solver_iter, &
                        s% model_number, fold
                     write(*,5) 'alam', k, iter, s% solver_iter, &
                        s% model_number, alam
                     write(*,5) 'alam2', k, iter, s% solver_iter, &
                        s% model_number, alam2
                     write(*,5) 'rhs1', k, iter, s% solver_iter, &
                        s% model_number, rhs1
                     write(*,5) 'rhs2', k, iter, s% solver_iter, &
                        s% model_number, rhs2
                     write(*,5) 'a1', k, iter, s% solver_iter, &
                        s% model_number, a1
                     write(*,5) 'a2', k, iter, s% solver_iter, &
                        s% model_number, a2
                  end if
                  if (a1 == 0) then
                     tmplam = -slope/(2*a2)
                  else
                     disc = a2*a2-3*a1*slope
                     if (disc < 0) then
                        tmplam = alam*alam_factor
                     else if (a2 <= 0) then
                        tmplam = (-a2+sqrt(disc))/(3*a1)
                     else
                        tmplam = -slope/(a2+sqrt(disc))
                     end if
                     if (dbg_adjust) then
                        write(*,5) 'disc', k, iter, s% solver_iter, &
                           s% model_number, disc
                     end if
                  end if
                  if (tmplam > alam*alam_factor) tmplam = alam*alam_factor
               end if

               alam2 = alam
               f2 = f
               alam = max(tmplam, alam*alam_factor, min_corr_coeff)

               if (dbg_adjust) then
                  write(*,5) 'tmplam', k, iter, s% solver_iter, &
                     s% model_number, tmplam
                  write(*,5) 'min_corr_coeff', k, iter, s% solver_iter, &
                     s% model_number, min_corr_coeff
                  write(*,5) 'alam_factor', k, iter, s% solver_iter, &
                     s% model_number, alam_factor
               end if

            end do search_loop

            do k=1,nz
               do i=1,nvar
                  s% solver_dx(i,k) = dxsave(i,k)
                  ddx(i,k) = ddxsave(i,k)
               end do
            end do

         end subroutine adjust_correction


         subroutine apply_coeff(nvar, nz, dxsave, soln, coeff, just_use_dx)
            integer, intent(in) :: nvar, nz
            real(dp), intent(in), dimension(:,:) :: dxsave, soln
            real(dp), intent(in) :: coeff
            logical, intent(in) :: just_use_dx
            integer :: i, k
            include 'formats'

            if (just_use_dx) then
               if (coeff == 1d0) then
                  do k=1,nz
                     do i=1,nvar
                        s% solver_dx(i,k) = s% solver_dx(i,k) + s% x_scale(i,k)*soln(i,k)
                     end do
                  end do
               else
                  do k=1,nz
                     do i=1,nvar
                        s% solver_dx(i,k) = s% solver_dx(i,k) + coeff*s% x_scale(i,k)*soln(i,k)
                     end do
                  end do
               end if
               return
            end if
            ! else use dxsave instead of dx
            if (coeff == 1d0) then
               do k=1,nz
                  do i=1,nvar
                     s% solver_dx(i,k) = dxsave(i,k) + s% x_scale(i,k)*soln(i,k)
                  end do
               end do
               return
            end if
            do k=1,nz
               do i=1,nvar
                  s% solver_dx(i,k) = dxsave(i,k) + coeff*s% x_scale(i,k)*soln(i,k)
               end do
            end do
         end subroutine apply_coeff


         logical function solve_equ()
            use star_utils, only: start_time, update_time
            integer ::  i, k, ierr
            real(dp) :: ferr, berr, total_time
            logical :: done
            include 'formats'
            ierr = 0
            solve_equ = .true.

            if (s% doing_timing) then
               call start_time(s, time0, total_time)
            end if
         
            !$omp simd
            do i=1,neq
               b1(i) = -equ1(i) ! b1 is rhs of matrix equation
            end do
            
            done = .false.
            if (s% x_logical_ctrl(19)) then ! testing GMRES
            
               ! preconditioning
               ! The special case of block tridiagonality
               ! http://www.netlib.org/linalg/html_templates/node72.html
            
               soln1(1:neq) = 0d0 ! simplest initial guess for now
               call solve_mtx_eqn_with_GMRES( &
                  s, nvar, nz, b1, soln1, lblk, dblk, ublk, &
                  lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, &
                  row_scale_factors1, col_scale_factors1, &
                  3*nvar*neq, save_blks, ipiv1, ierr)
            
               if (ierr /= 0) then
                  ierr = 0
               else
                  done = .true.
               end if
            
            end if
            
            if (.not. done) then
            
               if (s% use_DGESVX_in_bcyclic) then
                  !$omp simd
                  do i = 1, nvar*nvar*nz
                     save_ublk1(i) = ublk1(i)
                     save_dblk1(i) = dblk1(i)
                     save_lblk1(i) = lblk1(i)
                  end do
               end if
            
               call factor_mtx(ierr)
               if (ierr == 0) call solve_mtx(ierr)
            
               if (s% use_DGESVX_in_bcyclic) then
                  !$omp simd
                  do i = 1, nvar*nvar*nz
                     ublk1(i) = save_ublk1(i)
                     dblk1(i) = save_dblk1(i)
                     lblk1(i) = save_lblk1(i)
                  end do
               end if
            
            end if

            if (s% doing_timing) then
               call update_time(s, time0, total_time, s% time_solver_matrix)
            end if
            if (ierr /= 0) then
               solve_equ = .false.
               b(1:nvar,1:nz) = 0d0
            end if

         end function solve_equ


         subroutine factor_mtx(ierr)
            use star_bcyclic, only: bcyclic_factor
            integer, intent(out) :: ierr
            call bcyclic_factor( &
               s, nvar, nz, lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipiv1, &
               B1, row_scale_factors1, col_scale_factors1, &
               equed1, iter, ierr)
         end subroutine factor_mtx


         subroutine solve_mtx(ierr)
            use star_bcyclic, only: bcyclic_solve
            integer, intent(out) :: ierr
            call bcyclic_solve( &
               s, nvar, nz, lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipiv1, &
               B1, soln1, row_scale_factors1, col_scale_factors1, equed1, &
               iter, ierr)
         end subroutine solve_mtx


         subroutine test_cell_partials(s, k, save_dx, save_equ, ierr) 
            use star_utils, only: lookup_nameofvar, lookup_nameofequ
            type (star_info), pointer :: s
            integer, intent(in) :: k
            real(dp), pointer, dimension(:,:) :: save_dx, save_equ
            integer, intent(out) :: ierr
            integer :: i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index
            include 'formats'
            ierr = 0
            write(*,*)
            i_equ = lookup_nameofequ(s, s% solver_test_partials_equ_name)      
            if (i_equ == 0 .and. len_trim(s% solver_test_partials_equ_name) > 0) then
               if (s% solver_test_partials_equ_name == 'lnE') then ! testing eos
                  i_equ = -1
               else if (s% solver_test_partials_equ_name == 'eps_nuc') then
                  i_equ = -2
               else if (s% solver_test_partials_equ_name == 'opacity') then
                  i_equ = -3
               else if (s% solver_test_partials_equ_name == 'lnP') then
                  i_equ = -4
               else if (s% solver_test_partials_equ_name == 'non_nuc_neu') then
                  i_equ = -5
               else if (s% solver_test_partials_equ_name == 'gradT') then
                  i_equ = -6
               else if (s% solver_test_partials_equ_name == 'mlt_vc') then
                  i_equ = -7
               else if (s% solver_test_partials_equ_name == 'grad_ad') then
                  i_equ = -8
               end if 
            else if (i_equ /= 0) then
               write(*,1) 'equ name ' // trim(s% solver_test_partials_equ_name)
            end if
            i_var = lookup_nameofvar(s, s% solver_test_partials_var_name)            
            if (i_var /= 0) write(*,1) 'var name ' // trim(s% solver_test_partials_var_name)
            if (i_var > s% nvar_hydro) then ! get index in xa
               i_var_xa_index = i_var - s% nvar_hydro
            else
               i_var_xa_index = 0
            end if
            i_var_sink = lookup_nameofvar(s, s% solver_test_partials_sink_name)
            i_var_sink_xa_index  = 0
            if (i_var_sink > 0 .and. i_var > s% nvar_hydro) then
               write(*,1) 'sink name ' // trim(s% solver_test_partials_sink_name)
               if (i_var_sink > s% nvar_hydro) then ! get index in xa
                  i_var_sink_xa_index = i_var_sink - s% nvar_hydro
               else
                  write(*,*) 'ERROR: sink name must be a chem name for the current net'
                  ierr = -1
                  return
               end if
            end if
            if (s% solver_test_partials_equ_name == 'all') then
               do i_equ = 1, s% nvar_hydro
                  call test_equ_partials(s, &
                     i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                     k, save_dx, save_equ, ierr)   
                  if (ierr /= 0) stop 'failed solver_test_partials'
               end do
            else
               call test_equ_partials(s, &
                  i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  k, save_dx, save_equ, ierr)   
               if (ierr /= 0) stop 'failed solver_test_partials'
            end if     
         end subroutine test_cell_partials               


         subroutine test_equ_partials(s, &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               k, save_dx, save_equ, ierr)
            type (star_info), pointer :: s
            integer, intent(in) :: &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, k
            real(dp), pointer, dimension(:,:) :: save_dx, save_equ
            integer, intent(out) :: ierr
            real(dp) :: dvardx_0
            integer :: i, j_var_xa_index, j_var_sink_xa_index
            include 'formats'
            if (i_equ /= 0) then
               if (s% solver_test_partials_var_name == 'all') then
                  do i = 1, s% nvar_hydro
                     call test3_partials(s, &
                        i_equ, i, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                        k, save_dx, save_equ, ierr)
                     if (ierr /= 0) stop 'failed solver_test_partials'
                     write(*,*)
                  end do
               else if (i_var == 0) then
                  write(*,*) 'failed to recognize variable name'
               else
                  call test3_partials(s, &
                     i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                     k, save_dx, save_equ, ierr)
                  if (ierr /= 0) stop 'failed solver_test_partials'               
               end if
            else ! i_equ == 0
               if (i_var /= 0) then
                  call test1_partial(s, &
                     i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                     k, 0, s% solver_test_partials_dval_dx, save_dx, save_equ, ierr)
               else ! i_var == 0
                  if (s% solver_test_partials_var <= 0) then
                     write(*,2) 'need to set solver_test_partials_var', s% solver_test_partials_var
                     write(*,2) 'for solver_test_partials_k', s% solver_test_partials_k
                     stop 'failed solver_test_partials'
                  end if
                  if (s% solver_test_partials_var > s% nvar_hydro) then
                     j_var_xa_index = s% solver_test_partials_var - s% nvar_hydro
                     if (s% solver_test_partials_dx_sink > s% nvar_hydro) then
                        j_var_sink_xa_index = s% solver_test_partials_dx_sink - s% nvar_hydro
                     else
                        write(*,*) 'set solver_test_partials_dx_sink to variable index, not to xa index', &
                           s% solver_test_partials_dx_sink
                        stop 'failed solver_test_partials'
                     end if
                  else
                     j_var_xa_index = 0
                     j_var_sink_xa_index = 0
                  end if
                  call test1_partial(s, &
                     i_equ, s% solver_test_partials_var, s% solver_test_partials_dx_sink, &
                     j_var_xa_index, j_var_sink_xa_index, &                     
                     k, 0, s% solver_test_partials_dval_dx, save_dx, save_equ, ierr)
               end if
               if (ierr /= 0) stop 'failed solver_test_partials'
            end if               
            write(*,*)
         end subroutine test_equ_partials
         
         
         subroutine get_lnE_partials(s, &
               k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               dvardx0_m1, dvardx0_00, dvardx0_p1)
            use eos_def, only: i_lnE
            type (star_info), pointer :: s
            integer, intent(in) :: k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index
            real(dp), intent(out) :: dvardx0_m1, dvardx0_00, dvardx0_p1
            dvardx0_m1 = 0d0; dvardx0_00 = 0d0; dvardx0_p1 = 0d0
            if (i_var_xa_index > 0) then 
               dvardx0_00 = s% d_eos_dxa(i_lnE,i_var_xa_index,k) - &
                        s% d_eos_dxa(i_lnE,i_var_sink_xa_index,k)
            else if (i_var == s% i_lnd) then
               dvardx0_00 = s% dE_dRho_for_partials(k)*s% rho(k)/s% energy(k)
            else if (i_var == s% i_lnT) then
               dvardx0_00 = s% Cv_for_partials(k)*s% T(k)/s% energy(k)
            end if
         end subroutine get_lnE_partials
         
         
         subroutine get_lnP_partials(s, &
               k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               dvardx0_m1, dvardx0_00, dvardx0_p1)
            use eos_def, only: i_lnPgas
            type (star_info), pointer :: s
            integer, intent(in) :: k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index
            real(dp), intent(out) :: dvardx0_m1, dvardx0_00, dvardx0_p1
            dvardx0_m1 = 0d0; dvardx0_00 = 0d0; dvardx0_p1 = 0d0
            if (i_var_xa_index > 0) then 
               dvardx0_00 = s% Pgas(k)/s% Peos(k) * &
                  (s% d_eos_dxa(i_lnPgas,i_var_xa_index,k) - s% d_eos_dxa(i_lnPgas,i_var_sink_xa_index,k))
            else if (i_var == s% i_lnd) then
               dvardx0_00 = s% chiRho_for_partials(k)
            else if (i_var == s% i_lnT) then
               dvardx0_00 = s% chiT_for_partials(k)
            end if
         end subroutine get_lnP_partials
         
         
         subroutine get_grad_ad_partials(s, &
               k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               dvardx0_m1, dvardx0_00, dvardx0_p1)
            use eos_def, only: i_grad_ad
            type (star_info), pointer :: s
            integer, intent(in) :: k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index
            real(dp), intent(out) :: dvardx0_m1, dvardx0_00, dvardx0_p1
            dvardx0_m1 = 0d0; dvardx0_00 = 0d0; dvardx0_p1 = 0d0
            if (i_var_xa_index > 0) then 
               dvardx0_00 = &
                  (s% d_eos_dxa(i_grad_ad,i_var_xa_index,k) - s% d_eos_dxa(i_grad_ad,i_var_sink_xa_index,k))
            else if (i_var == s% i_lnd) then
               dvardx0_00 = s% d_eos_dlnd(i_grad_ad,k)
            else if (i_var == s% i_lnT) then
               dvardx0_00 = s% d_eos_dlnT(i_grad_ad,k)
            end if
         end subroutine get_grad_ad_partials
         
         
         subroutine get_eps_nuc_partials(s, &
               k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               dvardx0_m1, dvardx0_00, dvardx0_p1)
            type (star_info), pointer :: s
            integer, intent(in) :: k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index
            real(dp), intent(out) :: dvardx0_m1, dvardx0_00, dvardx0_p1
            dvardx0_m1 = 0d0; dvardx0_00 = 0d0; dvardx0_p1 = 0d0
            if (i_var > s% nvar_hydro) then 
               dvardx0_00 = s% d_epsnuc_dx(i_var_xa_index,k) - s% d_epsnuc_dx(i_var_sink_xa_index,k)
            else if (i_var == s% i_lnd) then
               dvardx0_00 = s% d_epsnuc_dlnd(k)
            else if (i_var == s% i_lnT) then
               dvardx0_00 = s% d_epsnuc_dlnT(k)
            end if
         end subroutine get_eps_nuc_partials
         
         
         subroutine get_non_nuc_neu_partials(s, &
               k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               dvardx0_m1, dvardx0_00, dvardx0_p1)
            type (star_info), pointer :: s
            integer, intent(in) :: k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index
            real(dp), intent(out) :: dvardx0_m1, dvardx0_00, dvardx0_p1
            dvardx0_m1 = 0d0; dvardx0_00 = 0d0; dvardx0_p1 = 0d0
            if (i_var > s% nvar_hydro) then 
               dvardx0_00 = 0d0
            else if (i_var == s% i_lnd) then
               dvardx0_00 = s% d_nonnucneu_dlnd(k)
            else if (i_var == s% i_lnT) then
               dvardx0_00 = s% d_nonnucneu_dlnT(k)
            end if
         end subroutine get_non_nuc_neu_partials
         
         
         subroutine get_gradT_partials(s, &
               k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               dvardx0_m1, dvardx0_00, dvardx0_p1)
            type (star_info), pointer :: s
            integer, intent(in) :: k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index
            real(dp), intent(out) :: dvardx0_m1, dvardx0_00, dvardx0_p1
            dvardx0_m1 = 0d0; dvardx0_00 = 0d0; dvardx0_p1 = 0d0
            if (i_var > s% nvar_hydro) then 
               dvardx0_00 = 0d0
            else if (i_var == s% i_lnd) then
               dvardx0_m1 = s% gradT_ad(k)%d1Array(i_lnd_m1)
               dvardx0_00 = s% gradT_ad(k)%d1Array(i_lnd_00)
            else if (i_var == s% i_lnT) then
               dvardx0_m1 = s% gradT_ad(k)%d1Array(i_lnT_m1)
               dvardx0_00 = s% gradT_ad(k)%d1Array(i_lnT_00)
            else if (i_var == s% i_lnR) then
               dvardx0_00 = s% gradT_ad(k)%d1Array(i_lnR_00)
            else if (i_var == s% i_lum) then
               dvardx0_00 = s% gradT_ad(k)%d1Array(i_L_00)
            else if (i_var == s% i_ln_cvpv0) then
               dvardx0_00 = s% gradT_ad(k)%d1Array(i_xtra1_00)
            else if (i_var == s% i_w_div_wc) then
               dvardx0_00 = s% gradT_ad(k)%d1Array(i_xtra2_00)
            end if
         end subroutine get_gradT_partials
         
         
         subroutine get_mlt_vc_partials(s, &
               k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               dvardx0_m1, dvardx0_00, dvardx0_p1)
            type (star_info), pointer :: s
            integer, intent(in) :: k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index
            real(dp), intent(out) :: dvardx0_m1, dvardx0_00, dvardx0_p1
            dvardx0_m1 = 0d0; dvardx0_00 = 0d0; dvardx0_p1 = 0d0
            if (i_var > s% nvar_hydro) then 
               dvardx0_00 = 0d0
            else if (i_var == s% i_lnd) then
               dvardx0_m1 = s% mlt_vc_ad(k)%d1Array(i_lnd_m1)
               dvardx0_00 = s% mlt_vc_ad(k)%d1Array(i_lnd_00)
            else if (i_var == s% i_lnT) then
               dvardx0_m1 = s% mlt_vc_ad(k)%d1Array(i_lnT_m1)
               dvardx0_00 = s% mlt_vc_ad(k)%d1Array(i_lnT_00)
            else if (i_var == s% i_lnR) then
               dvardx0_00 = s% mlt_vc_ad(k)%d1Array(i_lnR_00)
            else if (i_var == s% i_lum) then
               dvardx0_00 = s% mlt_vc_ad(k)%d1Array(i_L_00)
            end if
         end subroutine get_mlt_vc_partials
         
         
         subroutine get_opacity_partials(s, &
               k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               dvardx0_m1, dvardx0_00, dvardx0_p1)
            type (star_info), pointer :: s
            integer, intent(in) :: k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index
            real(dp), intent(out) :: dvardx0_m1, dvardx0_00, dvardx0_p1
            dvardx0_m1 = 0d0; dvardx0_00 = 0d0; dvardx0_p1 = 0d0
            if (i_var > s% nvar_hydro) then 
               dvardx0_00 = 0d0 ! s% d_opacity_dx(i_var_xa_index,k) - s% d_opacity_dx(i_var_sink_xa_index,k)
            else if (i_var == s% i_lnd) then
               dvardx0_00 = s% d_opacity_dlnd(k)
            else if (i_var == s% i_lnT) then
               dvardx0_00 = s% d_opacity_dlnT(k)
            end if
         end subroutine get_opacity_partials


         subroutine test3_partials(s, &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               k, save_dx, save_equ, ierr)
            type (star_info), pointer :: s
            integer, intent(in) :: &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, k
            real(dp), pointer, dimension(:,:) :: save_dx, save_equ
            integer, intent(out) :: ierr
            real(dp) :: dvardx0_m1, dvardx0_00, dvardx0_p1
            dvardx0_m1 = 0d0
            dvardx0_00 = 0d0
            dvardx0_p1 = 0d0
            if (i_equ > 0) then
               if (i_var > s% nvar_hydro) then ! testing abundance
                  if (k > 1) dvardx0_m1 = &
                     s% lblk(i_equ,i_var,k)/s% x_scale(i_var,k-1) - s% lblk(i_equ,i_var_sink,k)/s% x_scale(i_var_sink,k-1)
                  dvardx0_00 = &
                     s% dblk(i_equ,i_var,k)/s% x_scale(i_var,k) - s% dblk(i_equ,i_var_sink,k)/s% x_scale(i_var_sink,k)
                  if (k < s% nz) dvardx0_p1 = &
                     s% ublk(i_equ,i_var,k)/s% x_scale(i_var,k+1) - s% ublk(i_equ,i_var_sink,k)/s% x_scale(i_var_sink,k+1)
               else
                  if (k > 1) dvardx0_m1 = s% lblk(i_equ,i_var,k)/s% x_scale(i_var,k-1)
                  dvardx0_00 = s% dblk(i_equ,i_var,k)/s% x_scale(i_var,k)
                  if (k < s% nz) dvardx0_p1 = s% ublk(i_equ,i_var,k)/s% x_scale(i_var,k+1)
               end if
            else if (i_equ == -1) then ! 'lnE'
               call get_lnE_partials(s, &
                  k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  dvardx0_m1, dvardx0_00, dvardx0_p1)
            elseif (i_equ == -2) then ! 'eps_nuc'
               call get_eps_nuc_partials(s, &
                  k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  dvardx0_m1, dvardx0_00, dvardx0_p1)
            else if (i_equ == -3) then ! 'opacity'
               call get_opacity_partials(s, &
                  k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  dvardx0_m1, dvardx0_00, dvardx0_p1)
            else if (i_equ == -4) then ! 'lnP'
               call get_lnP_partials(s, &
                  k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  dvardx0_m1, dvardx0_00, dvardx0_p1) 
            else if (i_equ == -5) then ! 'non_nuc_neu'
               call get_non_nuc_neu_partials(s, &
                  k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  dvardx0_m1, dvardx0_00, dvardx0_p1) 
            else if (i_equ == -6) then ! 'gradT'
               call get_gradT_partials(s, &
                  k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  dvardx0_m1, dvardx0_00, dvardx0_p1) 
            else if (i_equ == -7) then ! 'mlt_vc'
               call get_mlt_vc_partials(s, &
                  k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  dvardx0_m1, dvardx0_00, dvardx0_p1) 
            else if (i_equ == -8) then ! 'grad_ad'
               call get_grad_ad_partials(s, &
                  k, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  dvardx0_m1, dvardx0_00, dvardx0_p1) 
            end if 
            if (k > 1) then
               call test1_partial(s, &
                  i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  k, -1, dvardx0_m1, save_dx, save_equ, ierr)
               if (ierr /= 0) stop 'test3_partials'
            end if
            call test1_partial(s, &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               k, 0, dvardx0_00, save_dx, save_equ, ierr)
            if (ierr /= 0) stop 'test3_partials'
            if (k < s% nz) then
               call test1_partial(s, &
                  i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  k, 1, dvardx0_p1, save_dx, save_equ, ierr)
               if (ierr /= 0) stop 'test3_partials'
            end if
         end subroutine test3_partials
         

         subroutine test1_partial(s, &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               k, k_off, dvardx_0, save_dx, save_equ, ierr)
            use chem_def, only: chem_isos
            type (star_info), pointer :: s
            integer, intent(in) :: &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, k, k_off
            real(dp), intent(in) :: dvardx_0
            real(dp), pointer, dimension(:,:) :: save_dx, save_equ
            character (len=3) :: k_off_str
            integer, intent(out) :: ierr 
            character (len = 32) :: equ_str
            real(dp) :: dx_0, err, dvardx, xdum, uncertainty
            include 'formats'
            ierr = 0

            if (i_var > s% nvar_hydro) then ! testing abundance
               dx_0 = s% solver_test_partials_dx_0 * &
                  max(abs(s% xa_start(i_var_xa_index,k) + s% solver_dx(i_var,k)), &
                      abs(s% xa_start(i_var_xa_index,k)), &
                      1d-99)
               write(*,1) 'var name ' // chem_isos% name(s% chem_id(i_var_xa_index))
               write(*,1) 'sink name ' // chem_isos% name(s% chem_id(i_var_sink_xa_index))
            else
               dx_0 = s% solver_test_partials_dx_0 * &
                  max(abs(s% xh_start(i_var,k) + s% solver_dx(i_var,k)), &
                      abs(s% xh_start(i_var,k)))
               if (dx_0 == 0d0) dx_0 = s% solver_test_partials_dx_0
            end if
            dvardx = dfridr(s, &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               k, k_off, dx_0, save_dx, err)
            if (dvardx == 0d0 .and. abs(dvardx_0) < 1d-14) then
               xdum = 0d0
            else if (dvardx_0 == 0d0 .and. abs(dvardx) < 1d-14) then
               xdum = 0d0
            else if (dvardx == 0d0 .or. dvardx_0 == 0d0) then
               xdum = 1d0
            else
               xdum = abs(dvardx - dvardx_0)/min(abs(dvardx),abs(dvardx_0))
            end if
            if (ierr /= 0) then
               write(*,*) 'test1_partial failed'
               stop 'setmatrix'
            end if
            if (i_equ /= 0) then
               if (k_off == 0) then
                  k_off_str = ')  '
               else if (k_off == -1) then
                  k_off_str = '-1)'
               else if (k_off == 1) then
                  k_off_str = '+1)'
               end if
               if (dvardx /= 0d0) then
                  uncertainty = abs(err/dvardx)
               else
                  uncertainty = 0d0
               end if
               if (xdum > 1d-5 .and. uncertainty < 1d-6) then
                  write(*, '(a5,1x)', advance='no') '*****'
               else if (uncertainty > 1d-7) then
                  write(*, '(a5,1x)', advance='no') '?????'
               else
                  write(*, '(6x)', advance='no') 
               end if
               if (i_equ > 0) then
                  equ_str = s% nameofequ(i_equ)
               else if (i_equ == -1) then
                  equ_str = 'lnE'
               else if (i_equ == -2) then
                  equ_str = 'eps_nuc'
               else if (i_equ == -3) then
                  equ_str = 'opacity'
               else if (i_equ == -4) then
                  equ_str = 'lnP'
               else if (i_equ == -5) then
                  equ_str = 'non_nuc_neu'
               else if (i_equ == -6) then
                  equ_str = 'gradT'
               else if (i_equ == -7) then
                  equ_str = 'mlt_vc'
               else if (i_equ == -8) then
                  equ_str = 'grad_ad'
               else
                  equ_str = 'unknown'
               end if
               write(*,'(a25,3x,i5,4x,a,f12.5,4x,a,f12.5,99(4x,a,1pe22.12))') &
                  'd_' // trim(equ_str) // '(k)/d_' // trim(s% nameofvar(i_var)) // &
                  '(k' // trim(k_off_str), &
                  k, 'lg rel diff', safe_log10(xdum), 'lg uncertainty', safe_log10(uncertainty), &
                  'Analytic', dvardx_0, 'Numeric', dvardx
!               write(*,'(a70,2x,i5,f10.3,3x,a,f10.3,99(3x,a,1pe26.16))') &
!                  'log dfridr rel_diff partials wrt  '  // trim(s% nameofvar(i_var)) // &
!                  '(k' // k_off_str // ' of ' // trim(equ_str) // '(k)', &
!                  k, safe_log10(xdum), 'log uncertainty', safe_log10(uncertainty), &
!                  'analytic', dvardx_0, 'numeric', dvardx, &
!                  'analytic/numeric', abs(dvardx_0)/max(1d-99,abs(dvardx))
                  
            else
               write(*,*)
               write(*,1) 'analytic and numeric partials wrt ' // trim(s% nameofvar(i_var)), &
                  dvardx_0, dvardx
               write(*,1) 'log dfridr relative uncertainty for numeric partial', &
                  safe_log10(err/max(1d-99,abs(dvardx)))
               if (dvardx_0 /= 0d0) write(*,1) 'rel_diff', xdum
            end if
         end subroutine test1_partial
            
            
         real(dp) function dfridr_func(s, &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               k, k_off, delta_x, save_dx) result(val)
            type (star_info), pointer :: s
            integer, intent(in) :: &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, k, k_off
            real(dp), intent(in) :: delta_x
            real(dp), pointer, dimension(:,:) :: save_dx
            include 'formats'
            s% solver_dx(i_var,k+k_off) = save_dx(i_var,k+k_off) + delta_x
            if (i_var_xa_index > 0) then ! changing abundance
               !write(*,2) 'new dx, x for abundance', i_var, &
               !     dx(i_var,k+k_off), s% xa(i_var - s% nvar_hydro,k+k_off)
               if (i_var_sink_xa_index <= 0 .or. i_var_sink_xa_index > s% species) then
                  write(*,2) 'bad i_var_sink_xa_index', i_var_sink_xa_index
                  stop 'star_solver dfridr_func'
               end if
               s% solver_dx(i_var_sink,k+k_off) = save_dx(i_var_sink,k+k_off) - delta_x
            end if
            call do_equations(ierr)
            if (ierr /= 0) then
               !exit
               write(*,3) 'call eval_equations failed in dfridr_func'
               stop 'setmatrix'
            end if
            if (i_equ > 0) then
               val = equ(i_equ,k) ! testing partial of residual for cell k equation
            else if (i_equ == 0) then
               val = s% solver_test_partials_val
            else if (i_equ == -1) then
               val = s% lnE(k)
            else if (i_equ == -2) then
               val = s% eps_nuc(k)
            else if (i_equ == -3) then
               val = s% opacity(k)
            else if (i_equ == -4) then
               val = s% lnPeos(k)
            else if (i_equ == -5) then
               val = s% non_nuc_neu(k)
            else if (i_equ == -6) then
               val = s% gradT(k)
            else if (i_equ == -7) then
               val = s% mlt_vc(k)
            else if (i_equ == -8) then
               val = s% grada(k)
            else
               val = 0d0
            end if
            s% solver_dx(i_var,k+k_off) = save_dx(i_var,k+k_off)
            if (i_var_sink > 0) & ! restore sink abundance
               s% solver_dx(i_var_sink,k+k_off) = save_dx(i_var_sink,k+k_off)
         end function dfridr_func


         real(dp) function dfridr(s, &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               k, k_off, hx, save_dx, err)
            type (star_info), pointer :: s
            integer, intent(in) :: &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, k, k_off
            real(dp), intent(in) :: hx
            real(dp), pointer, dimension(:,:) :: save_dx
            real(dp), intent(out) :: err
            !  this routine returns the first derivative of a function func(x)
            !  at the point x, by ridders method of polynomial extrapolation.
            !  value hx is the initial step size;
            !  it should be an increment for which func changes substantially.
            !  an estimate of the error in the first derivative is returned in err.
            integer, parameter :: ntab = 20
            integer :: i,j
            real(dp) :: x,errt,fac,hh,a(ntab,ntab),xdum,ydum,f1,f2
            real(dp), parameter :: con2=2d0, con=sqrt(con2), big=1d50, safe=2d0
            include 'formats'
            dfridr = 0d0
            hh = hx
            ! 2nd order central difference
            f1 = dfridr_func(s, &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               k, k_off, hh, save_dx)
            !write(*,2) 'f1', 1, f1, save_dx(i_var,k) + hh
            f2 = dfridr_func(s, &
               i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
               k, k_off, -hh, save_dx)
            !write(*,2) 'f2', 1, f2, save_dx(i_var,k) - hh
            a(1,1) = (f1 - f2)/(2d0*hh)
            !write(*,2) 'dfdx', 1, a(1,1), &
            !   hh, (save_dx(s% solver_test_partials_var,s% solver_test_partials_k) + hh)/ln10, &
            !   save_dx(s% solver_test_partials_var,s% solver_test_partials_k)/ln10
            err = big
            ! succesive columns in the neville tableu will go to smaller stepsizes
            ! and higher orders of extrapolation
            do i=2,ntab
               hh = hh/con
               f1 = dfridr_func(s, &
                  i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  k, k_off, hh, save_dx)
               !write(*,2) 'f1', i, f1, save_dx(i_var,k) + hh
               f2 = dfridr_func(s, &
                  i_equ, i_var, i_var_sink, i_var_xa_index, i_var_sink_xa_index, &
                  k, k_off, -hh, save_dx)
               !write(*,2) 'f2', i, f2, save_dx(i_var,k) - hh
               a(1,i) = (f1 - f2)/(2d0*hh)
               !write(*,2) 'dfdx', i, a(1,i), &
               !   hh, (save_dx(s% solver_test_partials_var,s% solver_test_partials_k) + hh)/ln10, &
               !   save_dx(s% solver_test_partials_var,s% solver_test_partials_k)/ln10
               ! compute extrapolations of various orders; the error stratagy is to compare
               ! each new extrapolation to one order lower but both at the same stepsize
               ! and at the previous stepsize
               fac = con2
               do j=2,i
                  a(j,i) = (a(j-1,i)*fac - a(j-1,i-1))/(fac-1d0)
                  fac = con2*fac
                  errt = max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
                  if (errt <= err) then
                     err = errt
                     dfridr = a(j,i)
                     !write(*,1) 'dfridr', err/dfridr
                  end if
               end do
               ! if higher order is worse by a significant factor safe, then bail
               if (abs(a(i,i) - a(i-1,i-1)) >= safe*err) then
                  !write(*,1) 'higher order is worse'
                  return
               end if
            end do
         end function dfridr


         subroutine set_xtras(x,num_xtra)
            real(dp) :: x(:,:)
            integer, intent(in) :: num_xtra
            integer :: k
            include 'formats'
            if (.not. s% u_flag) then
               x(1,1:nz) = 0
               x(2,1:nz) = 0
               return
            end if
            do k=1,nz
               x(1,k) = s% u_face_ad(k)%val
               if (is_bad_num(x(1,k))) then
                  write(*,2) 'exit_setmatrix x(1,k)', k, x(1,k)
                  stop
               end if
            end do
            do k=1,nz
               x(2,k) = s% P_face_ad(k)%val
               if (is_bad_num(x(2,k))) then
                  write(*,2) 'exit_setmatrix x(2,k)', k, x(2,k)
                  stop
               end if
            end do
         end subroutine set_xtras
            
         
         subroutine store_mix_type_str(str, integer_string, i, k)
            character (len=5) :: str
            character (len=10) :: integer_string
            integer, intent(in) :: i, k
            integer :: mix_type, j
            if (k < 1 .or. k > s% nz) then
               str(i:i) = 'x'
               return
            end if
            mix_type = s% mixing_type(k)
            if (mix_type < 10) then
               j = mix_type+1
               str(i:i) = integer_string(j:j)
            else
               str(i:i) = '?'
            end if
         end subroutine store_mix_type_str


         subroutine write_msg(msg)
            use const_def, only: secyer
            character(*)  :: msg
            
            integer :: k
            character (len=64) :: max_resid_str, max_corr_str
            character (len=5) :: max_resid_mix_type_str, max_corr_mix_type_str
            character (len=10) :: integer_string
            include 'formats'
            
            if (.not. dbg_msg) return
            
            if (max_resid_j < 0) then
               call sizequ(s, nvar, residual_norm, max_residual, max_resid_k, max_resid_j, ierr)
            end if
            
            if (max_resid_j > 0) then
               write(max_resid_str,*) 'max resid ' // trim(s% nameofequ(max_resid_j))
            else
               max_resid_str = ''
            end if
            
            if (max_corr_j < 0) then
               call sizeB(s, nvar, B, &
                  max_correction, correction_norm, max_corr_k, max_corr_j, ierr)
            end if
            
            if (max_corr_j > 0) then
               write(max_corr_str,*) 'max corr ' // trim(s% nameofvar(max_corr_j))
            else
               max_corr_str = ''
            end if
            
            integer_string = '0123456789'
            k = max_corr_k
            call store_mix_type_str(max_corr_mix_type_str, integer_string, 1, k-2)
            call store_mix_type_str(max_corr_mix_type_str, integer_string, 2, k-1)
            call store_mix_type_str(max_corr_mix_type_str, integer_string, 3, k)
            call store_mix_type_str(max_corr_mix_type_str, integer_string, 4, k+1)
            call store_mix_type_str(max_corr_mix_type_str, integer_string, 5, k+2)
            
            k = max_resid_k
            call store_mix_type_str(max_resid_mix_type_str, integer_string, 1, k-2)
            call store_mix_type_str(max_resid_mix_type_str, integer_string, 2, k-1)
            call store_mix_type_str(max_resid_mix_type_str, integer_string, 3, k)
            call store_mix_type_str(max_resid_mix_type_str, integer_string, 4, k+1)
            call store_mix_type_str(max_resid_mix_type_str, integer_string, 5, k+2)

  111       format(i6, i3, 2x, a, f7.4, &
               1x, a, 1x, e10.3, 2x, a18, 1x, i5, e13.5, a, &
               1x, a, 1x, e10.3, 2x, a16, 1x, i5, e13.5, a, &
               1x, a)
!  111       format(i6, 2x, i3, 2x, a, f8.4, &
!               2x, a, 1x, e10.3, 2x, a19, 1x, i5, e11.3, 2x, a, &
!               2x, a, 1x, e10.3, 2x, a14, 1x, i5, e11.3, 2x, a, &
!               2x, a)
            write(*,111) &
               s% model_number, iter, &
               'coeff', coeff,  &
               'avg resid', residual_norm,  &
!               '   avg resid', residual_norm,  &
               trim(max_resid_str), max_resid_k, max_residual, &
               ' mix type ' // trim(max_resid_mix_type_str),  &
               'avg corr', correction_norm,  &
!               'mix type ' // trim(max_resid_mix_type_str),  &
!               '   avg corr', correction_norm,  &
               trim(max_corr_str), max_corr_k, max_correction,  &
               ' mix type ' // trim(max_corr_mix_type_str),  &
               ' ' // trim(msg)
!               'mix type ' // trim(max_corr_mix_type_str),  &
!               '   ' // trim(msg)
               
            if (is_bad(slope)) stop 'write_msg'

         end subroutine write_msg


         subroutine pointers(ierr)
            integer, intent(out) :: ierr

            integer :: i, j
            character (len=strlen) :: err_msg

            ierr = 0

            i = 1
            A1(1:3*nvar*neq) => work(i:i+3*nvar*neq-1); i = i+3*nvar*neq
            s% equ1(1:neq) => work(i:i+neq-1); i = i+neq
            s% equ(1:nvar,1:nz) => s% equ1(1:neq)
            equ1 => s% equ1
            equ => s% equ
            
            dxsave1(1:neq) => work(i:i+neq-1); i = i+neq
            dxsave(1:nvar,1:nz) => dxsave1(1:neq)
            
            ddxsave1(1:neq) => work(i:i+neq-1); i = i+neq
            ddxsave(1:nvar,1:nz) => ddxsave1(1:neq)
            
            B1 => work(i:i+neq-1); i = i+neq
            B(1:nvar,1:nz) => B1(1:neq)
            
            soln1 => work(i:i+neq-1); i = i+neq
            soln(1:nvar,1:nz) => soln1(1:neq)
            
            grad_f1(1:neq) => work(i:i+neq-1); i = i+neq
            grad_f(1:nvar,1:nz) => grad_f1(1:neq)
            
            rhs1(1:neq) => work(i:i+neq-1); i = i+neq
            rhs(1:nvar,1:nz) => rhs1(1:neq)
            
            xder1(1:neq) => work(i:i+neq-1); i = i+neq
            xder(1:nvar,1:nz) => xder1(1:neq)
            
            ddx1(1:neq) => work(i:i+neq-1); i = i+neq
            ddx(1:nvar,1:nz) => ddx1(1:neq)
            
            row_scale_factors1(1:neq) => work(i:i+neq-1); i = i+neq
            row_scale_factors(1:nvar,1:nz) => row_scale_factors1(1:neq)
            
            col_scale_factors1(1:neq) => work(i:i+neq-1); i = i+neq
            col_scale_factors(1:nvar,1:nz) => col_scale_factors1(1:neq)
            
            save_blks => work(i:i+nvar*neq-1)
            save_ublk1(1:nvar*neq) => save_blks; i = i+nvar*neq
            save_dblk1(1:nvar*neq) => work(i:i+nvar*neq-1); i = i+nvar*neq
            save_lblk1(1:nvar*neq) => work(i:i+nvar*neq-1); i = i+nvar*neq

            if (i-1 > lwork) then
               ierr = -1
               write(*,*) 'use_DGESVX_in_bcyclic', s% use_DGESVX_in_bcyclic
               write(*,  &
                  '(a, i12, a, i12, e26.6)') 'solver: lwork is too small.  must be at least', i-1, &
                  '   but is only ', lwork, dble(i-1 - lwork)/(neq*nvar)
               return
            end if

            i = 1
            ipiv1(1:neq) => iwork(i:i+neq-1); i = i+neq
            if (i-1 > liwork) then
               ierr = -1
               write(*, '(a, i6, a, i6)')  &
                        'solver: liwork is too small.  must be at least', i,  &
                        '   but is only ', liwork
               return
            end if

            A(1:3*nvar,1:neq) => A1(1:3*nvar*neq)
            Acopy1 => A1
            Acopy => A

            ublk1(1:nvar*neq) => A1(1:nvar*neq)
            dblk1(1:nvar*neq) => A1(1+nvar*neq:2*nvar*neq)
            lblk1(1:nvar*neq) => A1(1+2*nvar*neq:3*nvar*neq)

            lblk(1:nvar,1:nvar,1:nz) => lblk1(1:nvar*neq)
            dblk(1:nvar,1:nvar,1:nz) => dblk1(1:nvar*neq)
            ublk(1:nvar,1:nvar,1:nz) => ublk1(1:nvar*neq)

            ublkF1(1:nvar*neq) => AF1(1:nvar*neq)
            dblkF1(1:nvar*neq) => AF1(1+nvar*neq:2*nvar*neq)
            lblkF1(1:nvar*neq) => AF1(1+2*nvar*neq:3*nvar*neq)

            lblkF(1:nvar,1:nvar,1:nz) => lblkF1(1:nvar*neq)
            dblkF(1:nvar,1:nvar,1:nz) => dblkF1(1:nvar*neq)
            ublkF(1:nvar,1:nvar,1:nz) => ublkF1(1:nvar*neq)

         end subroutine pointers


         real(dp) function eval_slope(nvar, nz, grad_f, B)
            integer, intent(in) :: nvar, nz
            real(dp), intent(in), dimension(:,:) :: grad_f, B
            integer :: k, i
            eval_slope = 0
            do i=1,nvar
               eval_slope = eval_slope + dot_product(grad_f(i,1:nz),B(i,1:nz))
            end do
         end function eval_slope


         real(dp) function eval_f(nvar, nz, equ)
            integer, intent(in) :: nvar, nz
            real(dp), intent(in), dimension(:,:) :: equ
            integer :: k, i
            real(dp) :: q
            include 'formats'
            eval_f = 0
            do k = 1, nz
               do i = 1, nvar
                  q = equ(i,k)
                  eval_f = eval_f + q*q
               end do
            end do
            eval_f = eval_f/2
         end function eval_f

      end subroutine do_solver


      subroutine solve_mtx_eqn_with_GMRES( &
            s, nvar, nz, rhs1, soln1, lblk, dblk, ublk, &
            lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, &
            row_scale_factors1, col_scale_factors1, &
            nwrk, wrk1, ipiv1, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar, nz, nwrk
         real(dp), intent(in), pointer :: rhs1(:) ! the right hand side of the linear system
         real(dp), intent(inout), pointer :: soln1(:)
            ! on input, an approximation to the solution. 
            ! on output, an improved approximation.
         real(dp), intent(in), dimension(:,:,:), pointer :: lblk, dblk, ublk ! (nvar,nz)
         real(dp), intent(out), dimension(:), pointer :: wrk1 ! (nwrk)
         real(dp), intent(out), dimension(:), pointer :: &
            lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, row_scale_factors1, col_scale_factors1
         integer, intent(out), dimension(:), pointer :: ipiv1 ! (nwrk)
         integer, intent(out) :: ierr
         integer :: neq, i, itr_max, mr
         real(dp) :: tol_abs, tol_rel
         real(dp), dimension(:), pointer :: &
            work, r, cs, g, sn, y, v1, h1, b1, prod1
         real(dp), dimension(:,:), pointer :: v, h, b, prod
         real(dp), dimension(:,:,:), pointer :: p_diag, p_g, lblkF, dblkF, ublkF
         integer, dimension(:,:), pointer :: ipiv
         character (len=nz) :: equed1
         include 'formats'
         ierr = 0
         itr_max = 20
         tol_abs = 1.0D-08
         tol_rel = 1.0D-08
         neq = nvar*nz
         mr = min(neq-1, 20) 
         work => wrk1
         i = 1
         r(1:neq) => work(i:i+neq-1); i = i+neq
         cs(1:mr) => work(i:i+mr-1); i = i+mr
         sn(1:mr) => work(i:i+mr-1); i = i+mr
         g(1:mr+1) => work(i:i+mr); i = i+mr+1
         y(1:mr+1) => work(i:i+mr); i = i+mr+1
         v1(1:neq*(mr+1)) => work(i:i+neq*(mr+1)-1); i = i+neq*(mr+1)
         v(1:neq,1:mr+1) => v1(1:neq*(mr+1))
         h1(1:(mr+1)*mr) => work(i:i+(mr+1)*mr-1); i = i+(mr+1)*mr
         h(1:mr+1,1:mr) => v1(1:(mr+1)*mr)
         b1(1:neq) => work(i:i+neq-1); i = i+neq
         b(1:nvar,1:nz) => b1(1:neq)
         prod1(1:neq) => work(i:i+neq-1); i = i+neq
         prod(1:nvar,1:nz) => prod1(1:neq)
         
         lblkF(1:nvar,1:nvar,1:nz) => lblkF1(1:nvar*neq)
         dblkF(1:nvar,1:nvar,1:nz) => dblkF1(1:nvar*neq)
         ublkF(1:nvar,1:nvar,1:nz) => ublkF1(1:nvar*neq)
         ipiv(1:nvar,1:nz) => ipiv1(1:neq)
         
         !X => lblkF
         p_diag => dblkF
         p_g => ublkF
         
         if (i > nwrk) then
            write(*,2) 'i', i
            write(*,2) 'nwrk', nwrk
            stop 'i > nwrk in solve_mtx_eqn_with_GMRES'
         end if
         
         call create_preconditioner(ierr)
         if (ierr /= 0) then
            stop 'create_preconditioner failed in solve_mtx_eqn_with_GMRES'
         end if
         
         !call test02 ( )
         call test02_ILU_CR ( )
         
         
         
         copy test02_ILU_CR here and try with mgmres
         
         
         
         
         stop 'solve_mtx_eqn_with_GMRES'
         !call test_preconditioner(ierr)
         !if (ierr /= 0) then
         !   stop 'test_preconditioner failed in solve_mtx_eqn_with_GMRES'
         !end if

         call mgmres ( &
            neq, matvec, apply, soln1, rhs1, itr_max, mr, tol_abs, tol_rel, &
            r, v, cs, g, h, sn, y, ierr )
         
         contains
   
         subroutine test_preconditioner(ierr) 
            integer, intent(out) :: ierr
            integer :: i
            include 'formats'
            write(*,*) 'test_preconditioner'
            soln1 = rhs1
            call apply(soln1) ! soln = precond(rhs)
            call matvec(soln1, b1) ! b = mtx(soln)
            ! for precond = inv(mtx), will have b1 == rhs1
            do i=1,neq
               if (abs(b1(i) - rhs1(i)) > 1d-11+1d-10*abs(rhs1(i))) then
                  write(*,2) 'err b rhs1', i, b1(i) - rhs1(i), b1(i), rhs1(i)
               end if
            end do
            stop 'test_preconditioner'
         end subroutine test_preconditioner
   
         subroutine create_preconditioner(ierr) 
            use star_bcyclic, only: bcyclic_factor
            ! p_diag = preconditioner diagonal blocks
            integer, intent(out) :: ierr
            integer :: i, j, k, op_err
            real(dp) :: work(nvar)
            real(dp), dimension(nvar,nvar) :: temp1, temp2, temp3
            include 'formats'

            do k=1,nz
               call blk_inverse(k, dblk(:,:,k), p_diag(:,:,k), ipiv(:,k), ierr)
            end do
            return




            ! p_diag(:,:,1) = inv(dblk(:,:,1))
            call blk_inverse(k, dblk(:,:,1), p_diag(:,:,1), ipiv(:,1), ierr)
               if (ierr /= 0) then
                  write(*,2) '1', 1
                  stop 'blk_inverse failed in create_preconditioner'
               end if
            ! p_g(:,:,1) = p_diag(:,:,1)*lblk(:,:,2)
            call mm_0(p_diag(:,:,1), lblk(:,:,2), p_g(:,:,1))
            do k=2,nz-1
               call set1_p_diag(k)
               ! p_g(:,:,k) = p_diag(:,:,k)*lblk(:,:,k+1)
               call mm_0(p_diag(:,:,k), lblk(:,:,k+1), p_g(:,:,k))
            end do
            call set1_p_diag(nz)
         end subroutine create_preconditioner
   
         subroutine apply(x) ! apply preconditioner to x
            use star_bcyclic, only: bcyclic_solve
            real(dp), intent(inout) :: x(:) ! (neq)
            integer :: i, k, i1, i2
            real(dp) :: c(neq), temp(nvar)
            include 'formats'
            ierr = 0
            
            do k=1,nz
               i1 = nvar*(k-1)+1; i2 = i1+nvar-1 ! X(:,k) is x(i1:i2)
               call mv_self(p_diag(:,:,k), x(i1:i2)) ! x = a*x
            end do
            return
            
            
            
            
            
            
            col_scale_factors1 = x
            call bcyclic_solve( &
               s, nvar, nz, lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipiv1, &
               col_scale_factors1, row_scale_factors1, row_scale_factors1, col_scale_factors1, equed1, &
               0, ierr)
            if (ierr /= 0) stop 'apply failed in solve_mtx'
            x = row_scale_factors1
            return
            
            
            
            
            
            
            
            ! C(:,1)=p_diag(:,:,1)*X(:,1);
            call mv_0(ublk(:,:,k), x(1:nvar), c(1:nvar))
            do k = 2, nz
               i1 = nvar*(k-1)+1; i2 = i1+nvar-1 ! X(:,k) is x(i1:i2)
               ! C(:,k)=p_diag(:,:,k)*( X(:,k) - ublk(:,:,k-1)*c(:,k-1) );
               temp(:) = x(i1:i2)
               call mv_minus(ublk(:,:,k-1), c(i1-nvar:i2-nvar), temp(:))
               call mv_0(p_diag(:,:,k), temp(:), c(i1:i2))
            end do
            ! x(:,nblk)=C(:,nblk);
            i1 = nvar*(nz-1)+1; i2 = i1+nvar-1 ! X(:,nz) is x(i1:i2)
            x(i1:i2) = c(i1:i2)
            do k = (nz-1), 1, -1
               i1 = nvar*(k-1)+1; i2 = i1+nvar-1 ! X(:,k) is x(i1:i2)
               ! x(:,k)=C(:,k) - p_g(:,:,k)*X(:,k+1);
               x(i1:i2) = c(i1:i2)
               call mv_minus(p_g(:,:,k), x(i1+nvar:i2+nvar), x(i1:i2))
            end do
         end subroutine apply
         
         subroutine blk_inverse(k, blk, blk_inv, ipiv, ierr)
            integer, intent(in) :: k
            real(dp), intent(in) :: blk(:,:)
            real(dp), intent(out) :: blk_inv(:,:)
            integer, intent(out) :: ipiv(:)
            integer, intent(out) :: ierr
            include 'formats'
            blk_inv = blk
            call DGETRF(nvar, nvar, blk_inv, nvar, ipiv, ierr)
            if (ierr /= 0) then
               write(*,3) 'DGETRF failed', k, ierr
               stop 'create_pcond'
            end if
            call DGETRI(nvar, blk_inv, nvar, ipiv, work, nvar, ierr)
            if (ierr /= 0) then
               write(*,3) 'DGETRI failed', k, ierr
               stop 'create_pcond'
            end if         
         end subroutine blk_inverse

         subroutine do_block_dble_mv(nvar, nz, lblk, dblk, ublk, b, prod)
            ! set prod = A*b with A = block tridiagonal given by lblk, dblk, ublk
            integer, intent(in) :: nvar, nz    
            real(dp), pointer, dimension(:,:,:), intent(in) :: lblk, dblk, ublk ! (nvar,nvar,nz)
            real(dp), pointer, dimension(:,:), intent(in) :: b ! (nvar,nz)
            real(dp), pointer, dimension(:,:), intent(inout) :: prod ! (nvar,nz)         
            integer :: k        
            !$OMP PARALLEL DO PRIVATE(k)
            do k = 1, nz
               ! prod(k) = dblk(k)*b(k)
               call mv_0(dblk(:,:,k), b(:,k), prod(:,k))
               if (k > 1) then
                  ! prod(k) = prod(k) + lblk(k)*b(k-1)
                  call mv_plus(lblk(:,:,k), b(:,k-1), prod(:,k))
               end if
               if (k < nz) then
                  ! prod(k) = prod(k) + ublk(k)*b(k+1)
                  call mv_plus(ublk(:,:,k), b(:,k+1), prod(:,k))
               end if
            end do      
            !$OMP END PARALLEL DO         
         end subroutine do_block_dble_mv                  
   
         subroutine matvec(x, r) ! set r = Jacobian*x
            real(dp), intent(in) :: x(:) ! (neq)
            real(dp), intent(out) :: r(:) ! (neq)
            integer :: i
            include 'formats'
            !$omp simd
            do i=1,neq
               b1(i) = x(i)
            end do   
            call do_block_dble_mv(nvar, nz, lblk, dblk, ublk, b, prod)      
            !$omp simd
            do i=1,neq
               r(i) = prod1(i)
            end do                  
         end subroutine matvec
         
         subroutine set1_p_diag(k)
            integer, intent(in) :: k
            integer :: ierr
            real(dp), dimension(nvar,nvar) :: temp
            include 'formats'
            ! temp(:,:) = dblk(:,:) - ublk(:,:,k-1)*p_g(:,:,k-1)
            temp(:,:) = dblk(:,:,k)
            call mm_minus(ublk(:,:,k-1), p_g(:,:,k-1), temp(:,:))      
            ! p_diag(:,:,k)=inv(temp(:,:))
            call blk_inverse(k, temp(:,:), p_diag(:,:,k), ipiv(:,k), ierr)
            if (ierr /= 0) then
               write(*,2) 'nz', nz
               stop 'blk_inverse failed in set1_p_diag'
            end if
         end subroutine set1_p_diag
            
         subroutine mm_0(a, b, c) ! c := a*b, c different than b
            real(dp), dimension(:,:) :: a, b, c ! (nvar,nvar)
            integer :: j, i
            include 'formats'
            do j=1,nvar
               do i=1,nvar
                  c(i,j) = 0d0
               end do
            end do
            call mm_plus(a, b, c)
         end subroutine mm_0      
      
         subroutine mm_plus(a, b, c) ! c := c + a*b
            real(dp), dimension(:,:) :: a, b, c ! (nvar,nvar)
            real(dp) :: tmp
            integer :: j, l, i
            do j = 1,nvar
               do l = 1,nvar
                  tmp = b(l,j)
                  if (tmp /= 0d0) then
                     do i = 1,nvar
                        c(i,j) = c(i,j) + tmp*a(i,l)
                     end do
                  end if
               end do
            end do      
         end subroutine mm_plus
      
         subroutine mm_minus(a, b, c) ! c := c - a*b
            real(dp), dimension(:,:) :: a, b, c ! (nvar,nvar)
            real(dp) :: tmp
            integer :: j, l, i
            do j = 1,nvar
               do l = 1,nvar
                  tmp = b(l,j)
                  if (tmp /= 0d0) then
                     do i = 1,nvar
                        c(i,j) = c(i,j) - tmp*a(i,l)
                     end do
                  end if
               end do
            end do      
         end subroutine mm_minus

         subroutine mv_self(a,x) ! x = a*x
            real(dp) :: a(:,:) ! (nvar,nvar)
            real(dp) :: x(:) ! (nvar)
            real(dp) :: tmp, temp(nvar)
            integer :: j
            !$omp simd
            do j = 1,nvar
               temp(j) = x(j)
            end do
            call mv_0(a,temp,x) 
         end subroutine mv_self

         subroutine mv_0(a,x,y) ! y = a*x, y different than x
            real(dp) :: a(:,:) ! (nvar,nvar)
            real(dp) :: x(:), y(:) ! (nvar)
            real(dp) :: tmp
            integer :: j
            do j = 1,nvar
               y(j) = 0d0
            end do
            call mv_plus(a,x,y)
         end subroutine mv_0

         subroutine mv_plus(a,x,y) ! y = y + a*x, y different than x
            real(dp) :: a(:,:) ! (nvar,nvar)
            real(dp) :: x(:), y(:) ! (nvar)
            real(dp) :: tmp
            integer :: j, i
            do j = 1,nvar
               tmp = x(j)
               if (tmp /= 0d0) then
                  !$omp simd
                  do i = 1,nvar
                     y(i) = y(i) + tmp*a(i,j)
                  end do
               end if
            end do
         end subroutine mv_plus

         subroutine mv_minus(a,x,y) ! y = y - a*x, y different than x
            real(dp) :: a(:,:) ! (nvar,nvar)
            real(dp) :: x(:), y(:) ! (nvar)
            real(dp) :: tmp
            integer :: j, i
            do j = 1,nvar
               tmp = x(j)
               if (tmp /= 0d0) then
                  !$omp simd
                  do i = 1,nvar
                     y(i) = y(i) - tmp*a(i,j)
                  end do
               end if
            end do
         end subroutine mv_minus
            
      end subroutine solve_mtx_eqn_with_GMRES


      subroutine get_solver_work_sizes(s, nvar, nz, lwork, liwork, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar, nz
         integer, intent(out) :: lwork, liwork, ierr
         integer :: neq
         ierr = 0
         neq = nvar*nz
         liwork = neq
         lwork = neq*(7*nvar + 10)
      end subroutine get_solver_work_sizes


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
            r, v, c, g, h, s, y, ierr )
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
         integer, intent(in) :: itr_max, mr
         real(dp), intent(in) :: tol_abs, tol_rel
         real(dp), intent(out), dimension(:) :: r, c, g, s, y
         real(dp), intent(out), dimension(:,:) :: v, h
         integer, intent(out) :: ierr

         real(dp) :: av, mu, rho, rho_tol, htmp
         integer :: i, itr, itr_used, j, k, k_copy
         real(dp), parameter :: delta = 1.0D-03
         logical, parameter :: verbose = .true.
         itr_used = 0
         ierr = 0
         if ( n < mr ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'MGMRES_ST - Fatal error!'
            write ( *, '(a)' ) '  N < MR.'
            write ( *, '(a,i8)' ) '  N = ', n
            write ( *, '(a,i8)' ) '  MR = ', mr
            stop 'MGMRES'
         end if
         do itr = 1, itr_max ! loop back to here for restarts
            call matvec ( x, r )
            !$omp simd
            do j=1,n
               r(j) = rhs(j) - r(j) ! residual = rhs - trial solution
            end do
            call psolve ( r ) ! apply to residual
            rho = sqrt ( dot_product ( r(1:n), r(1:n) ) )
            if ( verbose ) &
               write ( *, '(a,i8,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
            if (is_bad(rho)) stop 'MGMRES bad residual'
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
               call psolve ( v(1:n,k+1) ) ! apply to result of matvec
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
               call mult_givens ( c(k), s(k), k, g )
               rho = abs ( g(k+1) )
               itr_used = itr_used + 1
               if ( verbose ) then
                  write ( *, '(a,i8,a,g14.6)' ) '  K =   ', k, '  Residual = ', rho
               end if
               if (is_bad(rho)) stop 'MGMRES bad residual'
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
            stop 'MGMRES'
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
      
      
      !*****************************************************************************80
      ! 
      ! TESTING
      !
      !*****************************************************************************80
      
      subroutine atx_cr ( n, nz_num, ia, ja, a, x, w )

      !*****************************************************************************80
      !
      !! ATX_CR computes A'*x for a matrix stored in sparse compressed row form.
      !
      !  Discussion:
      !
      !    The Sparse Compressed Row storage format is used.
      !
      !    The matrix A is assumed to be sparse.  To save on storage, only
      !    the nonzero entries of A are stored.  The vector JA stores the
      !    column index of the nonzero value.  The nonzero values are sorted
      !    by row, and the compressed row vector IA then has the property that
      !    the entries in A and JA that correspond to row I occur in indices
      !    IA[I] through IA[I+1]-1.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    17 July 2007
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
      !    Input, integer ( kind = 4 ) N, the order of the system.
      !
      !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
      !
      !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
      !    indices of the matrix values.  The row vector has been compressed.
      !
      !    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
      !
      !    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A'.
      !
      !    Output, real ( kind = 8 ) W(N), the value of A'*X.
      !
        implicit none

        integer ( kind = 4 ) n
        integer ( kind = 4 ) nz_num

        real ( kind = 8 ) a(nz_num)
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ia(n+1)
        integer ( kind = 4 ) ja(nz_num)
        integer ( kind = 4 ) k1
        integer ( kind = 4 ) k2
        real ( kind = 8 ) w(n)
        real ( kind = 8 ) x(n)

        w(1:n) = 0.0D+00

        do i = 1, n
          k1 = ia(i)
          k2 = ia(i+1) - 1
          w(ja(k1:k2)) = w(ja(k1:k2)) + a(k1:k2) * x(i)
        end do

        return
      end subroutine atx_cr
      
      subroutine atx_st ( n, nz_num, ia, ja, a, x, w )

      !*****************************************************************************80
      !
      !! ATX_ST computes A'*x for a matrix stored in sparset triplet form.
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
      !    Input, integer ( kind = 4 ) N, the order of the system.
      !
      !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
      !
      !    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
      !    indices of the matrix values.
      !
      !    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
      !
      !    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A'.
      !
      !    Output, real ( kind = 8 ) W(N), the value of A'*X.
      !
        implicit none

        integer ( kind = 4 ) n
        integer ( kind = 4 ) nz_num

        real ( kind = 8 ) a(nz_num)
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ia(nz_num)
        integer ( kind = 4 ) j
        integer ( kind = 4 ) ja(nz_num)
        integer ( kind = 4 ) k
        real ( kind = 8 ) w(n)
        real ( kind = 8 ) x(n)

        w(1:n) = 0.0D+00

        do k = 1, nz_num
          i = ia(k)
          j = ja(k)
          w(j) = w(j) + a(k) * x(i)
        end do

        return
      end subroutine atx_st
      
      subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )

      !*****************************************************************************80
      !
      !! AX_CR computes A*x for a matrix stored in sparse compressed row form.
      !
      !  Discussion:
      !
      !    The Sparse Compressed Row storage format is used.
      !
      !    The matrix A is assumed to be sparse.  To save on storage, only
      !    the nonzero entries of A are stored.  The vector JA stores the
      !    column index of the nonzero value.  The nonzero values are sorted
      !    by row, and the compressed row vector IA then has the property that
      !    the entries in A and JA that correspond to row I occur in indices
      !    IA[I] through IA[I+1]-1.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    17 July 2007
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
      !    Input, integer ( kind = 4 ) N, the order of the system.
      !
      !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
      !
      !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
      !    indices of the matrix values.  The row vector has been compressed.
      !
      !    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
      !
      !    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
      !
      !    Output, real ( kind = 8 ) W(N), the value of A*X.
      !
        implicit none

        integer ( kind = 4 ) n
        integer ( kind = 4 ) nz_num

        real ( kind = 8 ) a(nz_num)
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ia(n+1)
        integer ( kind = 4 ) ja(nz_num)
        integer ( kind = 4 ) k1
        integer ( kind = 4 ) k2
        real ( kind = 8 ) w(n)
        real ( kind = 8 ) x(n)

        w(1:n) = 0.0D+00

        do i = 1, n
          k1 = ia(i)
          k2 = ia(i+1) - 1
          w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
        end do

        return
      end subroutine ax_cr
      
      subroutine ax_st ( n, nz_num, ia, ja, a, x, w )

      !*****************************************************************************80
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
      !    Input, integer ( kind = 4 ) N, the order of the system.
      !
      !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
      !
      !    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
      !    indices of the matrix values.
      !
      !    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
      !
      !    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
      !
      !    Output, real ( kind = 8 ) W(N), the value of A*X.
      !
        implicit none

        integer ( kind = 4 ) n
        integer ( kind = 4 ) nz_num

        real ( kind = 8 ) a(nz_num)
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ia(nz_num)
        integer ( kind = 4 ) j
        integer ( kind = 4 ) ja(nz_num)
        integer ( kind = 4 ) k
        real ( kind = 8 ) w(n)
        real ( kind = 8 ) x(n)

        w(1:n) = 0.0D+00

        do k = 1, nz_num
          i = ia(k)
          j = ja(k)
          w(i) = w(i) + a(k) * x(j)
        end do

        return
      end subroutine ax_st
      
      subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

      !*****************************************************************************80
      !
      !! DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
      !
      !  Discussion:
      !
      !    The matrix A is assumed to be stored in compressed row format.  Only
      !    the nonzero entries of A are stored.  The vector JA stores the
      !    column index of the nonzero value.  The nonzero values are sorted
      !    by row, and the compressed row vector IA then has the property that
      !    the entries in A and JA that correspond to row I occur in indices
      !    IA[I] through IA[I+1]-1.
      !
      !    The array UA can be used to locate the diagonal elements of the matrix.
      !
      !    It is assumed that every row of the matrix includes a diagonal element,
      !    and that the elements of each row have been ascending sorted.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    18 July 2007
      !
      !  Author:
      !
      !    Original C version by Lili Ju.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) N, the order of the system.
      !
      !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
      !
      !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
      !    indices of the matrix values.  The row vector has been compressed.
      !    On output, the order of the entries of JA may have changed because of
      !    the sorting.
      !
      !    Output, integer ( kind = 4 ) UA(N), the index of the diagonal element
      !    of each row.
      !
        implicit none

        integer ( kind = 4 ) n
        integer ( kind = 4 ) nz_num

        integer ( kind = 4 ) i
        integer ( kind = 4 ) ia(n+1)
        integer ( kind = 4 ) k
        integer ( kind = 4 ) ja(nz_num)
        integer ( kind = 4 ) ua(n)

        ua(1:n) = -1

        do i = 1, n
          do k = ia(i), ia(i+1) - 1
            if ( ja(k) == i ) then
              ua(i) = k
            end if
          end do
        end do

        return
      end subroutine diagonal_pointer_cr
      
      subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )

      !*****************************************************************************80
      !
      !! ILU_CR computes the incomplete LU factorization of a matrix.
      !
      !  Discussion:
      !
      !    The matrix A is assumed to be stored in compressed row format.  Only
      !    the nonzero entries of A are stored.  The vector JA stores the
      !    column index of the nonzero value.  The nonzero values are sorted
      !    by row, and the compressed row vector IA then has the property that
      !    the entries in A and JA that correspond to row I occur in indices
      !    IA(I) through IA(I+1)-1.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    27 July 2007
      !
      !  Author:
      !
      !    Original C version by Lili Ju.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) N, the order of the system.
      !
      !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
      !
      !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
      !    indices of the matrix values.  The row vector has been compressed.
      !
      !    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
      !
      !    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
      !    of each row.
      !
      !    Output, real ( kind = 8 ) L(NZ_NUM), the ILU factorization of A.
      !
        implicit none

        integer ( kind = 4 ) n
        integer ( kind = 4 ) nz_num

        real ( kind = 8 ) a(nz_num)
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ia(n+1)
        integer ( kind = 4 ) iw(n)
        integer ( kind = 4 ) j
        integer ( kind = 4 ) ja(nz_num)
        integer ( kind = 4 ) jj
        integer ( kind = 4 ) jrow
        integer ( kind = 4 ) jw
        integer ( kind = 4 ) k
        real ( kind = 8 ) l(nz_num)
        real ( kind = 8 ) tl
        integer ( kind = 4 ) ua(n)
      !
      !  Copy A.
      !
        l(1:nz_num) = a(1:nz_num)

        do i = 1, n
      !
      !  IW points to the nonzero entries in row I.
      !
          iw(1:n) = -1

          do k = ia(i), ia(i+1) - 1
            iw(ja(k)) = k
          end do

          do j = ia(i), ia(i+1) - 1
            jrow = ja(j)
            if ( i <= jrow ) then
              exit
            end if
            tl = l(j) * l(ua(jrow))
            l(j) = tl
            do jj = ua(jrow) + 1, ia(jrow+1) - 1
              jw = iw(ja(jj))
              if ( jw /= -1 ) then
                l(jw) = l(jw) - tl * l(jj)
              end if
            end do
          end do

          ua(i) = j

          if ( jrow /= i ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'ILU_CR - Fatal error!'
            write ( *, '(a)' ) '  JROW ~= I'
            write ( *, '(a,i8)' ) '  JROW = ', jrow
            write ( *, '(a,i8)' ) '  I    = ', i
            stop
          end if

          if ( l(j) == 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'ILU_CR - Fatal error!'
            write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
            write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
            stop
          end if

          l(j) = 1.0D+00 / l(j)

        end do

        l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

        return
      end subroutine ilu_cr
      
      subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )

      !*****************************************************************************80
      !
      !! LUS_CR applies the incomplete LU preconditioner.
      !
      !  Discussion:
      !
      !    The linear system M * Z = R is solved for Z.  M is the incomplete
      !    LU preconditioner matrix, and R is a vector supplied by the user.
      !    So essentially, we're solving L * U * Z = R.
      !
      !    The matrix A is assumed to be stored in compressed row format.  Only
      !    the nonzero entries of A are stored.  The vector JA stores the
      !    column index of the nonzero value.  The nonzero values are sorted
      !    by row, and the compressed row vector IA then has the property that
      !    the entries in A and JA that correspond to row I occur in indices
      !    IA(I) through IA(I+1)-1.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    18 July 2007
      !
      !  Author:
      !
      !    Original C version by Lili Ju.
      !    FORTRAN90 version by John Burkardt.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) N, the order of the system.
      !
      !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
      !
      !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column
      !    indices of the matrix values.  The row vector has been compressed.
      !
      !    Input, real ( kind = 8 ) L(NZ_NUM), the matrix values.
      !
      !    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
      !    of each row.
      !
      !    Input, real ( kind = 8 ) R(N), the right hand side.
      !
      !    Output, real ( kind = 8 ) Z(N), the solution of the system M * Z = R.
      !
        implicit none

        integer ( kind = 4 ) n
        integer ( kind = 4 ) nz_num

        integer ( kind = 4 ) i
        integer ( kind = 4 ) ia(n+1)
        integer ( kind = 4 ) j
        integer ( kind = 4 ) ja(nz_num)
        real ( kind = 8 ) l(nz_num)
        real ( kind = 8 ) r(n)
        integer ( kind = 4 ) ua(n)
        real ( kind = 8 ) w(n)
        real ( kind = 8 ) z(n)
      !
      !  Copy R in.
      !
        w(1:n) = r(1:n)
      !
      !  Solve L * w = w where L is unit lower triangular.
      !
        do i = 2, n
          do j = ia(i), ua(i) - 1
            w(i) = w(i) - l(j) * w(ja(j))
          end do
        end do
      !
      !  Solve U * w = w, where U is upper triangular.
      !
        do i = n, 1, -1
          do j = ua(i) + 1, ia(i+1) - 1
            w(i) = w(i) - l(j) * w(ja(j))
          end do
          w(i) = w(i) / l(ua(i))
        end do
      !
      !  Copy Z out.
      !
        z(1:n) = w(1:n)

        return
      end subroutine lus_cr
      
      subroutine mgmres_st ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, tol_abs, &
        tol_rel )

      !*****************************************************************************80
      !
      !! MGMRES_ST applies restarted GMRES to a sparse triplet matrix.
      !
      !  Discussion:
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
      !    Input, integer ( kind = 4 ) N, the order of the linear system.
      !
      !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero matrix values.
      !
      !    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
      !    indices of the matrix values.
      !
      !    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
      !
      !    Input/output, real ( kind = 8 ) X(N); on input, an approximation to
      !    the solution.  On output, an improved approximation.
      !
      !    Input, real ( kind = 8 ) RHS(N), the right hand side of the linear system.
      !
      !    Input, integer ( kind = 4 ) ITR_MAX, the maximum number of (outer)
      !    iterations to take.
      !
      !    Input, integer ( kind = 4 ) MR, the maximum number of (inner) iterations
      !    to take.  0 < MR <= N.
      !
      !    Input, real ( kind = 8 ) TOL_ABS, an absolute tolerance applied to the
      !    current residual.
      !
      !    Input, real ( kind = 8 ) TOL_REL, a relative tolerance comparing the
      !    current residual to the initial residual.
      !
        implicit none

        integer ( kind = 4 ) mr
        integer ( kind = 4 ) n
        integer ( kind = 4 ) nz_num

        real ( kind = 8 ) a(nz_num)
        real ( kind = 8 ) av
        real ( kind = 8 ) c(1:mr)
        real ( kind = 8 ), parameter :: delta = 1.0D-03
        real ( kind = 8 ) g(1:mr+1)
        real ( kind = 8 ) h(1:mr+1,1:mr)
        real ( kind = 8 ) htmp
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ia(nz_num)
        integer ( kind = 4 ) itr
        integer ( kind = 4 ) itr_max
        integer ( kind = 4 ) itr_used
        integer ( kind = 4 ) j
        integer ( kind = 4 ) ja(nz_num)
        integer ( kind = 4 ) k
        integer ( kind = 4 ) k_copy
        real ( kind = 8 ) mu
        real ( kind = 8 ) r(1:n)
        real ( kind = 8 ) rho
        real ( kind = 8 ) rho_tol
        real ( kind = 8 ) rhs(1:n)
        real ( kind = 8 ) s(1:mr)
        real ( kind = 8 ) tol_abs
        real ( kind = 8 ) tol_rel
        real ( kind = 8 ) v(1:n,1:mr+1)
        logical, parameter :: verbose = .true.
        real ( kind = 8 ) x(1:n)
        real ( kind = 8 ) y(1:mr+1)

        itr_used = 0

        if ( n < mr ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MGMRES_ST - Fatal error!'
          write ( *, '(a)' ) '  N < MR.'
          write ( *, '(a,i8)' ) '  N = ', n
          write ( *, '(a,i8)' ) '  MR = ', mr
          stop
        end if

        do itr = 1, itr_max

          call ax_st ( n, nz_num, ia, ja, a, x, r )

          r(1:n) = rhs(1:n) - r(1:n)

          rho = sqrt ( dot_product ( r(1:n), r(1:n) ) )

          if ( verbose ) then
            write ( *, '(a,i8,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
          end if

          if ( itr == 1 ) then
            rho_tol = rho * tol_rel
          end if

          v(1:n,1) = r(1:n) / rho

          g(1) = rho
          g(2:mr+1) = 0.0D+00

          h(1:mr+1,1:mr) = 0.0D+00

          do k = 1, mr

            k_copy = k

            call ax_st ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) )

            av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

            do j = 1, k
              h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
              v(1:n,k+1) = v(1:n,k+1) - h(j,k) * v(1:n,j)
            end do

            h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

            if ( av + delta * h(k+1,k) == av ) then

              do j = 1, k
                htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
                h(j,k) = h(j,k) + htmp
                v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
              end do

              h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

            end if

            if ( h(k+1,k) /= 0.0D+00 ) then
              v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
            end if

            if ( 1 < k ) then

              y(1:k+1) = h(1:k+1,k)

              do j = 1, k - 1
                call mult_givens ( c(j), s(j), j, y(1:k+1) )
              end do

              h(1:k+1,k) = y(1:k+1)

            end if

            mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )
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

            if ( rho <= rho_tol .and. rho <= tol_abs ) then
              exit
            end if

          end do

          k = k_copy - 1

          y(k+1) = g(k+1) / h(k+1,k+1)

          do i = k, 1, -1
            y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
          end do

          do i = 1, n
            x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
          end do

          if ( rho <= rho_tol .and. rho <= tol_abs ) then
            exit
          end if

        end do

        if ( verbose ) then
          write ( *, '(a)'       ) ' '
          write ( *, '(a)'       ) 'MGMRES_ST:'
          write ( *, '(a,i8)'    ) '  Iterations = ', itr_used
          write ( *, '(a,g14.6)' ) '  Final residual = ', rho
        end if

        return
      end subroutine mgmres_st
      
      subroutine mult_givens ( c, s, k, g )

      !*****************************************************************************80
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
      !    Input, real ( kind = 8 ) C, S, the cosine and sine of a Givens
      !    rotation.
      !
      !    Input, integer ( kind = 4 ) K, indicates the location of the first
      !    vector entry.
      !
      !    Input/output, real ( kind = 8 ) G(1:K+1), the vector to be modified.
      !    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
      !
        implicit none

        integer ( kind = 4 ) k

        real ( kind = 8 ) c
        real ( kind = 8 ) g(1:k+1)
        real ( kind = 8 ) g1
        real ( kind = 8 ) g2
        real ( kind = 8 ) s

        g1 = c * g(k) - s * g(k+1)
        g2 = s * g(k) + c * g(k+1)

        g(k)   = g1
        g(k+1) = g2

        return
      end subroutine mult_givens
      
      subroutine pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, &
        tol_abs, tol_rel )

      !*****************************************************************************80
      !
      !! PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
      !
      !  Discussion:
      !
      !    The matrix A is assumed to be stored in compressed row format.  Only
      !    the nonzero entries of A are stored.  The vector JA stores the
      !    column index of the nonzero value.  The nonzero values are sorted
      !    by row, and the compressed row vector IA then has the property that
      !    the entries in A and JA that correspond to row I occur in indices
      !    IA(I) through IA(I+1)-1.
      !
      !    This routine uses the incomplete LU decomposition for the
      !    preconditioning.  This preconditioner requires that the sparse
      !    matrix data structure supplies a storage position for each diagonal
      !    element of the matrix A, and that each diagonal element of the
      !    matrix A is not zero.
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
      !    28 August 2012
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
      !    Input, integer ( kind = 4 ) N, the order of the linear system.
      !
      !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero matrix values.
      !
      !    Input, integer ( kind = 4 ) IA(N+1), JA(NZ_NUM), the row and column indices
      !    of the matrix values.  The row vector has been compressed.
      !
      !    Input, real ( kind = 8 ) A(NZ_NUM), the matrix values.
      !
      !    Input/output, real ( kind = 8 ) X(N); on input, an approximation to
      !    the solution.  On output, an improved approximation.
      !
      !    Input, real ( kind = 8 ) RHS(N), the right hand side of the linear system.
      !
      !    Input, integer ( kind = 4 ) ITR_MAX, the maximum number of (outer) 
      !    iterations to take.
      !
      !    Input, integer ( kind = 4 ) MR, the maximum number of (inner) iterations 
      !    to take.  MR must be less than N.
      !
      !    Input, real ( kind = 8 ) TOL_ABS, an absolute tolerance applied to the
      !    current residual.
      !
      !    Input, real ( kind = 8 ) TOL_REL, a relative tolerance comparing the
      !    current residual to the initial residual.
      !
        implicit none

        integer ( kind = 4 ) mr
        integer ( kind = 4 ) n
        integer ( kind = 4 ) nz_num

        real ( kind = 8 ) a(nz_num)
        real ( kind = 8 ) av
        real ( kind = 8 ) c(mr+1)
        real ( kind = 8 ), parameter :: delta = 1.0D-03
        real ( kind = 8 ) g(mr+1)
        real ( kind = 8 ) h(mr+1,mr)
        real ( kind = 8 ) htmp
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ia(n+1)
        integer ( kind = 4 ) itr
        integer ( kind = 4 ) itr_max
        integer ( kind = 4 ) itr_used
        integer ( kind = 4 ) j
        integer ( kind = 4 ) ja(nz_num)
        integer ( kind = 4 ) k
        integer ( kind = 4 ) k_copy
        real ( kind = 8 ) l(ia(n+1)+1)
        real ( kind = 8 ) mu
        real ( kind = 8 ) r(n)
        real ( kind = 8 ) rho
        real ( kind = 8 ) rho_tol
        real ( kind = 8 ) rhs(n)
        real ( kind = 8 ) s(mr+1)
        real ( kind = 8 ) tol_abs
        real ( kind = 8 ) tol_rel
        integer ( kind = 4 ) ua(n)
        real ( kind = 8 ) v(n,mr+1);
        logical, parameter :: verbose = .true.
        real ( kind = 8 ) x(n)
        real ( kind = 8 ) y(mr+1)

        itr_used = 0

        call rearrange_cr ( n, nz_num, ia, ja, a )

        call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

        call ilu_cr ( n, nz_num, ia, ja, a, ua, l )

        if ( verbose ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PMGMRES_ILU_CR'
          write ( *, '(a,i4)' ) '  Number of unknowns = ', n
        end if

        do itr = 1, itr_max

          call ax_cr ( n, nz_num, ia, ja, a, x, r )

          r(1:n) = rhs(1:n) - r(1:n)

          call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )

          rho = sqrt ( dot_product ( r, r ) )

          if ( verbose ) then
            write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
          end if

          if ( itr == 1 ) then
            rho_tol = rho * tol_rel
          end if

          v(1:n,1) = r(1:n) / rho

          g(1) = rho
          g(2:mr+1) = 0.0D+00

          h(1:mr+1,1:mr) = 0.0D+00

          do k = 1, mr

            k_copy = k

            call ax_cr ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) ) 

            call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )

            av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

            do j = 1, k
              h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
              v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
            end do

            h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

            if ( ( av + delta * h(k+1,k)) == av ) then
              do j = 1, k
                htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
                h(j,k) = h(j,k) + htmp
                v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
              end do
              h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
            end if

            if ( h(k+1,k) /= 0.0D+00 ) then
              v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
            end if

            if ( 1 < k ) then
              y(1:k+1) = h(1:k+1,k)
              do j = 1, k - 1
                call mult_givens ( c(j), s(j), j, y )
              end do
              h(1:k+1,k) = y(1:k+1)
            end if

            mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

            c(k) = h(k,k) / mu
            s(k) = -h(k+1,k) / mu
            h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
            h(k+1,k) = 0.0D+00
            call mult_givens ( c(k), s(k), k, g )

            rho = abs ( g(k+1) )

            itr_used = itr_used + 1

            if ( verbose ) then
              write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
            end if

            if ( rho <= rho_tol .and. rho <= tol_abs ) then
              exit
            end if

          end do

          k = k_copy - 1

          y(k+1) = g(k+1) / h(k+1,k+1)

          do i = k, 1, -1
            y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
          end do

          do i = 1, n
            x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
          end do

          if ( rho <= rho_tol .and. rho <= tol_abs ) then
            exit
          end if

        end do

        if ( verbose ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
          write ( *, '(a,i6)' ) '  Iterations = ', itr_used
          write ( *, '(a,g14.6)' ) '  Final residual = ', rho
        end if

        return
      end subroutine pmgmres_ilu_cr
      
      subroutine r8vec_uniform_01 ( n, seed, r )

      !*****************************************************************************80
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
        implicit none

        integer ( kind = 4 ) n

        integer ( kind = 4 ) i
        integer ( kind = 4 ) k
        integer ( kind = 4 ) seed
        real ( kind = 8 ) r(n)

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

        return
      end subroutine r8vec_uniform_01
      
      subroutine rearrange_cr ( n, nz_num, ia, ja, a )

      !*****************************************************************************80
      !
      !! REARRANGE_CR sorts a sparse compressed row matrix.
      !
      !  Discussion:
      !
      !    This routine guarantees that the entries in the CR matrix
      !    are properly sorted.
      !
      !    After the sorting, the entries of the matrix are rearranged in such
      !    a way that the entries of each column are listed in ascending order
      !    of their column values.
      !
      !    The matrix A is assumed to be stored in compressed row format.  Only
      !    the nonzero entries of A are stored.  The vector JA stores the
      !    column index of the nonzero value.  The nonzero values are sorted
      !    by row, and the compressed row vector IA then has the property that
      !    the entries in A and JA that correspond to row I occur in indices
      !    IA(I) through IA(I+1)-1.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    17 July 2007
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
      !    Input, integer ( kind = 4 ) N, the order of the system.
      !
      !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
      !
      !    Input, integer ( kind = 4 ) IA(N+1), the compressed row indices.
      !
      !    Input/output, integer ( kind = 4 ) JA(NZ_NUM), the column indices.
      !    On output, these may have been rearranged by the sorting.
      !
      !    Input/output, real ( kind = 8 ) A(NZ_NUM), the matrix values.  On output,
      !    the matrix values may have been moved somewhat because of the sorting.
      !
        implicit none

        integer ( kind = 4 ) n
        integer ( kind = 4 ) nz_num

        real ( kind = 8 ) a(nz_num)
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ia(n+1)
        integer ( kind = 4 ) i4temp
        integer ( kind = 4 ) ja(nz_num)
        integer ( kind = 4 ) k
        integer ( kind = 4 ) l
        real ( kind = 8 ) r8temp

        do i = 1, n

          do k = ia(i), ia(i+1) - 2
            do l = k + 1, ia(i+1) - 1

              if ( ja(l) < ja(k) ) then
                i4temp = ja(l)
                ja(l)  = ja(k)
                ja(k)  = i4temp

                r8temp = a(l)
                a(l)   = a(k)
                a(k)   = r8temp
              end if

            end do
          end do

        end do

        return
      end subroutine rearrange_cr
      
      subroutine timestamp ( )

      !*****************************************************************************80
      !
      !! TIMESTAMP prints the current YMDHMS date as a time stamp.
      !
      !  Example:
      !
      !    31 May 2001   9:45:54.872 AM
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    18 May 2013
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    None
      !
        implicit none

        character ( len = 8 ) ampm
        integer ( kind = 4 ) d
        integer ( kind = 4 ) h
        integer ( kind = 4 ) m
        integer ( kind = 4 ) mm
        character ( len = 9 ), parameter, dimension(12) :: month = (/ &
          'January  ', 'February ', 'March    ', 'April    ', &
          'May      ', 'June     ', 'July     ', 'August   ', &
          'September', 'October  ', 'November ', 'December ' /)
        integer ( kind = 4 ) n
        integer ( kind = 4 ) s
        integer ( kind = 4 ) values(8)
        integer ( kind = 4 ) y

        call date_and_time ( values = values )

        y = values(1)
        m = values(2)
        d = values(3)
        h = values(5)
        n = values(6)
        s = values(7)
        mm = values(8)

        if ( h < 12 ) then
          ampm = 'AM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Noon'
          else
            ampm = 'PM'
          end if
        else
          h = h - 12
          if ( h < 12 ) then
            ampm = 'PM'
          else if ( h == 12 ) then
            if ( n == 0 .and. s == 0 ) then
              ampm = 'Midnight'
            else
              ampm = 'AM'
            end if
          end if
        end if

        write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
          d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

        return
      end subroutine timestamp
      
      subroutine test02 ( )

      !*****************************************************************************80
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

        integer ( kind = 4 ), parameter :: n = 9
        integer ( kind = 4 ), parameter :: nz_num = 23

        real ( kind = 8 ), dimension(nz_num) :: a = (/ &
           2.0D+00, -1.0D+00, &
           2.0D+00, -1.0D+00, &
          -1.0D+00,  2.0D+00, &
          -1.0D+00,  2.0D+00, -1.0D+00, &
          -1.0D+00,  2.0D+00, -1.0D+00, &
          -1.0D+00,  2.0D+00, -1.0D+00, &
          -1.0D+00,  2.0D+00, -1.0D+00, &
          -1.0D+00,  2.0D+00, -1.0D+00, &
          -1.0D+00,  2.0D+00 /)
        integer ( kind = 4 ) i
        integer ( kind = 4 ), dimension(nz_num) :: ia = (/ &
          1, 1, &
          2, 2, &
          3, 3, &
          4, 4, 4, &
          5, 5, 5, &
          6, 6, 6, &
          7, 7, 7, &
          8, 8, 8, &
          9, 9 /)
        integer ( kind = 4 ) itr_max
        integer ( kind = 4 ), dimension(nz_num) :: ja = (/ &
          1, 4, &
          2, 3, &
          2, 3, &
          1, 4, 5, &
          4, 5, 6, &
          5, 6, 7, &
          6, 7, 8, &
          7, 8, 9, &
          8, 9 /)
        integer ( kind = 4 ) mr
        real ( kind = 8 ), dimension(n) :: rhs = (/ &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00 /)
        integer ( kind = 4 ) :: seed = 123456789
        integer ( kind = 4 ) test
        real ( kind = 8 ) tol_abs
        real ( kind = 8 ) tol_rel
        real ( kind = 8 ) x_error
        real ( kind = 8 ) x_estimate(n)
        real ( kind = 8 ), dimension(n) :: x_exact = (/ &
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

          x_error = sqrt ( sum ( ( x_exact(1:n) - x_estimate(1:n) )**2 ) )

          write ( *, '(a,g14.6)' ) '  Before solving, X_ERROR = ', x_error

          itr_max = 20
          mr = n - 1
          tol_abs = 1.0D-08
          tol_rel = 1.0D-08

          call mgmres_st ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr, &
            tol_abs, tol_rel )
          !call pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, &
          !  mr, tol_abs, tol_rel )

          x_error = sqrt ( sum ( ( x_exact(1:n) - x_estimate(1:n) )**2 ) )

          write ( *, '(a,g14.6)' ) '  After solving, X_ERROR = ', x_error

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Final solution estimate:'
          write ( *, '(a)' ) ' '
          do i = 1, n
            write ( *, '(2x,i8,2x,g14.6)' ) i, x_estimate(i)
          end do

        end do

        return
      end subroutine test02
      
      
      subroutine test03 ( )

      !*****************************************************************************80
      !
      !! TEST03 tests PMGMRES_ILU_CR on the simple -1,2-1 matrix.
      !
      !  Discussion:
      !
      !    This is a very weak test, since the matrix has such a simple
      !    structure, is diagonally dominant (though not strictly), 
      !    and is symmetric.
      !
      !    To make the matrix bigger, simply increase the value of N.
      !
      !    Note that PGMRES_ILU_CR expects the matrix to be stored using the
      !    sparse compressed row format.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license. 
      !
      !  Modified:
      !
      !    28 August 2012
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 4 ), parameter :: n = 20
        integer ( kind = 4 ), parameter :: nz_num = ( 3 * n - 2 )

        real ( kind = 8 ) a(nz_num)
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ia(n+1)
        integer ( kind = 4 ) itr_max
        integer ( kind = 4 ) ja(nz_num)
        integer ( kind = 4 ) k
        integer ( kind = 4 ) mr
        real ( kind = 8 ) rhs(n)
        integer ( kind = 4 ) test
        real ( kind = 8 ) tol_abs
        real ( kind = 8 ) tol_rel
        real ( kind = 8 ) x_error
        real ( kind = 8 ) x_estimate(n)
        real ( kind = 8 ) x_exact(n)

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST03'
        write ( *, '(a)' ) '  Test PMGMRES_ILU_CR on the simple -1,2-1 matrix.'
      !
      !  Set the matrix.
      !  Note that we use 1-based index values in IA and JA.
      !
        k = 1
        ia(1) = 1

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4,a,i4)' ) '  ia(', 1, ') = ', ia(1)

        do i = 1, n

          ia(i+1) = ia(i)

          if ( 1 < i ) then
            ia(i+1) = ia(i+1) + 1
            ja(k) = i - 1
            a(k) = -1.0D+00
            k = k + 1
          end if

          ia(i+1) = ia(i+1) + 1
          ja(k) = i
          a(k) = 2.0D+00
          k = k + 1

          if ( i < n ) then
            ia(i+1) = ia(i+1) + 1
            ja(k) = i + 1
            a(k) = -1.0D+00
            k = k + 1
          end if
          write ( *, '(a,i4,a,i4)' ) '  ia(', i + 1, ') = ', ia(i+1)
        end do
      !
      !  Set the right hand side:
      !
        rhs(1:n-1) = 0.0D+00
        rhs(n) = real ( n + 1, kind = 8 )
      !
      !  Set the exact solution.
      !
        do i = 1, n
          x_exact(i) = real ( i, kind = 8 )
        end do

        do test = 1, 1 ! 3
      !
      !  Set the initial solution estimate.
      !
          x_estimate(1:n) = 0.0D+00
          x_error = 0.0D+00
          do i = 1, n
            x_error = x_error + ( x_exact(i) - x_estimate(i) ) ** 2
          end do
          x_error = sqrt ( x_error )

          if ( test == 1 ) then
            itr_max = 1
            mr = 20
          else if ( test == 2 ) then
            itr_max = 2
            mr = 10
          else if ( test == 3 ) then
            itr_max = 5
            mr = 4
          end if

          tol_abs = 1.0D-08
          tol_rel = 1.0D-08

          write ( *, '(a)' ) ' '
          write ( *, '(a,i4)' ) '  Test ', test
          write ( *, '(a,i4)' ) '  Matrix order N = ', n
          write ( *, '(a,i4)' ) '  Inner iteration limit = ', mr
          write ( *, '(a,i4)' ) '  Outer iteration limit = ', itr_max
          write ( *, '(a,g14.6)' ) '  Initial X_ERROR = ', x_error

          call pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, &
            mr, tol_abs, tol_rel )

          x_error = 0.0D+00
          do i = 1, n
            x_error = x_error + ( x_exact(i) - x_estimate(i) ) ** 2
          end do
          x_error = sqrt ( x_error )

          write ( *, '(a,g14.6)' ) '  Final X_ERROR = ', x_error

        end do

        return
      end subroutine test03
      
      
      subroutine test04 ( )

      !*****************************************************************************80
      !
      !! TEST04 tests PMGMRES_ILU_CR on a simple 5 by 5 matrix.
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license. 
      !
      !  Modified:
      !
      !    29 August 2012
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 4 ), parameter :: n = 5
        integer ( kind = 4 ), parameter :: nz_num = 9

        real ( kind = 8 ), dimension ( nz_num ) :: a = (/ &
           1.0, 2.0, 1.0, &
           2.0, &
           3.0, 3.0, &
           4.0, &
           1.0, 5.0 /)
        integer ( kind = 4 ) i
        integer ( kind = 4 ), dimension ( n + 1 ) :: ia = (/ 1, 4, 5, 7, 8, 10 /)
        integer ( kind = 4 ) itr_max
        integer ( kind = 4 ), dimension ( nz_num ) :: ja = (/ &
          1, 4, 5, &
          2, &
          1, 3, &
          4, &
          2, 5 /)
        integer ( kind = 4 ) mr
        real ( kind = 8 ), dimension ( n ) :: rhs = (/ &
          14.0, 4.0, 12.0, 16.0, 27.0 /)
        integer ( kind = 4 ) test
        real ( kind = 8 ) tol_abs
        real ( kind = 8 ) tol_rel
        real ( kind = 8 ) x_error
        real ( kind = 8 ) x_estimate(n)
        real ( kind = 8 ), dimension ( n ) :: x_exact = (/ 1.0, 2.0, 3.0, 4.0, 5.0 /)

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST04'
        write ( *, '(a)' ) '  Test PMGMRES_ILU_CR on a simple 5 x 5 matrix.'

        write ( *, '(a)' ) ' '
        do i = 1, n
          write ( *, '(a,i2,a,i2)' ) '  ia(', i, ') = ', ia(i)
        end do
 
        do test = 1, 3
      !
      !  Set the initial solution estimate.
      !
          x_estimate(1:n) = 0.0D+00
          x_error = 0.0D+00
          do i = 1, n
            x_error = x_error + ( x_exact(i) - x_estimate(i) ) ** 2
          end do
          x_error = sqrt ( x_error )

          if ( test == 1 ) then
            itr_max = 1
            mr = 20
          else if ( test == 2 ) then
            itr_max = 2
            mr = 10
          else if ( test == 3 ) then
            itr_max = 5
            mr = 4
          end if

          tol_abs = 1.0D-08
          tol_rel = 1.0D-08

          write ( *, '(a)' ) ' '
          write ( *, '(a,i4)' ) '  Test ', test
          write ( *, '(a,i4)' ) '  Matrix order N = ', n
          write ( *, '(a,i4)' ) '  Inner iteration limit = ', mr
          write ( *, '(a,i4)' ) '  Outer iteration limit = ', itr_max
          write ( *, '(a,g14.6)' ) '  Initial X_ERROR = ', x_error

          call pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, &
            mr, tol_abs, tol_rel )

          x_error = 0.0D+00
          do i = 1, n
            x_error = x_error + ( x_exact(i) - x_estimate(i) ) ** 2
          end do
          x_error = sqrt ( x_error )

          write ( *, '(a,g14.6)' ) '  Final X_ERROR = ', x_error

        end do

        return
      end subroutine test04
      
      
      subroutine test02_ILU_CR ( )
      !*****************************************************************************80
      !
      !! TEST02 modified to use PMGMRES_ILU_CR.    switch from ST to CR.
      !
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

        integer ( kind = 4 ), parameter :: n = 9
        integer ( kind = 4 ), parameter :: nz_num = 23
        ! nonzero values are sorted by row
        real ( kind = 8 ), dimension(nz_num) :: a = (/ &
           2.0D+00, -1.0D+00, & ! row 1
           2.0D+00, -1.0D+00, & ! row 2 
          -1.0D+00,  2.0D+00, & ! row 3
          -1.0D+00,  2.0D+00, -1.0D+00, & ! row 4
          -1.0D+00,  2.0D+00, -1.0D+00, & ! row 5
          -1.0D+00,  2.0D+00, -1.0D+00, & ! row 6
          -1.0D+00,  2.0D+00, -1.0D+00, & ! row 7
          -1.0D+00,  2.0D+00, -1.0D+00, & ! row 8
          -1.0D+00,  2.0D+00 /) ! row 9        
        ! JA stores the column index of the nonzero value
        integer ( kind = 4 ), dimension(nz_num) :: ja = (/ &
          1, 4, & ! row 1
          2, 3, & ! row 2
          2, 3, & ! row 3
          1, 4, 5, & ! row 4
          4, 5, 6, & ! row 5
          5, 6, 7, & ! row 6
          6, 7, 8, & ! row 7
          7, 8, 9, & ! row 8
          8, 9 /) ! row 9          
        ! the entries in A and JA that correspond to row I occur in indices
        ! IA[I] through IA[I+1]-1.
        integer ( kind = 4 ), dimension(n+1) :: ia = (/ &
          1, & ! row 1
          3, & ! row 2
          5, & ! row 3
          7, & ! row 4
          10, & ! row 5
          13, & ! row 6
          16, & ! row 7
          19, & ! row 8
          22, & ! row 9
          24  /)           
        integer ( kind = 4 ) i
        integer ( kind = 4 ) itr_max          
        integer ( kind = 4 ) mr
        real ( kind = 8 ), dimension(n) :: rhs = (/ &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00, &
          1.0D+00 /)
        integer ( kind = 4 ) :: seed = 123456789
        integer ( kind = 4 ) test
        real ( kind = 8 ) tol_abs
        real ( kind = 8 ) tol_rel
        real ( kind = 8 ) x_error
        real ( kind = 8 ) x_estimate(n)
        real ( kind = 8 ), dimension(n) :: x_exact = (/ &
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
        write ( *, '(a)' ) 'TEST02_ILU_CR'
        write ( *, '(a)' ) '  Test pmgmres_ilu_cr on a matrix that is not quite '
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

          x_error = sqrt ( sum ( ( x_exact(1:n) - x_estimate(1:n) )**2 ) )

          write ( *, '(a,g14.6)' ) '  Before solving, X_ERROR = ', x_error

          itr_max = 20
          mr = n - 1
          tol_abs = 1.0D-08
          tol_rel = 1.0D-08

          call pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr, &
            tol_abs, tol_rel )

          x_error = sqrt ( sum ( ( x_exact(1:n) - x_estimate(1:n) )**2 ) )

          write ( *, '(a,g14.6)' ) '  After solving, X_ERROR = ', x_error

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Final solution estimate:'
          write ( *, '(a)' ) ' '
          do i = 1, n
            write ( *, '(2x,i8,2x,g14.6)' ) i, x_estimate(i)
          end do

        end do

        return
      end subroutine test02_ILU_CR
      


      end module star_solver
