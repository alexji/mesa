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


      module abtilu

      use star_private_def
      use const_def, only: dp
      use utils_lib, only: is_bad

      implicit none

      private
      public :: test_abtilu, &
         solve_abtilu_with_Bi_CG_Stab, solve_abtilu_with_mgmres, &
         show_vec, write_MM_mxt, write_MM_vec


      contains      
      
      
      subroutine test_abtilu()
         include 'formats'
         integer :: j
         do j=1,1000
            !call test_Bi_CG_Stab
            call test_abtilu_nvar1(j)
            call test_abtilu_nvar2(j)
         end do
         write(*,*)
         stop 'done test_abtilu'
      end subroutine test_abtilu
      
      
      subroutine solve_abtilu_with_Bi_CG_Stab( &
            nvar, nz, A, ublk, dblk, lblk, rhs1, &
            equilibrate, exact, &
            num_sweeps_factor, num_sweeps_solve, & ! for abtilu
            max_iter, tol, & ! for Bi_CG_Stab
            soln1, verbose, ierr)
         real(dp), dimension(:,:), intent(in) :: A ! (neq,neq) for debugging
         integer, intent(in) :: nvar, nz, max_iter, &
            num_sweeps_factor, num_sweeps_solve
         real(dp), dimension(:,:,:), intent(inout) :: & !(nvar,nvar,nz)
            ublk, dblk, lblk
         real(dp), dimension(:), intent(inout) :: rhs1 ! (neq)
         logical, intent(in) :: equilibrate, exact, verbose
         real(dp), intent(in) :: tol
         real(dp), dimension(:), intent(inout) :: soln1 ! (neq)
            ! input: initial guess (can be 0)
            ! output: final approximation
         integer, intent(out) :: ierr
         
         real(dp), dimension(:,:,:), allocatable :: & !(nvar,nvar,nz)
            Dhat, invDhat, invDhat_lblk, invDhat_ublk
         real(dp), dimension(:,:), allocatable :: & ! (neq,neq) for debugging
            Acopy, AF_abtilu, AF_exact, A_temp, D, E, F, invD, F1
         real(dp), dimension(:), allocatable :: & ! (neq)
            prev1, w1, x1, y1, r, DR, DC, check1_abtilu, check1_exact
         integer, dimension(:,:), allocatable :: ipiv ! (nvar,nz)
         integer, dimension(:), allocatable :: & ! (nvar*nz)
            IPIV_AF_abtilu, IPIV_AF_exact
         integer :: neq, i, j, k, iter
         real(dp) :: error
         logical :: use_dummy_psolve, debug
         include 'formats'
         ierr = 0
         neq = nvar*nz
         debug = .false.

         allocate( &
            Dhat(nvar,nvar,nz), invDhat(nvar,nvar,nz), &
            invDhat_lblk(nvar,nvar,nz), invDhat_ublk(nvar,nvar,nz), &
            Acopy(neq,neq), prev1(neq), w1(neq), x1(neq), y1(neq), r(neq), &
            DR(neq), DC(neq), ipiv(nvar,nz))
         if (debug) allocate( &
            D(neq,neq), E(neq,neq), F(neq,neq), invD(neq,neq), &
            F1(neq,neq), AF_abtilu(neq,neq), IPIV_AF_abtilu(neq), &
            check1_exact(neq), check1_abtilu(neq), &
            A_temp(neq,neq), AF_exact(neq,neq), IPIV_AF_exact(neq))
            
         if (debug) then ! check that A matches ublk, dblk, lblk
            Acopy = 0d0
            call copy_upper_to_square(nvar,nz,ublk,Acopy)
            call copy_diag_to_square(nvar,nz,dblk,Acopy)
            call copy_lower_to_square(nvar,nz,lblk,Acopy)
            do i=1,neq
               do j=1,neq
                  if (abs(A(i,j) - Acopy(i,j)) > 1d-12*(1d-50 + abs(A(i,j)))) then
                     ierr = -1
                     write(*,3) 'A(i,j) Acopy(i,j) bad', i, j, &
                        (A(i,j) - Acopy(i,j))/(1d-50 + abs(A(i,j))), &
                        A(i,j), Acopy(i,j)
                  end if
               end do
            end do
            if (ierr /= 0) stop 'A fails to match ublk, dblk, lblk'
            if (verbose) write(*,*) 'checked A matches ublk, dblk, lblk'
         end if
         
         call create_preconditioner_abtilu(ierr)     
         if (ierr /= 0) stop 'failed in create_preconditioner_abtilu_Bi_CG_Stab'
         
         if (equilibrate) then ! scale rhs1 with DR
            if (verbose .and. nz < 20) then
               write(*,*) 'scale rhs1'
               call show_vec(nvar, nz, rhs1)
            end if
            !$omp simd
            do j=1,neq
               rhs1(j) = rhs1(j)*DR(j)
            end do
            if (verbose .and. nz < 20) call show_vec(nvar, nz, rhs1)
         end if
         
         if (.false.) then ! TESTING. 
            use_dummy_psolve = .true.
            call Bi_CG_Stab( &     
               matvec_abtilu_Bi_CG_Stab, solve_abtilu_Bi_CG_Stab, &
               rhs1, soln1, max_iter, tol, error, iter, verbose, ierr)
            if (ierr /= 0) then
               stop 'failed in Bi_CG_Stab'
            end if
         else
            call Bi_CG_Stab( &     
               matvec_abtilu_Bi_CG_Stab, solve_abtilu_Bi_CG_Stab, &
               rhs1, soln1, max_iter, tol, error, iter, verbose, ierr)
            if (ierr /= 0) then
               stop 'failed in Bi_CG_Stab'
            end if
         end if
         
         if (equilibrate) then ! scale soln1 with DC
            !$omp simd
            do j=1,neq
               soln1(j) = soln1(j)*DC(j)
            end do
         end if
         
         contains
         
         subroutine create_preconditioner_abtilu(ierr)
            integer, intent(out) :: ierr
            integer :: i, j
            real(dp) :: a_exact, a_abtilu
            logical :: okay
            include 'formats'
            if (debug) then
               F = 0d0; E = 0d0
               call copy_upper_to_square(nvar,nz,ublk,F)
               call copy_lower_to_square(nvar,nz,lblk,E)
            end if
            call equilibrate_and_factor_abtilu_right( &
                nvar, nz, num_sweeps_factor, equilibrate, exact, & ! input
                lblk, dblk, ublk, & ! input/output
                Dhat, ipiv, & ! work
                DR, DC, invDhat, invDhat_lblk, invDhat_ublk, &
                verbose .and. nz < 20, ierr) ! output
            if (ierr /= 0) then
               stop 'failed in equilibrate_and_factor_abtilu_right'
            end if
            if (debug .and. exact) then
               ! for exact ILU0 block tridiagonal, A = (Dhat + E)*invDhat*(Dhat + F)  5.8
               ! solve for y: (Dhat + E)*y = r  =>   y = invDhat*(r - E*y) = invDhat*r - invDhat*E*y
               ! solve for z: invDhat*(Dhat + F)*z = y  => z = y - invDhat*F*z
               ! save invDhat*E as invDhat_lblk and invDhat*F as invDhat_ublk
               
               D = 0d0; invD = 0d0; F1 = 0d0; A_temp = 0d0
               call copy_diag_to_square(nvar,nz,Dhat,D)
               call copy_diag_to_square(nvar,nz,invDhat,invD)
               A_temp = D + F
               F1 = matmul(invD,A_temp)
               A_temp = D + E
               AF_abtilu = matmul(A_temp,F1)
               ! for exact factoring, AF_abtilu should = A
               do i=1,neq
                  do j=1,neq  
                     a_exact = A(i,j)
                     a_abtilu = AF_abtilu(i,j)
                     if (abs(a_exact - a_abtilu) > 1d-12*(1d-50 + abs(a_exact))) then
                        ierr = -1
                        write(*,2) 'a_exact a_abtilu bad', j, &
                           (a_exact - a_abtilu)/(1d-50 + abs(a_exact)), &
                           a_exact, a_abtilu
                     end if
                  end do
               end do
               if (ierr /= 0) stop 'failed check in create_preconditioner_abtilu'
               if (verbose) write(*,*) 'checked factor_abtilu'
               ! now factor for use in checking solve_abtilu
               CALL DGETRF( neq, neq, AF_abtilu, neq, IPIV_AF_abtilu, ierr )
               if (ierr /= 0) then
                  stop 'failed in DGETRF for AF_abtilu'
               end if
               AF_exact = A
               CALL DGETRF( neq, neq, AF_exact, neq, IPIV_AF_exact, ierr )
               if (ierr /= 0) then
                  stop 'failed in DGETRF for AF_exact'
               end if
            end if            
         end subroutine create_preconditioner_abtilu
     
         subroutine solve_abtilu_Bi_CG_Stab(r1, z1) ! set z = Precond*r
            real(dp), intent(in) :: r1(:)
            real(dp), intent(out) :: z1(:)
            integer :: j
            logical :: okay
            include 'formats'
            z1 = r1
            if (debug .and. exact) then
               check1_exact = r1
            end if
            call solve_abtilu_right( &
                nvar, nz, num_sweeps_solve, &
                invDhat, invDhat_lblk, invDhat_ublk, lblk, ublk, r1, & ! input
                DR, DC, equilibrate, exact, & ! input
                prev1, w1, x1, y1, & ! work
                z1, verbose, ierr) ! output
            if (ierr /= 0) then
               stop 'failed in solve_abtilu_right'
            end if
            if (debug .and. exact) then
               ! compare solve_abtilu_right to DGETRS with AF_exact
               CALL DGETRS( 'N', neq, 1, AF_exact, neq, IPIV_AF_exact, check1_exact, neq, ierr )
               if (ierr /= 0) then
                  stop 'failed in DGETRS'
               end if
               okay = check_match(check1_exact, z1)
               if (.not. okay) then
                  !z1 = check1_exact      to show that can get good result from Bi_CG_Stab
                  !write(*,*) 'switch to check1_exact in place of z'
                  stop 'failed comparison test in solve_abtilu_Bi_CG_Stab'
               end if
               if (verbose) write(*,*) 'checked solve_abtilu'
            end if
         end subroutine solve_abtilu_Bi_CG_Stab
         
         logical function check_match(v1, v2) result(okay)
            real(dp), intent(in) :: v1(:), v2(:) ! (neq)
            integer :: j
            include 'formats'
            okay = .true.
            do j=1,neq
               if (abs(v1(j) - v2(j)) > 1d-12*(1d-50 + abs(v1(j)))) then
                  okay = .false.
                  write(*,2) 'v1 v2 bad', j, &
                     (v1(j) - v2(j))/(1d-50 + abs(v1(j))), &
                     v1(j), v2(j)
               end if
            end do
         end function check_match            

         subroutine matvec_abtilu_Bi_CG_Stab(b1, r1) ! set r = Jacobian*b
            real(dp), intent(in) :: b1(:) ! (neq)
            real(dp), intent(out) :: r1(:) ! (neq)
            integer :: k
            include 'formats'
            call block_tridiag_mv1(nvar, nz, lblk, dblk, ublk, b1, r1)
         end subroutine matvec_abtilu_Bi_CG_Stab
          
         subroutine matvec(x, r) ! set r = A*x
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: r(:)
            r = matmul(A,x)
         end subroutine matvec
         
         subroutine psolve(x,z) ! set z = Precond*x
            use const_def, only: dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: z(:)
            integer :: ipiv(neq), ierr
            real(dp) :: Acopy(neq,neq)
            z = x
            if (use_dummy_psolve) return
            ierr = 0
            Acopy = A
            call DGESV( neq, 1, Acopy, neq, ipiv, z, neq, ierr )
            if (ierr /= 0) then
               write(*,*) 'DGESV failed', ierr
               stop 'test_Bi_CG_Stab'
            end if
         end subroutine psolve
         
      end subroutine solve_abtilu_with_Bi_CG_Stab
      
      
      subroutine solve_abtilu_with_mgmres( &
            nvar, nz, ublk, dblk, lblk, rhs1, &
            equilibrate, exact, &
            num_sweeps_factor, num_sweeps_solve, & ! for abtilu
            itr_max, mr, tol_abs, tol_rel, & ! for mgmres
            soln1, verbose, ierr)
         integer, intent(in) :: nvar, nz, itr_max, mr, &
            num_sweeps_factor, num_sweeps_solve
         real(dp), dimension(:,:,:), intent(inout) :: & !(nvar,nvar,nz)
            ublk, dblk, lblk
         real(dp), dimension(:), intent(inout) :: rhs1 ! (neq)
         logical, intent(in) :: equilibrate, exact, verbose
         real(dp), intent(in) :: tol_abs, tol_rel
         real(dp), dimension(:), intent(inout) :: soln1 ! (neq)
            ! input: initial guess (can be 0)
            ! output: final approximation
         integer, intent(out) :: ierr
         
         real(dp), dimension(:,:,:), allocatable :: & !(nvar,nvar,nz)
            Dhat, invDhat, invDhat_lblk, invDhat_ublk
         real(dp), dimension(:), allocatable :: & ! (neq)
            prev1, w1, x1, y1, r, DR, DC
         real(dp), allocatable :: v(:,:) ! (neq,mr+1)
         integer, allocatable :: ipiv(:,:) ! (nvar,nz)
         integer :: neq, j
         include 'formats'
         ierr = 0
         neq = nvar*nz

         allocate( &
            Dhat(nvar,nvar,nz), invDhat(nvar,nvar,nz), &
            invDhat_lblk(nvar,nvar,nz), invDhat_ublk(nvar,nvar,nz), &
            prev1(neq), w1(neq), x1(neq), y1(neq), r(neq), &
            DR(neq), DC(neq), v(neq,mr+1), ipiv(nvar,nz))
         call create_preconditioner_abtilu(ierr)     
         if (ierr /= 0) stop 'failed in create_preconditioner_abtilu_mgmres'
         
         if (equilibrate) then ! scale rhs1 with DR
            if (verbose .and. nz < 20) then
               write(*,*) 'scale rhs1'
               call show_vec(nvar, nz, rhs1)
            end if
            !$omp simd
            do j=1,neq
               rhs1(j) = rhs1(j)*DR(j)
            end do
            if (verbose .and. nz < 20) call show_vec(nvar, nz, rhs1)
         end if
         
         call mgmres( &     
            neq, matvec_abtilu_mgmres, solve_abtilu_mgmres, &
            soln1, rhs1, r, v, itr_max, mr, &
            tol_abs, tol_rel, verbose )
         
         if (equilibrate) then ! scale soln1 with DC
            !$omp simd
            do j=1,neq
               soln1(j) = soln1(j)*DC(j)
            end do
         end if
         
         if (verbose .and. nz < 20) then
            write(*,*)
            write(*,*) 'solve_abtilu_with_mgmres solution'
            call show_vec(nvar, nz, soln1)
         end if
         
         contains
         
         subroutine create_preconditioner_abtilu(ierr)
            integer, intent(out) :: ierr
            include 'formats'
            call equilibrate_and_factor_abtilu_left( &
                nvar, nz, num_sweeps_factor, equilibrate, exact, & ! input
                lblk, dblk, ublk, & ! input/output
                Dhat, ipiv, & ! work
                DR, DC, invDhat, invDhat_lblk, invDhat_ublk, &
                verbose .and. nz < 20, ierr) ! output
         end subroutine create_preconditioner_abtilu
     
         subroutine solve_abtilu_mgmres(v1)
            real(dp), intent(inout) :: v1(:) ! (neq)
            include 'formats'
            call solve_abtilu_left( &
                nvar, nz, num_sweeps_solve, &
                invDhat, invDhat_lblk, invDhat_ublk, lblk, ublk, v1, & ! input
                DR, DC, equilibrate, exact, & ! input
                prev1, w1, x1, y1, & ! work
                v1, verbose, ierr) ! output
         end subroutine solve_abtilu_mgmres

         subroutine matvec_abtilu_mgmres(b1, r1) ! set r = Jacobian*b
            real(dp), intent(in) :: b1(:) ! (neq)
            real(dp), intent(out) :: r1(:) ! (neq)
            integer :: k
            include 'formats'
            call block_tridiag_mv1(nvar, nz, lblk, dblk, ublk, b1, r1)
         end subroutine matvec_abtilu_mgmres
         
      end subroutine solve_abtilu_with_mgmres
         
         
!*****************************************************************************
!*****************************************************************************
!
!  support routines for abtilu
!
!*****************************************************************************
!*****************************************************************************


      subroutine solve_abtilu_right( &
            nvar, nz, num_sweeps, & ! input
            invDhat, invDhat_lblk, invDhat_ublk, lblk, ublk, &
            b1, DR, DC, equilibrate, exact, & ! input
            prev1, w1, x1, y1, & ! work
            z1, verbose, ierr) ! output
         ! kashi phd thesis, pg 122, redone for block tridiagonal
         integer, intent(in) :: nvar, nz, num_sweeps
         logical, intent(in) :: equilibrate, exact, verbose
         real(dp), dimension(:,:,:), intent(in) :: & ! input (nvar,nvar,nz)
            invDhat, invDhat_lblk, invDhat_ublk, lblk, ublk
         real(dp), dimension(:), intent(in) :: b1, DR, DC ! input (nvar*nz)
         real(dp), dimension(:), intent(out) :: & ! (nvar*nz)
            prev1, w1, x1, y1, & ! work
            z1 ! output
         integer, intent(out) :: ierr
         ierr = 0
         call solve_abtilu_left( &
            nvar, nz, num_sweeps, & ! input
            invDhat, invDhat_lblk, invDhat_ublk, lblk, ublk, &
            b1, DR, DC, equilibrate, exact, & ! input
            prev1, w1, x1, y1, & ! work
            z1, verbose, ierr) ! output
      end subroutine solve_abtilu_right
      

      subroutine equilibrate_and_factor_abtilu_right( &
            nvar, nz, num_sweeps, equilibrate, exact, & ! input
            lblk, dblk, ublk, & ! input/output
            Dhat, ipiv, & ! work
            DR, DC, invDhat, invDhat_lblk, invDhat_ublk, &
            verbose, ierr) ! output
         integer, intent(in) :: nvar, nz, num_sweeps
         logical, intent(in) :: equilibrate, exact, verbose
         real(dp), dimension(:,:,:), intent(out) :: Dhat ! work (nvar,nvar,nz)
         integer, intent(out) :: ipiv(:,:) ! work (nvar,nz)
         real(dp), dimension(:,:,:), intent(inout) :: & ! input (nvar,nvar,nz)
            lblk, dblk, ublk
         real(dp), dimension(:), intent(out) :: DR, DC ! output (nvar*nz)
         real(dp), dimension(:,:,:), intent(out) :: & ! output (nvar,nvar,nz)
            invDhat, invDhat_lblk, invDhat_ublk
         integer, intent(out) :: ierr
         ierr = 0
         call equilibrate_and_factor_abtilu_left( &
            nvar, nz, num_sweeps, equilibrate, exact, & ! input
            lblk, dblk, ublk, & ! input/output
            Dhat, ipiv, & ! work
            DR, DC, invDhat, invDhat_lblk, invDhat_ublk, &
            verbose, ierr)
      end subroutine equilibrate_and_factor_abtilu_right
      
      
      subroutine solve_abtilu_left( &
            nvar, nz, num_sweeps, & ! input
            invDhat, invDhat_lblk, invDhat_ublk, lblk, ublk, &
            b1, DR, DC, equilibrate, exact, & ! input
            prev1, w1, x1, y1, & ! work
            z1, verbose, ierr) ! output
         ! kashi phd thesis, pg 122, redone for block tridiagonal
         integer, intent(in) :: nvar, nz, num_sweeps
         logical, intent(in) :: equilibrate, exact, verbose
         real(dp), dimension(:,:,:), intent(in) :: & ! input (nvar,nvar,nz)
            invDhat, invDhat_lblk, invDhat_ublk, lblk, ublk
         real(dp), dimension(:), intent(in) :: b1, DR, DC ! input (nvar*nz)
         real(dp), dimension(:), intent(out) :: & ! (nvar*nz)
            prev1, w1, x1, y1, & ! work
            z1 ! output
         integer, intent(out) :: ierr
         logical :: incomplete
         logical, parameter :: debug = .false.
         integer :: swp, j, k, neq
         include 'formats'
         ierr = 0
         incomplete = .not. exact
         neq = nvar*nz
         
         !$OMP PARALLEL DO PRIVATE(k)
         do k=1,nz ! initialize y(k) = invDhat(k)*b(k)
            call set_y_to_inDhat_b(k)
            call copy_y_to_w(k) ! save invDhat*b in w
         end do
         !$OMP END PARALLEL DO  
         
         do swp=1,num_sweeps
            ! analogous to forward elimination phase for tridiagonal solve
            !$omp simd
            do j=1,neq
               prev1(j) = y1(j) ! initialize prev to y
            end do
            !$OMP PARALLEL DO PRIVATE(k) IF(incomplete)
            do k=1,nz
               call set_y(k)
            end do
            !$OMP END PARALLEL DO  
            if (exact) exit             
         end do
         
         !$omp simd
         do j=1,neq
            z1(j) = y1(j) ! initialize z to y
         end do
         
         do swp=1,num_sweeps               
            ! analogous to backward substitution phase for tridiagonal solve
            !$omp simd
            do j=1,neq
               prev1(j) = z1(j) ! initialize prev to z
            end do
            !$OMP PARALLEL DO PRIVATE(k) IF(incomplete)
            do k=nz,1,-1 ! ok parallel
               call set_z(k)
            end do
            !$OMP END PARALLEL DO   
            if (exact) exit             
         end do
         
         contains
                  
         subroutine set_y_to_inDhat_b(k) ! y(k) = invDhat(k)*b(k)
            integer, intent(in) :: k
            integer :: s00
            s00 = (k-1)*nvar 
            call mv_0(invDhat(:,:,k),b1(s00+1:s00+nvar),y1(s00+1:s00+nvar))
         end subroutine set_y_to_inDhat_b
         
         subroutine copy_y_to_w(k)
            integer, intent(in) :: k
            integer :: j, s00
            s00 = (k-1)*nvar 
            !$omp simd
            do j=1,nvar
               w1(s00+j) = y1(s00+j)
            end do
         end subroutine copy_y_to_w
                  
         subroutine set_y(k) 
            ! y(k) = invDhat(k)*(b(k) - lblk(k)*y(k-1)) eq 5.9
            ! y(k) = invDhat(k)*b(k) - invDhat(k)*lblk(k)*y(k-1)
            ! y(k) = w(k) - invDhat_lblk(k)*y(k-1)
            integer, intent(in) :: k
            integer :: s00, sm1, j
            s00 = (k-1)*nvar 
            sm1 = s00 - nvar
            if (debug) then ! for debug, do not use invDhat_lblk.  use z as temp instead.
               if (k == 1) then
                  !$omp simd
                  do j=1,nvar ! y = w = invDhat*b
                     y1(s00+j) = w1(s00+j)
                  end do
                  return
               end if
               ! z(k) = lblk(k)*prev1(k-1)
               call mv_0(lblk(:,:,k),prev1(sm1+1:sm1+nvar),z1(s00+1:s00+nvar))
               ! y(k) = invDhat(k)*z(k)
               call mv_0(invDhat(:,:,k),z1(s00+1:s00+nvar),y1(s00+1:s00+nvar))
               ! y(k) = w(j) - y(k)
               !$omp simd
               do j=1,nvar ! y = w = invDhat*b
                  y1(s00+j) = w1(s00+j) - y1(s00+j)
               end do
               return
            end if
            !$omp simd
            do j=1,nvar ! y = w = invDhat*b
               y1(s00+j) = w1(s00+j)
            end do
            if (k == 1) return
            call mv_minus( & ! y(k) = y(k) - a(k)*prev1(k-1), y different than x
               invDhat_lblk(:,:,k),prev1(sm1+1:sm1+nvar),y1(s00+1:s00+nvar))
         end subroutine set_y
         
         subroutine set_z(k) 
            ! z(k) = y(k) - invDhat(k)*ublk(k)*z(k+1) eq 5.10
            ! z(k) = y(k) - invDhat_ublk(k)*z(k+1)
            ! z(k) = y(k) - invDhat_ublk(k)*prev(k+1)
            integer, intent(in) :: k
            integer :: s00, sp1, j
            s00 = (k-1)*nvar
            sp1 = s00 + nvar
            if (debug) then ! for debug, do not use invDhat_ublk.  use w1 as temp instead.
               if (k == nz) then
                  !$omp simd
                  do j=1,nvar ! z = y
                     z1(s00+j) = y1(s00+j)
                  end do
                  return
               end if
               ! w(k) = ublk(k)*prev(k+1)
               call mv_0(ublk(:,:,k),prev1(sp1+1:sp1+nvar),w1(s00+1:s00+nvar))
               ! z(k) = invDhat(k)*w(k)
               call mv_0(invDhat(:,:,k),w1(s00+1:s00+nvar),z1(s00+1:s00+nvar))
               ! z(k) = y(k) - z(k)
               !$omp simd
               do j=1,nvar ! z = y - z
                  z1(s00+j) = y1(s00+j) - z1(s00+j)
               end do
               return
            end if
            !$omp simd
            do j=1,nvar ! z = y
               z1(s00+j) = y1(s00+j)
            end do
            if (k == nz) return
            call mv_minus( & ! z(k) = z(k) - a(k)*prev1(k+1)
               invDhat_ublk(:,:,k),prev1(sp1+1:sp1+nvar),z1(s00+1:s00+nvar))
         end subroutine set_z
         
      end subroutine solve_abtilu_left


      subroutine equilibrate_and_factor_abtilu_left( &
            nvar, nz, num_sweeps, equilibrate, exact, & ! input
            lblk, dblk, ublk, & ! input/output
            Dhat, ipiv, & ! work
            DR, DC, invDhat, invDhat_lblk, invDhat_ublk, &
            verbose, ierr) ! output
         integer, intent(in) :: nvar, nz, num_sweeps
         logical, intent(in) :: equilibrate, exact, verbose
         real(dp), dimension(:,:,:), intent(out) :: Dhat ! work (nvar,nvar,nz)
         integer, intent(out) :: ipiv(:,:) ! work (nvar,nz)
         real(dp), dimension(:,:,:), intent(inout) :: & ! input (nvar,nvar,nz)
            lblk, dblk, ublk
         real(dp), dimension(:), intent(out) :: DR, DC ! output (nvar*nz)
         real(dp), dimension(:,:,:), intent(out) :: & ! output (nvar,nvar,nz)
            invDhat, invDhat_lblk, invDhat_ublk
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         if (equilibrate) then
            if (verbose) call show_mtx(nvar, nz, lblk, dblk, ublk)
            call get_scaling_vectors( &
               nvar, nz, lblk, dblk, ublk, & ! input
               DR, DC, ierr) ! output
            if (ierr /= 0) stop 'failed in get_scaling_vectors'
            call apply_scaling_vectors( &
               nvar, nz, DR, DC, & ! input
               lblk, dblk, ublk, & ! input/output
               ierr)
            if (ierr /= 0) stop 'failed in apply_scaling_vectors'
            if (verbose) then
               write(*,*)
               call show_mtx(nvar, nz, lblk, dblk, ublk)
               write(*,*) 'DR'
               call show_vec(nvar, nz, DR)
               write(*,*) 'DC'
               call show_vec(nvar, nz, DC)
               write(*,*)
            end if
            !stop 'equilibrate_and_factor_abtilu'
         end if
         call factor_abtilu_left( &
            nvar, nz, lblk, dblk, ublk, & ! input
            num_sweeps, exact, & ! input
            Dhat, ipiv, & ! work
            invDhat, invDhat_lblk, invDhat_ublk, ierr) ! output
         if (ierr /= 0) stop 'failed in factor_abtilu'
      end subroutine equilibrate_and_factor_abtilu_left


      subroutine apply_scaling_vectors( &
            nvar, nz, DR, DC, & ! input
            lblk, dblk, ublk, & ! input/output
            ierr)
         integer, intent(in) :: nvar, nz
         real(dp), dimension(:), intent(in) :: DR, DC ! (nvar*nz)
         real(dp), dimension(:,:,:), intent(inout) :: & ! (nvar,nvar,nz)
            lblk, dblk, ublk
         integer, intent(out) :: ierr
         integer :: iter, k, i
         include 'formats'
         ierr = 0
         
         !$OMP PARALLEL DO PRIVATE(k,i)
         do k=1,nz
            do i=1,nvar
               call row_scale(i,k)
            end do
         end do
         !$OMP END PARALLEL DO
         
         !$OMP PARALLEL DO PRIVATE(k,i)
         do k=1,nz
            do i=1,nvar
               call col_scale(i,k)
            end do
         end do
         !$OMP END PARALLEL DO
            
         contains
         
         subroutine row_scale(i,k)
            integer, intent(in) :: i, k
            real(dp) :: scale
            integer :: j
            include 'formats'
            scale = DR((k-1)*nvar + i)
            if (k > 1) then
               do j=1,nvar
                  lblk(i,j,k) = scale*lblk(i,j,k)
               end do
            end if
            do j=1,nvar
               dblk(i,j,k) = scale*dblk(i,j,k)
            end do
            if (k < nz) then
               do j=1,nvar
                  ublk(i,j,k) = scale*ublk(i,j,k)
               end do
            end if
         end subroutine row_scale
         
         subroutine col_scale(j,k) ! including effect of DR
            integer, intent(in) :: j, k
            real(dp) :: scale
            integer :: i
            include 'formats'
            scale = DC((k-1)*nvar + j)
            if (k > 1) then
               do i=1,nvar
                  ublk(i,j,k-1) = ublk(i,j,k-1)*scale
               end do
            end if
            do i=1,nvar
               dblk(i,j,k) = dblk(i,j,k)*scale
            end do
            if (k < nz) then
               do i=1,nvar
                  lblk(i,j,k+1) = lblk(i,j,k+1)*scale
               end do
            end if
         end subroutine col_scale
            
      end subroutine apply_scaling_vectors


      subroutine get_scaling_vectors( & ! from lapack dgeequ
            nvar, nz, lblk, dblk, ublk, & ! input
            DR, DC, ierr) ! output
         integer, intent(in) :: nvar, nz
         real(dp), dimension(:,:,:), intent(in) :: & ! input (nvar,nvar,nz)
            lblk, dblk, ublk
         real(dp), dimension(:), intent(out) :: DR, DC ! output (nvar*nz)
         integer, intent(out) :: ierr
         integer :: k, i
         include 'formats'
         ierr = 0
        
         ! compute row scale factors
         !$OMP PARALLEL DO PRIVATE(k,i)
         do k=1,nz
            do i=1,nvar
               DR((k-1)*nvar+i) = 1d0/row_max_abs(i,k)
            end do
         end do
         !$OMP END PARALLEL DO
         
         ! compute column scale factors assuming effect of DR
         !$OMP PARALLEL DO PRIVATE(k,i)
         do k=1,nz
            do i=1,nvar
               DC((k-1)*nvar+i) = 1d0/col_max_abs(i,k)
            end do
         end do
         !$OMP END PARALLEL DO
            
         contains
         
         real(dp) function row_max_abs(i,k)
            integer, intent(in) :: i, k
            real(dp) :: lmax, dmax, umax
            integer :: j
            include 'formats'
            lmax = 0d0
            dmax = 0d0
            umax = 0d0
            if (k > 1) then
               do j=1,nvar
                  lmax = max(lmax, abs(lblk(i,j,k)))
                  !write(*,4) 'lblk lmax', i, j, k, lblk(i,j,k), lmax
               end do
            end if
            do j=1,nvar
               dmax = max(dmax, abs(dblk(i,j,k)))
               !write(*,4) 'dblk dmax', i, j, k, dblk(i,j,k), dmax
            end do
            if (k < nz) then
               do j=1,nvar
                  umax = max(umax, abs(ublk(i,j,k)))
                  !write(*,4) 'ublk umax', i, j, k, ublk(i,j,k), umax
               end do
            end if
            row_max_abs = max(lmax, dmax, umax)
            !write(*,3) 'row_max_abs', i, k, row_max_abs, lmax, dmax, umax
         end function row_max_abs
         
         real(dp) function col_max_abs(j,k) ! including effect of DR
            integer, intent(in) :: j, k
            real(dp) :: lmax, dmax, umax
            integer :: i, s
            include 'formats'
            lmax = 0d0
            dmax = 0d0
            umax = 0d0
            if (k > 1) then
               s = (k-2)*nvar
               do i=1,nvar
                  umax = max(umax, DR(s+i)*abs(ublk(i,j,k-1)))
               end do
            end if
            s = (k-1)*nvar
            do i=1,nvar
               dmax = max(dmax, DR(s+i)*abs(dblk(i,j,k)))
            end do
            if (k < nz) then
               s = k*nvar
               do i=1,nvar
                  lmax = max(lmax, DR(s+i)*abs(lblk(i,j,k+1)))
               end do
            end if
            col_max_abs = max(lmax, dmax, umax)
            !write(*,3) 'col_max_abs', i, k, col_max_abs, lmax, dmax, umax
         end function col_max_abs
            
      end subroutine get_scaling_vectors
      
      
      subroutine factor_abtilu_left( &
            nvar, nz, lblk, dblk, ublk, & ! input
            num_sweeps, exact, & ! input
            Dhat, ipiv, & ! work
            invDhat, invDhat_lblk, invDhat_ublk, ierr) ! output
         ! kashi phd thesis, pg 121, Algorithm 12, redone for block tridiagonal
         integer, intent(in) :: nvar, nz, num_sweeps
         real(dp), dimension(:,:,:), intent(in) :: & ! input (nvar,nvar,nz)
            lblk, dblk, ublk
         logical, intent(in) :: exact
         real(dp), dimension(:,:,:), intent(out) :: Dhat ! work (nvar,nvar,nz)
         integer, intent(out) :: ipiv(:,:) ! work (nvar,nz)
         real(dp), dimension(:,:,:), intent(out) :: & ! output (nvar,nvar,nz)
            invDhat, invDhat_lblk, invDhat_ublk
         integer, intent(out) :: ierr
         integer :: i, j, k, swp, op_err
         logical :: incomplete
         include 'formats'
         ierr = 0
         op_err = 0
         incomplete = .not. exact
         
         ! initialize Dhat and invDhat (1:nz)
         !$OMP PARALLEL DO PRIVATE(k,op_err)
         do k = 1, nz
            call copy_dblk_to_Dhat(k)
            call set_invDhat(k,op_err)
            if (op_err /= 0) ierr = op_err
         end do      
         !$OMP END PARALLEL DO       
         if (ierr /= 0) then
            write(*,2) 'factor_abtilu: initial m_inverse failed', k
            return
         end if  

         ! iteratively refine Dhat and invDhat (2:nz)
         do swp=1,num_sweeps 
            !$OMP PARALLEL DO PRIVATE(k,op_err) IF(incomplete)
            do k = 2, nz
               call set_Dhat(k)
               if (exact) then
                  call set_invDhat(k,op_err)
                  if (op_err /= 0) ierr = op_err
               end if
            end do
            !$OMP END PARALLEL DO 
            if (exact) exit
            !$OMP PARALLEL DO PRIVATE(k,op_err)
            do k = 2, nz  ! update invDhat for new Dhat
               call set_invDhat(k,op_err)
               if (op_err /= 0) ierr = op_err
            end do
            !$OMP END PARALLEL DO       
            if (ierr /= 0) then
               write(*,3) 'factor_abtilu: sweep m_inverse failed', swp, k
               return
            end if  
         end do         
         
         ! create matrix products that are useful for solve
         !$OMP PARALLEL DO PRIVATE(k)
         do k = 1, nz 
            if (k > 1) then
               call set_invDhat_lblk(k)
            else
               invDhat_lblk(:,:,k) = 0d0
            end if
            if (k < nz) then
               call set_invDhat_ublk(k)
            else
               invDhat_ublk(:,:,k) = 0d0
            end if
         end do
         !$OMP END PARALLEL DO  
         
         contains
         
         subroutine set_invDhat_ublk(k)
            integer, intent(in) :: k
            ! invDhat_ublk(k) = invDhat(k)*ublk(k)
            call mm_0(invDhat(:,:,k), ublk(:,:,k), invDhat_ublk(:,:,k)) ! c := a*b
         end subroutine set_invDhat_ublk
         
         subroutine set_invDhat_lblk(k)
            integer, intent(in) :: k
            ! invDhat_ublk(k) = invDhat(k)*lblk(k)
            call mm_0(invDhat(:,:,k), lblk(:,:,k), invDhat_lblk(:,:,k)) ! c := a*b
         end subroutine set_invDhat_lblk
         
         subroutine set_Dhat(k)
            integer, intent(in) :: k
            call copy_dblk_to_Dhat(k)
            if (k == 1) return
            call set_invDhat_ublk(k-1)
            ! Dhat(k) = dblk(k) - lblk(k)*invDhat_ublk(k-1)
            call mm_minus( &
               lblk(:,:,k), invDhat_ublk(:,:,k-1), Dhat(:,:,k)) ! c := c - a*b            
         end 
         
         subroutine copy_dblk_to_Dhat(k)
            integer, intent(in) :: k
            integer :: i, j
            do j=1,nvar
               !$omp simd
               do i=1,nvar
                  Dhat(i,j,k) = dblk(i,j,k)
               end do
            end do
         end subroutine copy_dblk_to_Dhat
         
         subroutine set_invDhat(k,ierr)
            integer, intent(in) :: k
            integer, intent(out) :: ierr
            include 'formats'
            ierr = 0
            call m_inverse(k, nvar, Dhat(:,:,k), invDhat(:,:,k), ipiv(:,k), ierr)
         end subroutine set_invDhat
                  
      end subroutine factor_abtilu_left

      
      subroutine show_vec(nvar, nz, v)
         integer, intent(in) :: nvar, nz
         real(dp), dimension(:), intent(in) :: v ! neq
         character (len=256) :: fmt
         integer :: i, k, neq
         include 'formats'
         fmt = '(1x,1pe11.4)'
         do k=1,nz
            do i=1,nvar
               write(*,fmt,advance='no') v((k-1)*nvar + i)
            end do
            write(*,'(2x)',advance='no') 
         end do
         write(*,*)
      end subroutine show_vec
      
      
      subroutine show_mtx(nvar, nz, lblk, dblk, ublk)
         integer, intent(in) :: nvar, nz
         real(dp), dimension(:,:,:), intent(in) :: lblk, dblk, ublk
         integer :: i, j, k, n
         character (len=256) :: fmt
         fmt = '(1x,1pe11.4)'
         include 'formats'
         write(*,*) 'show_mtx'
         do k=1,nz ! block row k
            do i=1,nvar ! row i of block k
               if (k > 2) then
                  do n=1,k-2  
                     do j=1,nvar
                        write(*,'(12x)',advance='no') 
                     end do
                     write(*,'(2x)',advance='no') 
                  end do
               end if
               if (k > 1) then ! show lblk row i of k
                  do j=1,nvar
                     write(*,fmt,advance='no') lblk(i,j,k)
                  end do
                  write(*,'(2x)',advance='no') 
               end if
               ! show dblk row i of k
               do j=1,nvar
                  write(*,fmt,advance='no') dblk(i,j,k)
               end do
               write(*,'(2x)',advance='no') 
               if (k < nz) then ! show ublk row i of k
                  do j=1,nvar
                     write(*,fmt,advance='no') ublk(i,j,k)
                  end do
               end if
               write(*,*) ! end of row j
            end do
            write(*,*)
         end do      
      end subroutine show_mtx

      
      subroutine copy_lower_to_square(nvar,nz,blk,s)
         ! k > 1
         ! ii = (k-2)*nvar + i
         ! jj = (k-1)*nvar + j
         ! s(ii,jj) = blk(i,j,k)
         integer, intent(in) :: nvar, nz
         real(dp), intent(in) :: blk(:,:,:) ! (nvar,nvar,nz)
         real(dp), intent(inout) :: s(:,:) ! (neq,neq)
         integer :: k, i, j, ii, jj
         do k=2,nz
            ii = (k-1)*nvar
            do i=1,nvar
               jj = (k-2)*nvar
               do j=1,nvar
                  s(ii+i,jj+j) = blk(i,j,k)
               end do
            end do
         end do
      end subroutine copy_lower_to_square

      subroutine copy_diag_to_square(nvar,nz,blk,s)         
         ! ii = (k-1)*nvar + i
         ! jj = (k-1)*nvar + j
         ! s(ii,jj) = blk(i,j,k)
         integer, intent(in) :: nvar, nz
         real(dp), intent(in) :: blk(:,:,:) ! (nvar,nvar,nz)
         real(dp), intent(inout) :: s(:,:) ! (neq,neq)
         integer :: k, i, j, ii, jj
         do k=1,nz
            ii = (k-1)*nvar
            do i=1,nvar
               jj = (k-1)*nvar
               do j=1,nvar
                  s(ii+i,jj+j) = blk(i,j,k)
               end do
            end do
         end do
      end subroutine copy_diag_to_square
      
      subroutine copy_upper_to_square(nvar,nz,blk,s)
         ! k < nz
         ! ii = k*nvar + i
         ! jj = (k-1)*nvar + j
         ! s(ii,jj) = blk(i,j,k)
         integer, intent(in) :: nvar, nz
         real(dp), intent(in) :: blk(:,:,:) ! (nvar,nvar,nz)
         real(dp), intent(inout) :: s(:,:) ! (neq,neq)
         integer :: k, i, j, ii, jj
         do k=1,nz
            ii = (k-1)*nvar
            do i=1,nvar
               jj = k*nvar
               do j=1,nvar
                  s(ii+i,jj+j) = blk(i,j,k)
               end do
            end do
         end do
      end subroutine copy_upper_to_square
      
      
!*****************************************************************************
!*****************************************************************************
!
!  matrix support routines
!
!*****************************************************************************
!*****************************************************************************
      
      
      subroutine block_tridiag_mv1(nvar, nz, lblk, dblk, ublk, b1, r1)
         ! set r = A*b with A = block tridiagonal given by lblk, dblk, ublk
         integer, intent(in) :: nvar, nz    
         real(dp), dimension(:,:,:), intent(in) :: lblk, dblk, ublk ! (nvar,nvar,nz)
         real(dp), dimension(:), intent(in) :: b1 ! (nvar*nz)
         real(dp), dimension(:), intent(out) :: r1 ! (nvar*nz)         
         integer :: k,s
         !$OMP PARALLEL DO PRIVATE(k,s)
         do k = 1, nz
            ! prod(k) = dblk(k)*b(k)
            s = (k-1)*nvar
            call mv_0(dblk(:,:,k), b1(s+1:s+nvar), r1(s+1:s+nvar))
            if (k > 1) then
               ! prod(k) = prod(k) + lblk(k)*b(k-1)
               call mv_plus(lblk(:,:,k), b1(s+1-nvar:s), r1(s+1:s+nvar))
            end if
            if (k < nz) then
               ! prod(k) = prod(k) + ublk(k)*b(k+1)
               call mv_plus(ublk(:,:,k), b1(s+1+nvar:s+2*nvar), r1(s+1:s+nvar))
            end if
         end do      
         !$OMP END PARALLEL DO         
      end subroutine block_tridiag_mv1                  
      
      subroutine m_inverse(k, nvar, blk, blk_inv, ipiv, ierr)
         use star_bcyclic, only: my_getf2
         integer, intent(in) :: k, nvar
         real(dp), intent(in) :: blk(:,:) ! (nvar,nvar)
         real(dp), intent(out) :: blk_inv(:,:) ! (nvar,nvar)
         integer, intent(out) :: ipiv(:), ierr
         integer :: i, j
         include 'formats'
         do j=1,nvar
            !$omp simd
            do i=1,nvar
               blk_inv(i,j) = blk(i,j)
            end do
         end do
         call my_getf2(nvar, blk_inv, nvar, ipiv, ierr)
         if (ierr /= 0) then
            write(*,3) 'my_getf2 failed', k, ierr
            stop 'm_inverse'
         end if
         call my_getri(nvar, blk_inv, ipiv, ierr)
         if (ierr /= 0) then
            write(*,3) 'my_getri failed', k, ierr
            stop 'm_inverse'
         end if         
      end subroutine m_inverse
      
      subroutine my_getri(nvar, blk_inv, ipiv, ierr)
         integer, intent(in) :: nvar
         real(dp), intent(inout) :: blk_inv(:,:) ! (nvar,nvar)
         integer, intent(in) :: ipiv(:)
         integer, intent(out) :: ierr
         real(dp) :: work(nvar), binv(nvar,nvar)
         integer :: i, j, ip(nvar)
         !$omp simd
         do j=1,nvar
            ip(j) = ipiv(j)
         end do
         do j=1,nvar
            !$omp simd
            do i=1,nvar
               binv(i,j) = blk_inv(i,j)
            end do
         end do
         call DGETRI(nvar, binv, nvar, ip, work, nvar, ierr)
         do j=1,nvar
            !$omp simd
            do i=1,nvar
               blk_inv(i,j) = binv(i,j)
            end do
         end do
      end subroutine my_getri
         
      subroutine mm_0(a, b, c) ! c := a*b, c different than b
         real(dp), dimension(:,:) :: a, b, c ! (nvar,nvar)
         integer :: j, i, nvar
         include 'formats'
         nvar=size(a,dim=1)
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
         integer :: j, l, i, nvar
         nvar=size(a,dim=1)
         do j = 1,nvar
            do l = 1,nvar
               tmp = b(l,j)
               if (tmp /= 0d0) then
                  !$omp simd
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
         integer :: j, l, i, nvar
         nvar=size(a,dim=1)
         do j = 1,nvar
            do l = 1,nvar
               tmp = b(l,j)
               if (tmp /= 0d0) then
                  !$omp simd
                  do i = 1,nvar
                     c(i,j) = c(i,j) - tmp*a(i,l)
                  end do
               end if
            end do
         end do      
      end subroutine mm_minus

      subroutine mv_self(a,x,wrk) ! x = a*x
         real(dp) :: a(:,:) ! (nvar,nvar)
         real(dp) :: x(:), wrk(:) ! (nvar)
         integer :: j, nvar
         nvar=size(a,dim=1)
         !$omp simd
         do j = 1,nvar
            wrk(j) = x(j)
         end do
         call mv_0(a,wrk,x) 
      end subroutine mv_self

      subroutine mv_0(a,x,y) ! y = a*x, y different than x
         real(dp) :: a(:,:) ! (nvar,nvar)
         real(dp) :: x(:), y(:) ! (nvar)
         integer :: j, nvar
         nvar=size(a,dim=1)
         do j = 1,nvar
            y(j) = 0d0
         end do
         call mv_plus(a,x,y)
      end subroutine mv_0

      subroutine mv_plus(a,x,y) ! y = y + a*x, y different than x
         real(dp) :: a(:,:) ! (nvar,nvar)
         real(dp) :: x(:), y(:) ! (nvar)
         real(dp) :: tmp
         integer :: j, i, nvar
         nvar=size(a,dim=1)
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
         integer :: j, i, nvar
         nvar=size(a,dim=1)
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


!*****************************************************************************
!*****************************************************************************
!
!  MGMRES
!
!*****************************************************************************
!*****************************************************************************

      !*****************************************************************************
      !
      !  MGMRES applies restarted GMRES.   derived from Burkardt's MGMRES_ST
      !
      !  Discussion:
      !
      !    The linear system A*X=B is solved iteratively.
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
      subroutine mgmres( &
            n, matvec, psolve, x, rhs, &
            r, v, & ! work
            itr_max, mr, tol_abs, tol_rel, verbose )
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
        real(dp), intent(inout) :: x(:) ! (n)
        real(dp), intent(in) :: rhs(:) ! (n)
        real(dp), intent(out) :: r(:) ! (n)
        real(dp), intent(out) :: v(:,:) ! (n,mr+1)
        integer, intent(in) :: itr_max, mr
        real(dp), intent(in) :: tol_abs, tol_rel
        logical, intent(in) :: verbose
         
        ! locals
        real(dp) :: av
        real(dp) :: c(mr+1)
        real(dp), parameter :: delta = 1.0D-03
        real(dp) :: g(mr+1)
        real(dp) :: h(mr+1,mr)
        real(dp) :: htmp
        integer :: i
        integer :: itr
        integer :: itr_used
        integer :: j
        integer :: k
        integer :: k_copy
        real(dp) :: mu
        real(dp) :: rho
        real(dp) :: rho_tol
        real(dp) :: s(mr+1)
        real(dp) :: y(mr+1)

        itr_used = 0
        if ( verbose ) then
          !write ( *, '(a,i4)' ) '  mgmres number of unknowns = ', n
        end if
        k_copy = 0
        do itr = 1, itr_max
          call matvec ( x, r )
          r(1:n) = rhs(1:n) - r(1:n)
          call psolve ( r ) ! apply to residual
          rho = sqrt ( dot_product ( r, r ) )
          if (rho == 0d0) exit ! can happen when form exact LU and exact apply
          if ( itr == 1 ) then
            rho_tol = rho * tol_rel
          end if
          if ( verbose ) then
            write ( *, '(a,i4,a,3g14.6)' ) &
               '  ITR = ', itr, '  Residual = ', rho, rho_tol, tol_abs
          end if
          v(1:n,1) = r(1:n) / rho
          g(1) = rho
          g(2:mr+1) = 0.0D+00
          h(1:mr+1,1:mr) = 0.0D+00
          do k = 1, mr
            k_copy = k
            call matvec ( v(1:n,k), v(1:n,k+1) )
            call psolve ( v(1:n,k+1) ) ! apply to result of matvec
            av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
            if (is_bad(av)) then
               write(*,*) 'bad result from matvec/psolve: is_bad(av)', av
               stop 'mgmres'
            end if
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
              write ( *, '(a,i4,a,g14.6)' ) '    K = ', k, '  Residual = ', rho
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

        contains
  
        subroutine mult_givens ( c, s, k, g )
          integer :: k
          real(dp) :: c
          real(dp) :: g(1:k+1)
          real(dp) :: g1
          real(dp) :: g2
          real(dp) :: s

          g1 = c * g(k) - s * g(k+1)
          g2 = s * g(k) + c * g(k+1)

          g(k)   = g1
          g(k+1) = g2

          return
        end subroutine mult_givens
        
      end subroutine mgmres
      
      
!*****************************************************************************
!*****************************************************************************
!
!  testing abtilu
!
!*****************************************************************************
!*****************************************************************************

      
      subroutine write_MM_mxt( &
            nvar, nz, ublk, dblk, lblk, filename)
         integer, intent(in) :: nvar, nz
         real(dp), dimension(:,:,:), intent(in) :: & ! (nvar,nvar,nz)
            ublk, dblk, lblk
         character (len=*), intent(in) :: filename
         integer :: neq, non_zeros, iounit
         neq = nvar*nz
         non_zeros = 0d0
         call for_each_nonzero(.false.)
         
         open(newunit=iounit, file=trim(filename), action='write', status='replace')
         write(iounit, '(a)') '%%MatrixMarket matrix coordinate double general'
         write(iounit, '(a)') '%from mesa/star abtilu'
         write(iounit, '(3i8)') neq, neq, non_zeros
         call for_each_nonzero(.true.)
         close(iounit)
         
         contains
         
         subroutine for_each_nonzero(write_flag)
            logical, intent(in) :: write_flag
            integer :: i, j, k, r, c
            do k=1,nz
               do i=1,nvar 
                  r = (k-1)*nvar + i
                  if (k > 1) then
                     c = (k-2)*nvar
                     do j=1,nvar
                        call do1(r, c+j, lblk(i,j,k), write_flag)
                     end do
                  end if
                  c = (k-1)*nvar
                  do j=1,nvar
                     call do1(r, c+j, dblk(i,j,k), write_flag)
                  end do
                  if (k < nz) then
                     c = k*nvar
                     do j=1,nvar
                        call do1(r, c+j, ublk(i,j,k), write_flag)
                     end do
                  end if
               end do
            end do
         end subroutine for_each_nonzero
         
         subroutine do1(r,c,v,write_flag)
            integer, intent(in) :: r, c
            real(dp), intent(in) :: v
            logical, intent(in) :: write_flag
            if (v == 0d0) return
            if (write_flag) then
               write(iounit, '(2i8,1e26.16)') r, c, v
            else
               non_zeros = non_zeros + 1
            end if
         end subroutine do1
         
      end subroutine write_MM_mxt

      
      subroutine write_MM_vec(nvar, nz, b1, filename)
         integer, intent(in) :: nvar, nz
         real(dp), dimension(:), intent(in) :: b1 ! (nvar*nz)
         character (len=*), intent(in) :: filename
         integer :: neq, non_zeros, iounit
         neq = nvar*nz
         non_zeros = 0d0
         call for_each_nonzero(.false.)
         
         open(newunit=iounit, file=trim(filename), action='write', status='replace')
         write(iounit, '(a)') '%%MatrixMarket matrix coordinate double general'
         write(iounit, '(a)') '%from mesa/star abtilu'
         write(iounit, '(3i8)') 1, neq, non_zeros
         call for_each_nonzero(.true.)
         close(iounit)
         
         contains
         
         subroutine for_each_nonzero(write_flag)
            logical, intent(in) :: write_flag
            integer :: j, k, c
            do k=1,nz
               if (k > 1) then
                  c = (k-2)*nvar
                  do j=1,nvar
                     call do1(c+j, write_flag)
                  end do
               end if
               c = (k-1)*nvar
               do j=1,nvar
                  call do1(c+j, write_flag)
               end do
               if (k < nz) then
                  c = k*nvar
                  do j=1,nvar
                     call do1(c+j, write_flag)
                  end do
               end if
            end do
         end subroutine for_each_nonzero
         
         subroutine do1(c,write_flag)
            integer, intent(in) :: c
            logical, intent(in) :: write_flag
            real(dp) :: v
            v = b1(c)
            if (v == 0d0) return
            if (write_flag) then
               write(iounit, '(2i8,1e26.16)') 1, c, v
            else
               non_zeros = non_zeros + 1
            end if
         end subroutine do1
         
      end subroutine write_MM_vec
      
      
      subroutine get_lapack_solution( &
            neq, A, rhs1, actual_soln1, verbose, ierr)
         integer, intent(in) :: neq
         real(dp), intent(in) :: A(:,:)
         real(dp), intent(inout) :: rhs1(:), actual_soln1(:)
         logical, intent(in) :: verbose
         integer, intent(out) :: ierr
          character (len=1) :: fact, trans, equed
          real(dp) :: rcond
          integer, parameter :: nrhs = 1
          real(dp) :: ain(neq,neq), af(neq,neq), b(neq,nrhs), x(neq,nrhs), &
             r(neq), c(neq), ferr(nrhs), berr(nrhs), work(4*neq)
          integer :: ipiv(neq), iwork(neq), i, j
          fact = 'E' ! equilibrated and factored
          trans = 'N' ! no transpose
          equed = 'B' ! both row and column
          ierr = 0
          ain(1:neq,1:neq) = A(1:neq,1:neq)
          b(1:neq,1) = rhs1(1:neq)
          if (verbose) then
            write(*,*) 'A'
            do i=1,neq
               do j=1,neq
                  write(*,'(1x,f11.5)',advance='no') ain(i,j)
               end do
               write(*,*)
            end do
          end if
          call DGESVX(fact, trans, neq, nrhs, ain, neq, af, neq, ipiv, &
                      equed, r, c, b, neq, x, neq, rcond, ferr, berr, &
                      work, iwork, ierr)
         if (ierr /= 0) return
         actual_soln1(1:neq) = x(1:neq,1)
         if (verbose) then
            write(*,*) 'DGESVX equilibrated A'
            do i=1,neq
               do j=1,neq
                  write(*,'(1x,f11.5)',advance='no') r(i)*ain(i,j)*c(j)
               end do
               write(*,*)
            end do
            write(*,*) 'DGESVX soln1'
            do j=1,neq
               write(*,'(1x,f26.16)',advance='no') actual_soln1(j)
            end do
            write(*,*)
            write(*,*) 'DGESVX r'
            do j=1,neq
               write(*,'(1x,f26.16)',advance='no') r(j)
            end do
            write(*,*)
            write(*,*) 'DGESVX c'
            do j=1,neq
               write(*,'(1x,f26.16)',advance='no') c(j)
            end do
            write(*,*)
            write(*,*) 'DGESVX rcond', rcond
            write(*,*)
         end if
      end subroutine get_lapack_solution
      

      subroutine test_abtilu_nvar1(round)
         use utils_lib, only: fill_with_NaNs, fill_with_NaNs_2D, fill_with_NaNs_3D
         integer, intent(in) :: round
         integer, parameter :: nvar = 1, nz = 4, neq = nvar*nz
         real(dp), dimension(:,:), allocatable :: A
         real(dp), dimension(:,:,:), allocatable :: &
            ublk, dblk, lblk
         real(dp), dimension(:), allocatable :: &
            soln1, rhs1, actual_soln1
         integer :: test, i, k, mr, itr_max, shft, ierr, &
            num_sweeps_factor, num_sweeps_solve, iter
         logical :: equilibrate, exact, verbose, write_mm, &
            use_dummy_psolve, use_Bi_CG_Stab
         real(dp) :: tol_abs, tol_rel, &
            x_error, soln_error, error
         include 'formats'
         
         mr = neq - 1
         allocate(A(neq,neq), &
            ublk(nvar,nvar,nz), dblk(nvar,nvar,nz), lblk(nvar,nvar,nz), &
            soln1(neq), rhs1(neq), actual_soln1(neq))
         
         call fill_with_NaNs_2D(A)
         call fill_with_NaNs(soln1)
         call fill_with_NaNs(rhs1)
         call fill_with_NaNs(actual_soln1)
         call fill_with_NaNs_3D(ublk)
         call fill_with_NaNs_3D(dblk)
         call fill_with_NaNs_3D(lblk)
         
         A(1,:) = (/  5d0, -3d0,  0d0,  0d0 /)
         A(2,:) = (/ -1d0,  4d0,  2d0,  0d0 /)
         A(3,:) = (/  0d0, -5d0,  2d0, -2d0 /)
         A(4,:) = (/  0d0,  0d0, -1d0,  6d0  /)
         
         ublk(1,1,1) = A(1,2) 
         ublk(1,1,2) = A(2,3)
         ublk(1,1,3) = A(3,4)
         ublk(1,1,4) = 0d0
         
         dblk(1,1,1) = A(1,1)
         dblk(1,1,2) = A(2,2)
         dblk(1,1,3) = A(3,3)
         dblk(1,1,4) = A(4,4)
         
         lblk(1,1,1) = 0d0
         lblk(1,1,2) = A(2,1)
         lblk(1,1,3) = A(3,2)
         lblk(1,1,4) = A(4,3)
             
         itr_max = 50
         tol_abs = 1d-8
         tol_rel = 1d-8

         num_sweeps_factor = 1
         num_sweeps_solve = 3
         
         !exact = .true.
         exact = .false.
         
         verbose = .false.
         !verbose = .true.
         
         !equilibrate = .false.
         equilibrate = .true.
         
         write_mm = .false.
         !write_mm = .true.

         use_Bi_CG_Stab = .false.
         !use_Bi_CG_Stab = .true.

         rhs1 = 1d0

         if (write_mm) then
            call write_MM_mxt( &
               nvar, nz, ublk, dblk, lblk, 'abtilu_1_mtx.txt')
            call write_MM_vec(nvar, nz, rhs1, 'abtilu_1_rhs.txt')
         end if
                  
         use_dummy_psolve = .false.
         call psolve(rhs1,actual_soln1)
         rhs1 = 1d0
         soln1 = 0d0
         x_error = norm2_of_diff(neq, actual_soln1, soln1)
         
         if (.false.) then  ! this works
            !use_dummy_psolve = .false.
            use_dummy_psolve = .true.
            call Bi_CG_Stab( &
              matvec, psolve, rhs1, soln1, itr_max, &
              tol_rel, error, iter, verbose, ierr)
            if (ierr /= 0) then
               write(*,*) 'Bi_CG_Stab ierr', ierr
               return
            end if
         else if (use_Bi_CG_Stab) then  ! this fails
            write(*,*) 'call solve_abtilu_with_Bi_CG_Stab'
            call solve_abtilu_with_Bi_CG_Stab( &
               nvar, nz, A, ublk, dblk, lblk, rhs1, &
               equilibrate, exact, &
               num_sweeps_factor, num_sweeps_solve, &
               itr_max, tol_rel, &
               soln1, verbose, ierr)
         else  ! this works
            !write(*,*) 'call solve_abtilu_with_mgmres'
            call solve_abtilu_with_mgmres( &
               nvar, nz, ublk, dblk, lblk, rhs1, &
               equilibrate, exact, &
               num_sweeps_factor, num_sweeps_solve, &
               itr_max, mr, tol_abs, tol_rel, &
               soln1, verbose, ierr)
         end if

         if (write_mm) then
            call write_MM_vec(nvar, nz, soln1, 'abtilu_1_soln.txt')
            call write_MM_vec(nvar, nz, actual_soln1, 'abtilu_1_DGESVX_soln.txt')
         end if
         
         if (verbose) then
            write(*,*) 'actual and abtilu solutions'
            call show_vec(nvar, nz, actual_soln1)
            call show_vec(nvar, nz, soln1)
         end if
         x_error = norm2_of_diff(neq, actual_soln1, soln1)
         if (x_error > 1d-12) then
            write ( *, '(a,i4,g19.11)' ) '  test_abtilu_nvar1 x=soln error ', round, x_error
            stop 'bad x_error test_abtilu_nvar1'
         else if (use_Bi_CG_Stab) then
            write ( *, '(a,i4,g19.11)' ) '  Bi_CG_Stab test_abtilu_nvar1 x=soln good ', round, x_error
         else
            write ( *, '(a,i4,g19.11)' ) '  mgmres test_abtilu_nvar1 x=soln good ', round, x_error
         end if

        contains
        
        real(dp) function norm2_of_diff(neq,a,b) result(v)
           integer, intent(in) :: neq
           real(dp), intent(in), dimension(:) :: a, b
           v = sqrt(sum(pow2(a(1:neq) - b(1:neq))))
        end function norm2_of_diff
          
         subroutine matvec(x, r) ! set r = A*x
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: r(:)
            r = matmul(A,x)
         end subroutine matvec
         
         subroutine psolve(x,z) ! set z = Precond*x
            use const_def, only: dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: z(:)
            integer :: ipiv(neq), ierr
            real(dp) :: Acopy(neq,neq)
            z = x
            if (use_dummy_psolve) return
            ierr = 0
            Acopy = A
            call DGESV( neq, 1, Acopy, neq, ipiv, z, neq, ierr )
            if (ierr /= 0) then
               write(*,*) 'DGESV failed', ierr
               stop 'test_Bi_CG_Stab'
            end if
         end subroutine psolve
        
      end subroutine test_abtilu_nvar1      
      

      subroutine test_abtilu_nvar2(round)
         use utils_lib, only: fill_with_NaNs, fill_with_NaNs_2D, fill_with_NaNs_3D
         integer, intent(in) :: round
         integer, parameter :: nvar = 2, nz = 4, neq = nvar*nz
         real(dp), dimension(:,:), allocatable :: A(:,:)
         real(dp), dimension(:,:,:), allocatable :: &
            ublk, dblk, lblk
         real(dp), dimension(:), allocatable :: &
            soln1, rhs1, actual_soln1
         integer :: test, i, k, mr, itr_max, shft, ierr, &
            num_sweeps_factor, num_sweeps_solve
         logical :: equilibrate, exact, verbose, write_mm, use_Bi_CG_Stab
         real(dp) :: tol_abs, tol_rel, &
            x_error, soln_error
         include 'formats'
         
         mr = neq - 1
         allocate(A(neq,neq), &
            ublk(nvar,nvar,nz), dblk(nvar,nvar,nz), lblk(nvar,nvar,nz), &
            soln1(neq), rhs1(neq), actual_soln1(neq))
         
         call fill_with_NaNs(soln1)
         call fill_with_NaNs(rhs1)
         call fill_with_NaNs(actual_soln1)
         call fill_with_NaNs_2D(A)
         call fill_with_NaNs_3D(ublk)
         call fill_with_NaNs_3D(dblk)
         call fill_with_NaNs_3D(lblk)

         A(1,:) = (/  2d0,  3d0,  0d0, -1d0,  0d0,  0d0,  0d0,  0d0 /)
         A(2,:) = (/  0d0,  2d0, -1d0,  0d0,  0d0,  0d0,  0d0,  0d0 /)
         A(3,:) = (/  0d0, -1d0,  2d0,  4d0,  0d0,  0d0,  0d0,  0d0 /)
         A(4,:) = (/ -1d0,  0d0,  0d0,  2d0, -1d0,  0d0,  0d0,  0d0 /)
         A(5,:) = (/  0d0,  0d0,  0d0, -1d0,  3d0, -1d0,  0d0,  0d0 /)
         A(6,:) = (/  0d0,  0d0,  0d0,  0d0, -1d0,  3d0, -2d0,  0d0 /)
         A(7,:) = (/  0d0,  0d0,  0d0,  0d0,  0d0, -1d0,  2d0, -1d0 /)
         A(8,:) = (/  0d0,  0d0,  0d0,  0d0,  0d0,  0d0, -2d0,  2d0 /)
         
         ublk(1,1:2,1) = A(1,3:4)
         ublk(2,1:2,1) = A(2,3:4)
         ublk(1,1:2,2) = A(3,5:6)
         ublk(2,1:2,2) = A(4,5:6)
         ublk(1,1:2,3) = A(5,7:8)
         ublk(2,1:2,3) = A(6,7:8)
         ublk(1,1:2,4) = 0d0
         ublk(2,1:2,4) = 0d0

         dblk(1,1:2,1) = A(1,1:2)
         dblk(2,1:2,1) = A(2,1:2)
         dblk(1,1:2,2) = A(3,3:4)
         dblk(2,1:2,2) = A(4,3:4)
         dblk(1,1:2,3) = A(5,5:6)
         dblk(2,1:2,3) = A(6,5:6)
         dblk(1,1:2,4) = A(7,7:8)
         dblk(2,1:2,4) = A(8,7:8)

         lblk(1,1:2,1) = 0d0
         lblk(2,1:2,1) = 0d0
         lblk(1,1:2,2) = A(3,1:2)
         lblk(2,1:2,2) = A(4,1:2)
         lblk(1,1:2,3) = A(5,3:4)
         lblk(2,1:2,3) = A(6,3:4)
         lblk(1,1:2,4) = A(7,5:6)
         lblk(2,1:2,4) = A(8,5:6)
             
         itr_max = 50
         tol_abs = 1.0D-08
         tol_rel = 1.0D-08

         num_sweeps_factor = 1
         num_sweeps_solve = 3
         !exact = .true.
         exact = .false.
         verbose = .false.
         !verbose = .true.
         equilibrate = .true.
         write_mm = .true.
         !use_Bi_CG_Stab = .true.
         use_Bi_CG_Stab = .false.
         
         rhs1 = 1d0

         if (write_mm) then
            call write_MM_mxt( &
               nvar, nz, ublk, dblk, lblk, 'abtilu_2_mtx.txt')
            call write_MM_vec(nvar, nz, rhs1, 'abtilu_2_rhs.txt')
         end if
         
         ! get solution from DGESVX
         call get_lapack_solution(neq, A, rhs1, actual_soln1, verbose, ierr)
         if (ierr /= 0) stop 'test_abtilu_nvar2 failed in DGESVX'
         rhs1 = 1d0
         soln1 = 0d0
         x_error = norm2_of_diff(neq, actual_soln1, soln1)
         
         if (use_Bi_CG_Stab) then
            call solve_abtilu_with_Bi_CG_Stab( &
               nvar, nz, A, ublk, dblk, lblk, rhs1, &
               equilibrate, exact, &
               num_sweeps_factor, num_sweeps_solve, &
               itr_max, tol_rel, &
               soln1, verbose, ierr)
         else
            call solve_abtilu_with_mgmres( &
               nvar, nz, ublk, dblk, lblk, rhs1, &
               equilibrate, exact, &
               num_sweeps_factor, num_sweeps_solve, &
               itr_max, mr, tol_abs, tol_rel, &
               soln1, verbose, ierr)
         end if
         
         !if (verbose) then
         !   do i = 1, neq
         !      write ( *, '(2x,i8,2x,e26.16)' ) i, soln1(i)
         !   end do
         !end if

         if (write_mm) then
            call write_MM_vec(nvar, nz, soln1, 'abtilu_2_soln.txt')
            call write_MM_vec(nvar, nz, actual_soln1, 'abtilu_2_DGESVX_soln.txt')
         end if

         x_error = norm2_of_diff(neq, actual_soln1, soln1)
         if (x_error > 1d-12) then
            write ( *, '(a,i4,g19.11)' ) '  test_abtilu_nvar2 x=soln error ', round, x_error
            stop 'bad x_error test_abtilu_nvar1'
         else if (use_Bi_CG_Stab) then
            write ( *, '(a,i4,g19.11)' ) '  Bi_CG_Stab test_abtilu_nvar2 x=soln good ', round, x_error
         else
            write ( *, '(a,i4,g19.11)' ) '  mgmres test_abtilu_nvar2 x=soln good ', round, x_error
         end if

        contains
        
        real(dp) function norm2_of_diff(neq,a,b) result(v)
           integer, intent(in) :: neq
           real(dp), intent(in), dimension(:) :: a, b
           v = sqrt(sum(pow2(a(1:neq) - b(1:neq))))
        end function norm2_of_diff
        
      end subroutine test_abtilu_nvar2     


      subroutine Bi_CG_Stab( & ! based on http://www.netlib.org/templates/matlab/bicgstab.m
            matvec, psolve, b, x, max_iter, tol, error, iter, verbose, ierr)
         interface
            subroutine matvec(x,r) ! set r = A*x
               use const_def, only: dp
               real(dp), intent(in) :: x(:)
               real(dp), intent(out) :: r(:)
            end subroutine matvec
            subroutine psolve(x,z) ! set z = Precond*x
               use const_def, only: dp
               real(dp), intent(in) :: x(:)
               real(dp), intent(out) :: z(:)
            end subroutine psolve
         end interface
         real(dp), intent(in) :: b(:)
         real(dp), intent(inout) :: x(:) ! in: initial guess, out: soln
         integer, intent(in) :: max_iter
         real(dp), intent(in) :: tol
         real(dp), intent(out) :: error
         logical, intent(in) :: verbose
         integer, intent(out) :: iter, ierr         
         real(dp), dimension(:), allocatable :: &
            r, r_tld, v, p, s, t, p_hat, s_hat
         real(dp), parameter :: e = 1d-33
         real(dp) :: rho, rho_1, alpha, omega, beta, &
            rnrm2, bnrm2, snrm2 
         integer :: neq
         include 'formats'
         neq = size(b,dim=1)
         if(neq /= size(x, dim=1)) &
            stop "Error: Improper dimension of vector x in Bi_CG_Stab."
         allocate( &
            r(neq), r_tld(neq), v(neq), p(neq), s(neq), t(neq), &
            p_hat(neq), s_hat(neq))
         
         iter = 0
         bnrm2 = sqrt(dot_product(b,b))   
         if (bnrm2 == 0d0) bnrm2 = 1d0
         call matvec(x,r) ! r = A*x
         r = b - r ! r is residual
         rnrm2 = sqrt(dot_product(r,r))          
         error = rnrm2/bnrm2
         if (error < tol) return
         error = 1d0
         omega = 1d0  
         rho_1 = 1d0
         rho = 1d0
         r_tld = r                                         
         v = 0d0; p = 0d0                     
         do iter = 1, max_iter                   
            rho = dot_product(r_tld,r) 
            if (rho == 0d0) then
               ierr = 2
               return
            end if
            if (iter > 1) then                       
               beta = (rho/rho_1) * (alpha/omega)           
               p = r + beta * (p - omega*v) 
            else
               p = r
            end if
            call psolve(p,p_hat)   
            call matvec(p_hat,v)                        
            alpha = rho/dot_product(r_tld,v)                    
            s = r - alpha*v   
            snrm2 = sqrt(dot_product(s,s))      
            if (snrm2 < tol) then
               x = x + alpha*p_hat
               if (verbose) write(*,2) 'snrm2 < tol', iter, snrm2, tol
               ierr = 0
               return
            end if      
            call psolve(s,s_hat)                
            call matvec(s_hat,t)                
            omega = dot_product(t,s)/dot_product(t,t)        
            x = x + alpha*p + omega*s_hat                   
            r = s - omega*t                              
            rnrm2 = sqrt(dot_product(r,r))                  
            error = rnrm2/bnrm2
            if (error < tol) then
               if (verbose) write(*,2) 'error < tol', iter, error, tol
               ierr = 0
               return
            end if
            if (omega == 0d0) then
               if (verbose) write(*,2) 'omega == 0d0', iter
               ierr = -2
               return
            end if
            rho_1 = rho
            if (verbose) write(*,2) 'iter error rho', iter, error, rho
         end do                                                    
         ierr = 1
      end subroutine Bi_CG_Stab     


      subroutine test_Bi_CG_Stab
          integer, parameter :: m=4, n=4
          real(dp), dimension(1:m,1:n) :: A
          real(dp), dimension(1:m) :: x_calculated, x_actual, b
          real(dp) :: tol, error
          integer :: max_iter, iter, ierr
          logical :: verbose, use_dummy_psolve
          ierr = 0
          A(1,:) = (/  5d0, -3d0,  0d0,  0d0 /)
          A(2,:) = (/ -1d0,  4d0,  2d0,  0d0 /)
          A(3,:) = (/  0d0, -5d0,  2d0, -2d0 /)
          A(4,:) = (/  0d0,  0d0, -1d0,  6d0 /)

          b = 1d0
          use_dummy_psolve = .false.
          call psolve(b,x_actual)
          !use_dummy_psolve = .true.
          
          !x_actual(:) = [ 5d0, 8d0, 10d0, 12d0]
          !b = matmul(A,x_actual)
          
          x_calculated = 0d0 ! initial guess
          max_iter = 100
          tol = 1d-8
          verbose = .true.
          
          call Bi_CG_Stab( &
            matvec, psolve, b, x_calculated, max_iter, &
            tol, error, iter, verbose, ierr)
          if (ierr /= 0) then
             write(*,*) 'Bi_CG_Stab ierr', ierr
             return
          end if
          
          print*, "iterations required =", iter
          print*, "analytical solution =", x_actual
          print*, "  computed solution =", x_calculated
          print*, "     relative error =", abs(x_actual-x_calculated)/x_actual
          write(*,*)
          
          contains
          
         subroutine matvec(x, r) ! set r = A*x
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: r(:)
            r = matmul(A,x)
         end subroutine matvec
         
         subroutine psolve(x,z) ! set z = Precond*x
            use const_def, only: dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: z(:)
            integer :: ipiv(n), ierr
            real(dp) :: Acopy(n,n)
            z = x
            if (use_dummy_psolve) return
            ierr = 0
            Acopy = A
            call DGESV( n, 1, Acopy, m, ipiv, z, n, ierr )
            if (ierr /= 0) then
               write(*,*) 'DGESV failed', ierr
               stop 'test_Bi_CG_Stab'
            end if
         end subroutine psolve
          

      end subroutine test_Bi_CG_Stab



      
      end module abtilu
