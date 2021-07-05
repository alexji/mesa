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
      public :: test_abtilu, solve_with_mgmres, &
         factor_abtilu, solve_abtilu, block_tridiag_mv1, mgmres


      contains      
      
      
      subroutine test_abtilu()
         include 'formats'
         integer :: j
         do j=1,1
            !write(*,2) 'testing round', j
            !call test_abtilu_mgmres_nvar1(j)
            call test_abtilu_mgmres_nvar2(j)
         end do
         write(*,*)
         stop 'done test_abtilu'
      end subroutine test_abtilu
      
      subroutine solve_with_mgmres( &
            nvar, nz, ublk, dblk, lblk, rhs1, &
            itr_max, mr, exact, tol_abs, tol_rel, &
            num_sweeps_factor, num_sweeps_solve, &
            soln1, verbose, ierr)
         integer, intent(in) :: nvar, nz, itr_max, mr, &
            num_sweeps_factor, num_sweeps_solve
         real(dp), dimension(:,:,:), intent(in) :: & !(nvar,nvar,nz)
            ublk, dblk, lblk
         real(dp), dimension(:), intent(in) :: rhs1 ! (neq)
         logical, intent(in) :: exact, verbose
         real(dp), intent(in) :: tol_abs, tol_rel
         real(dp), dimension(:), intent(inout) :: soln1 ! (neq)
            ! input: initial guess (can be 0)
            ! output: final approximation
         integer, intent(out) :: ierr
         
         real(dp), dimension(:,:,:), allocatable :: & !(nvar,nvar,nz)
            Dhat, invDhat, invDhat_lblk, invDhat_ublk
         real(dp), dimension(:), allocatable :: & ! (neq)
            prev1, w1, x1, y1, r
         real(dp), allocatable :: v(:,:) ! (neq,mr+1)
         integer, allocatable :: ipiv(:,:) ! (nvar,nz)
         integer :: neq
         include 'formats'
         ierr = 0
         neq = nvar*nz
         allocate( &
            Dhat(nvar,nvar,nz), invDhat(nvar,nvar,nz), &
            invDhat_lblk(nvar,nvar,nz), invDhat_ublk(nvar,nvar,nz), &
            prev1(neq), w1(neq), x1(neq), y1(neq), r(neq), v(neq,mr+1), ipiv(nvar,nz))
         call create_preconditioner_abtilu(ierr)     
         if (ierr /= 0) stop 'failed in create_preconditioner_abtilu_mgmres'
         call mgmres ( &     
            neq, matvec_abtilu_mgmres, solve_abtilu_mgmres, &
            soln1, rhs1, r, v, itr_max, mr, tol_abs, tol_rel, verbose )
         
         contains
         
         subroutine create_preconditioner_abtilu(ierr)
            integer, intent(out) :: ierr
            include 'formats'
            call factor_abtilu( &
                nvar, nz, lblk, dblk, ublk, & ! input
                num_sweeps_factor, exact, & ! input
                Dhat, ipiv, & ! work
                invDhat, invDhat_lblk, invDhat_ublk, ierr) ! output
         end subroutine create_preconditioner_abtilu
     
         subroutine solve_abtilu_mgmres(v1)
            real(dp), intent(inout) :: v1(:) ! (neq)
            include 'formats'
            call solve_abtilu( &
                nvar, nz, invDhat, invDhat_lblk, invDhat_ublk, v1, & ! input
                num_sweeps_solve, exact, & ! input
                prev1, w1, x1, y1, & ! work
                v1, ierr) ! output
         end subroutine solve_abtilu_mgmres

         subroutine matvec_abtilu_mgmres(b1, r1) ! set r = Jacobian*b
            real(dp), intent(in) :: b1(:) ! (neq)
            real(dp), intent(out) :: r1(:) ! (neq)
            integer :: k
            include 'formats'
            call block_tridiag_mv1(nvar, nz, lblk, dblk, ublk, b1, r1)
         end subroutine matvec_abtilu_mgmres
         
      end subroutine solve_with_mgmres
         
         
!*****************************************************************************
!*****************************************************************************
!
!  support routines for abtilu
!
!*****************************************************************************
!*****************************************************************************


      subroutine get_diagonal_scaling_vectors( &
            nvar, nz, lblk, dblk, ublk, & ! input
            max_iters, eps, & ! input
            DR, DC, & ! work
            D1, D2, ierr) ! output
         ! see Amestoy et al, 2008, A Parallel Matrix Scaling Algorithm
         !   to solve A*x = b, 
         !      instead solve (D1*A*D2)(D2^-1*x) = (D1*b)
         !      Ahat = D1*A*D2
         !      bhat = D1*b, 
         !      solve Ahat*xhat = bhat
         !      then x = D2*xhat
         integer, intent(in) :: nvar, nz, max_iters
         real(dp), intent(in) :: eps
         real(dp), dimension(:,:,:), intent(in) :: & ! input (nvar,nvar,nz)
            lblk, dblk, ublk
         real(dp), dimension(:), intent(out) :: DR, DC ! work (nvar*nz)
         real(dp), dimension(:), intent(out) :: D1, D2 ! output (nvar*nz)
         integer, intent(out) :: ierr
         integer :: iter, k, i, j, neq
         logical :: row_converged, col_converged
         include 'formats'
         ierr = 0
         neq = nvar*nz
         D1(:) = 1d0
         D2(:) = 1d0
         do iter = 1, max_iters
            !$OMP PARALLEL DO PRIVATE(k,i,j)
            do k=1,nz
               do i=1,nvar
                  j = diag_index(i,k)
                  DR(j) = row_max_abs(i,k)
                  DC(j) = col_max_abs(i,k)
                  D1(j) = D1(j)/sqrt(DR(j))
                  D2(j) = D2(j)/sqrt(DC(j))
               end do
            end do
            !$OMP END PARALLEL DO       
            row_converged = abs(1d0 - maxval(D1(1:neq))) <= eps
            if (row_converged) then
               col_converged = abs(1d0 - maxval(D2(1:neq))) <= eps
               if (col_converged) exit
            end if
         end do
         
         contains
         
         integer function diag_index(i,k) result(j)
            integer, intent(in) :: i, k
            j = (k-1)*nvar + i
         end function diag_index
         
         real(dp) function row_max_abs(i,k)
            integer, intent(in) :: i, k
            integer :: s
            real(dp) :: lmax = 0d0, dmax = 0d0, umax = 0d0
            if (k > 1) then
               s = (k-2)*nvar
               do j=1,nvar
                  lmax = max(lmax, abs(lblk(i,j,k))*D2(s+j))
               end do
            end if
            if (k < nz) then
               s = k*nvar
               do j=1,nvar
                  umax = max(umax, abs(ublk(i,j,k))*D2(s+j))
               end do
            end if
            s = (k-1)*nvar
            do j=1,nvar
               dmax = max(dmax, abs(dblk(i,j,k))*D2(s+j))
            end do
            row_max_abs = D1(s+i)*max(lmax, dmax, umax)
         end function row_max_abs
         
         real(dp) function col_max_abs(i,k)
            integer, intent(in) :: i, k
            integer :: s
            real(dp) :: lmax = 0d0, dmax = 0d0, umax = 0d0
            if (k > 1) then
               s = (k-2)*nvar
               do j=1,nvar
                  lmax = max(lmax, D1(s+j)*abs(lblk(i,j,k)))
               end do
            end if
            if (k < nz) then
               s = k*nvar
               do j=1,nvar
                  umax = max(umax, D1(s+j)*abs(ublk(i,j,k)))
               end do
            end if
            s = (k-1)*nvar
            do j=1,nvar
               dmax = max(dmax, D1(s+j)*abs(dblk(i,j,k)))
            end do
            col_max_abs = max(lmax, dmax, umax)*D2(s+i)
         end function col_max_abs
            
      end subroutine get_diagonal_scaling_vectors
      
      
      subroutine factor_abtilu( &
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
         
         
         ! optionally provide initial guess for Dhat.  chow-patel, section 2.3
         ! e.g., Dhat from previous star solver newton iteration.
         
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
            call mm_minus(lblk(:,:,k), invDhat_ublk(:,:,k-1), Dhat(:,:,k)) ! c := c - a*b            
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
                  
      end subroutine factor_abtilu
      
      
      subroutine solve_abtilu( &
            nvar, nz, & ! input
            invDhat, invDhat_lblk, invDhat_ublk, b1, & ! input
            num_sweeps, exact, & ! input
            prev1, w1, x1, y1, & ! work
            z1, ierr) ! output
         ! kashi phd thesis, pg 122, redone for block tridiagonal
         integer, intent(in) :: nvar, nz, num_sweeps
         logical, intent(in) :: exact
         real(dp), dimension(:,:,:), intent(in) :: & ! input (nvar,nvar,nz)
            invDhat, invDhat_lblk, invDhat_ublk
         real(dp), intent(in) :: b1(:) ! input (nvar*nz)
         real(dp), dimension(:), intent(out) :: & ! (nvar*nz)
            prev1, w1, x1, y1, & ! work
            z1 ! output
         integer, intent(out) :: ierr
         logical :: incomplete
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
            integer, intent(in) :: k
            integer :: s00, sp1, j
            s00 = (k-1)*nvar
            sp1 = s00 + nvar
            !$omp simd
            do j=1,nvar ! z = y
               z1(s00+j) = y1(s00+j)
            end do
            if (k == nz) return
            call mv_minus( & ! z(k) = z(k) - a(k)*prev1(k+1)
               invDhat_ublk(:,:,k),prev1(sp1+1:sp1+nvar),z1(s00+1:s00+nvar))
         end subroutine set_z
         
      end subroutine solve_abtilu


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
            call matvec ( v(1:n,k), v(1:n,k+1) )
            call psolve ( v(1:n,k+1) ) ! apply to result of matvec
            av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
            if (is_bad(av)) then
               write(*,*) 'bad result from matvec/psolve: is_bad(av)'
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
      

      subroutine test_abtilu_mgmres_nvar1(round)
      !
      !    A = 
      !      2  -1   0  0 
      !     -1  2   -1  0 
      !
      !      0 -1    2 -1
      !      0  0   -1  2 
      !         
         use utils_lib, only: fill_with_NaNs, fill_with_NaNs_2D, fill_with_NaNs_3D
         integer, intent(in) :: round
         integer, parameter :: nvar = 1, nz = 4, neq = nvar*nz
         real(dp), dimension(:,:,:), allocatable :: &
            ublk, dblk, lblk
         real(dp), dimension(:), allocatable :: &
            soln1, rhs1, actual_soln1
         integer :: test, i, k, mr, itr_max, shft, ierr, &
            num_sweeps_factor, num_sweeps_solve
         logical :: exact, verbose
         real(dp) :: tol_abs, tol_rel, x_error, soln_error
         include 'formats'
         
         mr = neq - 1
         allocate( &
            ublk(nvar,nvar,nz), dblk(nvar,nvar,nz), lblk(nvar,nvar,nz), &
            soln1(neq), rhs1(neq), actual_soln1(neq))
         
         call fill_with_NaNs(soln1)
         call fill_with_NaNs(rhs1)
         call fill_with_NaNs(actual_soln1)
         call fill_with_NaNs_3D(ublk)
         call fill_with_NaNs_3D(dblk)
         call fill_with_NaNs_3D(lblk)
         
         ublk = -1d0     
         dblk = 2d0
         lblk = -1d0
         
         rhs1 = 1d0
         soln1 = 0d0
         actual_soln1 = (/ 2d0, 3d0, 3d0, 2d0 /)
             
         itr_max = 1 ! 20
         tol_abs = 1.0D-08
         tol_rel = 1.0D-08

         num_sweeps_factor = 1
         num_sweeps_solve = 3
         !exact = .true.
         exact = .false.
         verbose = .false.

         !write ( *, '(a)' ) ' '
         !write ( *, '(a)' ) 'test_abtilu_mgmres_nvar1'
         x_error = norm2_of_diff(neq, actual_soln1, soln1)
         !write ( *, '(a,g14.6)' ) '  Before solving, X_ERROR = ', x_error
         
         call solve_with_mgmres( &
            nvar, nz, ublk, dblk, lblk, rhs1, &
            itr_max, mr, exact, tol_abs, tol_rel, &
            num_sweeps_factor, num_sweeps_solve, &
            soln1, verbose, ierr)

         x_error = norm2_of_diff(neq, actual_soln1, soln1)
         write ( *, '(a,i4,g19.11)' ) '  test_abtilu_mgmres_nvar1 x=soln error ', round, x_error
         if (x_error > 1d-12) stop 'bad x_error test_abtilu_mgmres_nvar1'
         !write ( *, '(a)' ) '  x:'
         !do i = 1, neq
         !   write ( *, '(2x,i8,2x,g14.6)' ) i, soln1(i)
         !end do

        contains
        
        real(dp) function norm2_of_diff(neq,a,b) result(v)
           integer, intent(in) :: neq
           real(dp), intent(in), dimension(:) :: a, b
           v = sqrt(sum(pow2(a(1:neq) - b(1:neq))))
        end function norm2_of_diff
        
      end subroutine test_abtilu_mgmres_nvar1      
      

      subroutine test_abtilu_mgmres_nvar2(round)
      !
      !    A = 
      !      2  0    0 -1    0  0    0  0   
      !      0  2   -1  0    0  0    0  0  
      !
      !      0 -1    2  0    0  0    0  0 
      !     -1  0    0  2   -1  0    0  0 
      !
      !      0  0    0 -1    2 -1    0  0 
      !      0  0    0  0   -1  2   -1  0 
      !
      !      0  0    0  0    0 -1    2 -1 
      !      0  0    0  0    0  0   -1  2 
      !         
         use utils_lib, only: fill_with_NaNs, fill_with_NaNs_2D, fill_with_NaNs_3D
         integer, intent(in) :: round
         integer, parameter :: nvar = 2, nz = 4, neq = nvar*nz
         real(dp), dimension(:,:,:), allocatable :: &
            ublk, dblk, lblk
         real(dp), dimension(:), allocatable :: &
            soln1, rhs1, actual_soln1
         integer :: test, i, k, mr, itr_max, shft, ierr, &
            num_sweeps_factor, num_sweeps_solve
         logical :: exact, verbose
         real(dp) :: tol_abs, tol_rel, x_error, soln_error
         include 'formats'
         
         mr = neq - 1
         allocate( &
            ublk(nvar,nvar,nz), dblk(nvar,nvar,nz), lblk(nvar,nvar,nz), &
            soln1(neq), rhs1(neq), actual_soln1(neq))
         
         call fill_with_NaNs(soln1)
         call fill_with_NaNs(rhs1)
         call fill_with_NaNs(actual_soln1)
         call fill_with_NaNs_3D(ublk)
         call fill_with_NaNs_3D(dblk)
         call fill_with_NaNs_3D(lblk)
         
         ublk(1:2,1,1) = (/ 0d0, -1d0 /)
         ublk(1:2,2,1) = (/ -1d0, 0d0 /)
         ublk(1:2,1,2) = (/ 0d0, -1d0 /)
         ublk(1:2,2,2) = (/ 0d0, 0d0 /)
         ublk(1:2,1,3) = (/ 0d0, -1d0 /)
         ublk(1:2,2,3) = (/ 0d0, 0d0 /)
         ublk(:,:,4) = 0

         dblk(1:2,1,1) = (/ 2d0, 0d0 /)
         dblk(1:2,2,1) = (/ 0d0, 2d0 /)
         dblk(1:2,1,2) = (/ 2d0, 0d0 /)
         dblk(1:2,2,2) = (/ 0d0, 2d0 /)
         dblk(1:2,1,3) = (/ 2d0, -1d0 /)
         dblk(1:2,2,3) = (/ -1d0, 2d0 /)
         dblk(1:2,1,4) = (/ 2d0, -1d0 /)
         dblk(1:2,2,4) = (/ -1d0, 2d0 /)

         lblk(:,:,1) = 0d0
         lblk(1:2,1,2) = (/ 0d0, -1d0 /)
         lblk(1:2,2,2) = (/ -1d0, 0d0 /)
         lblk(1:2,1,3) = (/ 0d0, 0d0 /)
         lblk(1:2,2,3) = (/ -1d0, 0d0 /)
         lblk(1:2,1,4) = (/ 0d0, 0d0 /)
         lblk(1:2,2,4) = (/ -1d0, 0d0 /)
      
         rhs1 = 1d0
         soln1 = 0d0
         actual_soln1 = (/ 3d0, 1d0, 1d0, 5d0, 6d0, 6d0, 5d0, 3d0 /)
             
         itr_max = 1 ! 20
         mr = 3 ! nvar*nz - 1
         tol_abs = 1.0D-08
         tol_rel = 1.0D-08

         num_sweeps_factor = 1
         num_sweeps_solve = 3
         !exact = .true.
         exact = .false.
         verbose = .false.

         x_error = norm2_of_diff(neq, actual_soln1, soln1)

         ierr = 0
         
         call solve_with_mgmres( &
            nvar, nz, ublk, dblk, lblk, rhs1, &
            itr_max, mr, exact, tol_abs, tol_rel, &
            num_sweeps_factor, num_sweeps_solve, &
            soln1, verbose, ierr)

         x_error = norm2_of_diff(neq, actual_soln1, soln1)
         write ( *, '(a,i4,g19.11)' ) '  test_abtilu_mgmres_nvar2 x=soln error ', round, x_error
         if (x_error > 1d-12) stop 'bad x_error test_abtilu_mgmres_nvar2'

         if (round == 1) &
            call write_MM( &
               nvar, nz, ublk, dblk, lblk, rhs1, &
               soln1, 'test_abtilu_mgmres_nvar2', ierr)

        contains
        
        real(dp) function norm2_of_diff(neq,a,b) result(v)
           integer, intent(in) :: neq
           real(dp), intent(in), dimension(:) :: a, b
           v = sqrt(sum(pow2(a(1:neq) - b(1:neq))))
        end function norm2_of_diff
        
      end subroutine test_abtilu_mgmres_nvar2     
      
      subroutine write_MM( &
            nvar, nz, ublk, dblk, lblk, rhs1, &
            soln1, filename, ierr)
         integer, intent(in) :: nvar, nz
         real(dp), dimension(:,:,:), intent(in) :: & !(nvar,nvar,nz)
            ublk, dblk, lblk
         real(dp), dimension(:), intent(in) :: rhs1 ! (neq)
         real(dp), dimension(:), intent(in) :: soln1 ! (neq)
         character (len=*), intent(in) :: filename
         integer, intent(out) :: ierr
         
         write(*,*) 'write_MM'
         ierr = 0
         
      end subroutine write_MM

      
      end module abtilu
