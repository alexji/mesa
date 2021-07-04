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
      public :: test_abtilu


      contains      
      
      
      subroutine test_abtilu()
         call test_ILU_CR_mgmres_nvar1()
         call test_abtilu_mgmres_nvar1()
         write(*,*)
         stop 'done test_abtilu'
      end subroutine test_abtilu


!*****************************************************************************
!*****************************************************************************
!
!  support routines for abtilu
!
!*****************************************************************************
!*****************************************************************************
      
      
      subroutine factor_abtilu( &
            nvar, nz, lblk, dblk, ublk, & ! input
            num_sweeps, exact, & ! input
            Dhat, ipiv, & ! work
            invDhat, invDhat_lblk, invDhat_ublk, ierr) ! output
         ! kashi pg 121, Algorithm 12, redone for block tridiagonal
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
            ! Dhat(k) = dblk(k) - lblk(k-1)*invDhat(k-1)*ublk(k-1)
            !$OMP PARALLEL DO PRIVATE(k) IF(incomplete)
            do k = 2, nz
               call set_Dhat(k)
               if (exact) then
                  call set_invDhat(k,op_err)
                  if (op_err /= 0) ierr = op_err
               end if
            end do
            !$OMP END PARALLEL DO 
            !write(*,2) 'Dhat', swp, Dhat
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
         
         ! create matrix products that are useful for apply
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
         
         subroutine set_Dhat(k)
            integer, intent(in) :: k
            call set_invDhat_ublk(k-1)
            ! Dhat(k) = dblk(k) - lblk(k)*invDhat_ublk(k-1)
            call copy_dblk_to_Dhat(k)
            call mm_minus(lblk(:,:,k), invDhat_ublk(:,:,k-1), Dhat(:,:,k)) ! c := c - a*b
         end 
         
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
            nvar, nz, invDhat, lblk, ublk, b1, & ! input
            num_sweeps, exact, & ! input
            x1, y1, & ! work
            z1, ierr) ! output
         ! kashi pg 122, redone for block tridiagonal
         integer, intent(in) :: nvar, nz, num_sweeps
         logical, intent(in) :: exact
         real(dp), dimension(:,:,:), intent(in) :: & ! input (nvar,nvar,nz)
            invDhat, lblk, ublk
         real(dp), intent(in) :: b1(:) ! input (nvar*nz)
         real(dp), dimension(:), intent(out) :: & ! (nvar*nz)
            x1, y1, z1 ! output (nvar*nz)
         integer, intent(out) :: ierr
         logical :: incomplete
         integer :: swp, j, k, neq
         include 'formats'
         ierr = 0
         incomplete = .not. exact
         neq = nvar*nz
         
         y1(:) = 0d0 ! initialize y to 0
         
         do swp=1,num_sweeps         
            ! analogous to forward elimination phase for tridiagonal solve
            !$OMP PARALLEL DO PRIVATE(k) IF(incomplete)
            do k=1,nz
               call set_y(k)
            end do
            !$OMP END PARALLEL DO  
            !write(*,2) 'y1', swp, y1  
            if (exact) exit             
         end do
         
         !$omp simd
         do j=1,neq
            z1(j) = y1(j) ! initialize z to y
         end do
         
         do swp=1,num_sweeps               
            ! analogous to backward substitution phase for tridiagonal solve
            !$OMP PARALLEL DO PRIVATE(k) IF(incomplete)
            do k=nz,1,-1
               call set_z(k)
            end do
            !$OMP END PARALLEL DO   
            !write(*,2) 'z1', swp, z1                
            if (exact) exit             
         end do
         
         contains
                  
         subroutine set_y(k) ! y(k) = invDhat(k)*(b(k) - lblk(k)*y(k-1)) eq 5.9
            integer, intent(in) :: k
            integer :: s00, sm1, j
            s00 = (k-1)*nvar 
            sm1 = s00 - nvar
            if (k == 1) then ! lblk(1) = 0
               !$omp simd
               do j=1,nvar
                  x1(j) = b1(j)
               end do
            else
               call mv_0(lblk(:,:,k),y1(sm1+1:sm1+nvar),x1(s00+1:s00+nvar))
               !$omp simd
               do j=1,nvar
                  x1(s00+j) = b1(s00+j) - x1(s00+j)
               end do
            end if
            call mv_0(invDhat(:,:,k),x1(s00+1:s00+nvar),y1(s00+1:s00+nvar))
         end subroutine set_y
         
         subroutine set_z(k) ! z(k) = y(k) - invDhat(k)*ublk(k)*z(k+1) eq 5.10
            integer, intent(in) :: k
            integer :: j, s00, sp1
            s00 = (k-1)*nvar
            sp1 = s00 + nvar
            if (k == nz) then
               !$omp simd
               do j=1,nvar
                  z1(s00+j) = y1(s00+j)
               end do
            else
               ! x(k) = ublk(k)*z(k+1); x(k) = invDhat*x(k); z(k) = y(k) - x(k)
               call mv_0(ublk(:,:,k),z1(sp1+1:sp1+nvar),x1(s00+1:s00+nvar))
               call mv_self(invDhat(:,:,k),x1(s00+1:s00+nvar))
               !$omp simd
               do j=1,nvar
                  z1(s00+j) = y1(s00+j) - x1(s00+j)
               end do
            end if
         end subroutine set_z
         
      end subroutine solve_abtilu
      
      
!*****************************************************************************
!*****************************************************************************
!
!  testing abtilu
!
!*****************************************************************************
!*****************************************************************************
      

      subroutine test_abtilu_mgmres_nvar1()
      !
      !    A = 
      !      2  -1   0  0 
      !     -1  2   -1  0 
      !
      !      0 -1    2 -1
      !      0  0   -1  2 
      !         
         integer, parameter :: nvar = 1, nz = 4, neq = nvar*nz
         real(dp), dimension(nvar,nvar,nz) :: &
            ublk, dblk, lblk, Dhat, invDhat, &
            Lhat, invDhat_lblk, invDhat_ublk
         real(dp), dimension(nvar,nz) :: vec1, vec2, vec3
         integer :: ipiv(nvar,nz)
         real(dp), dimension(neq) :: w1, x1, y1, xmgmres1, rhs1, soln1
         integer :: seed = 123456789
         integer :: test, i, k, mr, itr_max, shft, ierr, &
            num_sweeps_factor, num_sweeps_solve
         logical :: exact
         real(dp) :: tol_abs, tol_rel, x_error, soln_error
         include 'formats'
         
         ublk = -1d0     
         dblk = 2d0
         lblk = -1d0
         
         rhs1 = 1d0
         soln1 = (/ 2d0, 3d0, 3d0, 2d0 /)
             
          itr_max = 1 ! 20
          mr = neq - 1
          tol_abs = 1.0D-08
          tol_rel = 1.0D-08
          
          exact = .true.
          num_sweeps_factor = nz ! for exact, need num_sweeps_factor = nz
          num_sweeps_solve = nz ! 3
             
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'test_abtilu_mgmres_nvar1'
          xmgmres1 = 0d0
          x_error = norm2_of_diff(neq, soln1, xmgmres1)
          write ( *, '(a,g14.6)' ) '  Before solving, X_ERROR = ', x_error
         
         ierr = 0
         call create_preconditioner_abtilu(ierr)     
         if (ierr /= 0) stop 'failed in create_preconditioner_abtilu_mgmres'
         call mgmres ( &     
            neq, matvec_abtilu_mgmres, solve_abtilu_mgmres, &
            xmgmres1, rhs1, itr_max, mr, tol_abs, tol_rel )

          x_error = norm2_of_diff(neq, soln1, xmgmres1)
          write ( *, '(a,g14.6)' ) '  x=soln error ', x_error
          write ( *, '(a)' ) '  x:'
          do i = 1, neq
             write ( *, '(2x,i8,2x,g14.6)' ) i, xmgmres1(i)
          end do

        contains
        
        real(dp) function norm2_of_diff(neq,a,b) result(v)
           integer, intent(in) :: neq
           real(dp), intent(in), dimension(:) :: a, b
           v = sqrt(sum(pow2(a(1:neq) - b(1:neq))))
        end function norm2_of_diff
        
        subroutine create_preconditioner_abtilu(ierr)
           integer, intent(out) :: ierr
           include 'formats'
           call factor_abtilu( &
               nvar, nz, lblk, dblk, ublk, & ! input
               num_sweeps_factor, exact, & ! input
               Dhat, ipiv, & ! work
               invDhat, invDhat_lblk, invDhat_ublk, ierr) ! output
            if (nvar == 1) then
               write(*,1) 'Dhat', Dhat
            end if
        end subroutine create_preconditioner_abtilu
     
        subroutine solve_abtilu_mgmres(v1)
           real(dp), intent(inout) :: v1(:) ! (neq)
           include 'formats'
           write(*,1) 'solve input', v1(1:neq)
           call solve_abtilu( &
               nvar, nz, invDhat, lblk, ublk, v1, & ! input
               num_sweeps_solve, exact, & ! input
               x1, y1, & ! work
               v1, ierr) ! output
           write(*,1) 'solve output', v1(1:neq)
        end subroutine solve_abtilu_mgmres

        subroutine matvec_abtilu_mgmres(b1, r1) ! set r = Jacobian*b
           real(dp), intent(in) :: b1(:) ! (neq)
           real(dp), intent(out) :: r1(:) ! (neq)
           integer :: k
           include 'formats'
           write(*,1) 'matvec input', b1(1:neq)
           call BTD_mv1(nvar, nz, lblk, dblk, ublk, b1, r1)
           write(*,1) 'matvec output', r1(1:neq)
        end subroutine matvec_abtilu_mgmres
        
      end subroutine test_abtilu_mgmres_nvar1      
      

      subroutine test_abtilu_mgmres_nvar2()
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
         integer, parameter :: nvar = 2, nz = 4, neq = nvar*nz
         real(dp), dimension(nvar,nvar,nz) :: &
            ublk, dblk, lblk, Dhat, &
            Lhat, invDhat, invDhat_lblk, invDhat_ublk
         real(dp), dimension(nvar,nz) :: vec1, vec2, vec3
         integer :: ipiv(nvar,nz)
         real(dp), dimension(neq) :: w1, x1, y1, rhs1, soln1
         integer :: seed = 123456789
         integer :: test, i, k, mr, itr_max, shft, ierr, &
            num_sweeps_factor, num_sweeps_solve
         logical :: exact
         real(dp) :: tol_abs, tol_rel, x_error, soln_error
         include 'formats'
         
         if (.false.) then ! identity matrix
            
            ublk = 0d0
            lblk = 0d0
            
            dblk(1:2,1,1) = (/ 1d0, 0d0 /)
            dblk(1:2,2,1) = (/ 0d0, 1d0 /)
            dblk(1:2,1,2) = (/ 1d0, 0d0 /)
            dblk(1:2,2,2) = (/ 0d0, 1d0 /)
            dblk(1:2,1,3) = (/ 1d0, 0d0 /)
            dblk(1:2,2,3) = (/ 0d0, 1d0 /)
            dblk(1:2,1,4) = (/ 1d0, 0d0 /)
            dblk(1:2,2,4) = (/ 0d0, 1d0 /)

            rhs1 = 1d0
            soln1 = 1d0
         
         else
         
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
            soln1 = (/ &
                3d0, 1d0, &
                1d0, 5d0, &
                6d0, 6d0, &
                5d0, 3d0 /)
             
         end if
             
          itr_max = 1 ! 20
          mr = 3 ! nvar*nz - 1
          tol_abs = 1.0D-08
          tol_rel = 1.0D-08
          
          num_sweeps_factor = 1
          num_sweeps_solve = 3
          exact = .true.
             
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'test_abtilu_mgmres_nvar2'
          x1 = 0d0
          x_error = norm2_of_diff(neq, soln1, x1)
          write ( *, '(a,g14.6)' ) '  Before solving, X_ERROR = ', x_error
         
         ierr = 0
         call create_preconditioner_abtilu(ierr)     
         if (ierr /= 0) stop 'failed in create_preconditioner_abtilu_mgmres'
         write(*,*) ' call mgmres'
         call mgmres ( &     
            neq, matvec_abtilu_mgmres, solve_abtilu_mgmres, &
            x1, rhs1, itr_max, mr, tol_abs, tol_rel )

          x_error = norm2_of_diff(neq, soln1, x1)
          write ( *, '(a,g14.6)' ) '  x=soln error ', x_error
          write ( *, '(a)' ) '  x:'
          do i = 1, neq
             write ( *, '(2x,i8,2x,g14.6)' ) i, x1(i)
          end do

        contains
        
        real(dp) function norm2_of_diff(neq,a,b) result(v)
           integer, intent(in) :: neq
           real(dp), intent(in), dimension(:) :: a, b
           v = sqrt(sum(pow2(a(1:neq) - b(1:neq))))
        end function norm2_of_diff
        
        subroutine create_preconditioner_abtilu(ierr)
           integer, intent(out) :: ierr
           call factor_abtilu( &
               nvar, nz, lblk, dblk, ublk, & ! input
               num_sweeps_factor, exact, & ! input
               Dhat, ipiv, & ! work
               invDhat, invDhat_lblk, invDhat_ublk, ierr) ! output
        end subroutine create_preconditioner_abtilu
     
        subroutine solve_abtilu_mgmres(v1)
           real(dp), intent(inout) :: v1(:) ! (neq)
           include 'formats'
           write(*,1) 'solve input', v1(1:neq)
           call solve_abtilu( &
               nvar, nz, invDhat, lblk, ublk, v1, & ! input
               num_sweeps_solve, exact, & ! input
               x1, y1, & ! work
               v1, ierr) ! output
           write(*,1) 'solve output', x1(1:neq)
        end subroutine solve_abtilu_mgmres

        subroutine matvec_abtilu_mgmres(b1, r1) ! set r = Jacobian*b
           real(dp), intent(in) :: b1(:) ! (neq)
           real(dp), intent(out) :: r1(:) ! (neq)
           include 'formats'
           write(*,1) 'matvec input', b1(1:neq)
           call BTD_mv1(nvar, nz, lblk, dblk, ublk, b1, r1)
           write(*,1) 'matvec output', r1(1:neq)
        end subroutine matvec_abtilu_mgmres
        
      end subroutine test_abtilu_mgmres_nvar2     


!*****************************************************************************
!*****************************************************************************
!
!  testing sparse arrays
!
!*****************************************************************************
!*****************************************************************************
      

      subroutine test_ILU_CR_mgmres_nvar1()
      !
      !    A = 
      !      2 -1    0  0    
      !     -1  2   -1  0    
      !
      !      0 -1    2 -1    
      !      0  0   -1  2   
      !      

        integer, parameter :: n = 4
        integer, parameter :: nz_num = 10
        ! nonzero values are sorted by row
        real(dp), dimension(nz_num) :: a = (/ &
           2.0D+00,  -1.0D+00, & ! row 1
           -1.0D+00,  2.0D+00, -1.0D+00, & ! row 2 
           -1.0D+00,  2.0D+00, -1.0D+00, & ! row 3
           -1.0D+00,  2.0D+00 /) ! row 4
        ! JA stores the column index of the nonzero value
        integer, dimension(nz_num) :: ja = (/ &
          1, 2, & ! row 1
          1, 2, 3, & ! row 2
          2, 3, 4, & ! row 3
          3, 4 /) ! row 4
        ! the entries in A and JA that correspond to row I occur in indices
        ! IA[I] through IA[I+1]-1.
        integer, dimension(n+1) :: ia = (/ &
          1, & ! row 1
          3, & ! row 2
          6, & ! row 3
          9, & ! row 4
          11  /)
        integer :: i
        integer :: itr_max          
        integer :: mr
        real(dp), dimension(n) :: rhs
        integer :: seed = 123456789
        integer :: test
        real(dp) :: tol_abs
        real(dp) :: tol_rel
        real(dp) :: x_error, soln_error
        real(dp) :: x_estimate(n)
        real(dp), dimension(n) :: x_exact = (/ &
          2d0, 3d0, 3d0, 2d0 /)
        integer :: ua(n)
        real(dp) :: l(nz_num+2) ! (ia(n+1)+1)
        real(dp) :: Ax_for_soln(n)
        
        include 'formats'
        
        rhs = 1d0
        
        if (nz_num+2 /= ia(n+1)+1) then
           stop 'bad counting test_ILU_CR_mgmres'
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'test_ILU_CR_mgmres_nvar1'

        do test = 1, 1 ! 2

          if ( test == 1 ) then

            x_estimate(1:n) = 0.0D+00

          else

            write ( *, '(a)' ) '  Second try, use random initial vector:'

            call r8vec_uniform_01 ( n, seed, x_estimate )

          end if

          x_error = sqrt ( sum ( ( x_exact(1:n) - x_estimate(1:n) )**2 ) )

          write ( *, '(a,g14.6)' ) '  Before solving, X_ERROR = ', x_error

          itr_max = 20
          mr = n - 1
          tol_abs = 1.0D-08
          tol_rel = 1.0D-08
                    
         call create_preconditioner_ILU_CR_mgmres()            
         !write(*,*) ' call mgmres'
         call mgmres ( &     
            n, matvec_ILU_CR_mgmres, apply_ILU_CR_mgmres, &
            x_estimate, rhs, itr_max, mr, tol_abs, tol_rel )
         call ax_cr ( n, nz_num, ia, ja, a, x_estimate, Ax_for_soln )

          soln_error = sqrt ( sum ( ( rhs(1:n) - Ax_for_soln(1:n) )**2 ) )
          x_error = sqrt ( sum ( ( x_exact(1:n) - x_estimate(1:n) )**2 ) )

          write ( *, '(a,g14.6)' ) '  x=soln error ', x_error
          write ( *, '(a)' ) '  x:'
          do i = 1, n
            write ( *, '(2x,i8,2x,g14.6)' ) i, x_estimate(i)
          end do

        end do

        contains
        
        subroutine create_preconditioner_ILU_CR_mgmres()
           call rearrange_cr ( n, nz_num, ia, ja, a )
           call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )
           call ilu_cr ( n, nz_num, ia, ja, a, ua, l ) 
        end subroutine create_preconditioner_ILU_CR_mgmres
     
        subroutine apply_ILU_CR_mgmres(x) ! apply preconditioner to x
           real(dp), intent(inout) :: x(:) ! (neq)
           include 'formats'
           write(*,1) 'solve input', x
           call lus_cr ( n, nz_num, ia, ja, l, ua, x, x ) 
           write(*,1) 'solve output', x
        end subroutine apply_ILU_CR_mgmres

        subroutine matvec_ILU_CR_mgmres(x, r) ! set r = Jacobian*x
           real(dp), intent(in) :: x(:) ! (neq)
           real(dp), intent(out) :: r(:) ! (neq)
           include 'formats'
           write(*,1) 'matvec input', x(1:n)
           call ax_cr ( n, nz_num, ia, ja, a, x, r ) 
           write(*,1) 'matvec output', r(1:n)
        end subroutine matvec_ILU_CR_mgmres
        
      end subroutine test_ILU_CR_mgmres_nvar1
      

      subroutine test_ILU_CR_mgmres_nvar2()
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

        integer, parameter :: n = 8
        integer, parameter :: nz_num = 20
        ! nonzero values are sorted by row
        real(dp), dimension(nz_num) :: a = (/ &
           2.0D+00, -1.0D+00, & ! row 1
           2.0D+00, -1.0D+00, & ! row 2 
          -1.0D+00,  2.0D+00, & ! row 3
          -1.0D+00,  2.0D+00, -1.0D+00, & ! row 4
          -1.0D+00,  2.0D+00, -1.0D+00, & ! row 5
          -1.0D+00,  2.0D+00, -1.0D+00, & ! row 6
          -1.0D+00,  2.0D+00, -1.0D+00, & ! row 7
          -1.0D+00,  2.0D+00 /) ! row 8
        ! JA stores the column index of the nonzero value
        integer, dimension(nz_num) :: ja = (/ &
          1, 4, & ! row 1
          2, 3, & ! row 2
          2, 3, & ! row 3
          1, 4, 5, & ! row 4
          4, 5, 6, & ! row 5
          5, 6, 7, & ! row 6
          6, 7, 8, & ! row 7
          7, 8 /) ! row 8
        ! the entries in A and JA that correspond to row I occur in indices
        ! IA[I] through IA[I+1]-1.
        integer, dimension(n+1) :: ia = (/ &
          1, & ! row 1
          3, & ! row 2
          5, & ! row 3
          7, & ! row 4
          10, & ! row 5
          13, & ! row 6
          16, & ! row 7
          19, & ! row 8
          21  /)
        integer :: i
        integer :: itr_max          
        integer :: mr
        real(dp), dimension(n) :: rhs = (/ &
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
        real(dp) :: x_error, soln_error
        real(dp) :: x_estimate(n)
        real(dp), dimension(n) :: x_exact = (/ &
          3D+00, &
          1D+00, &
          1D+00, &
          5D+00, &
          6D+00, &
          6D+00, &
          5D+00, &
          3D+00 /)
        integer :: ua(n)
        real(dp) :: l(nz_num+2) ! (ia(n+1)+1)
        real(dp) :: Ax_for_soln(n)
        
        include 'formats'
        
        if (nz_num+2 /= ia(n+1)+1) then
           stop 'bad counting test_ILU_CR_mgmres'
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'test_ILU_CR_mgmres_nvar2'

        do test = 1, 1 ! 2

          if ( test == 1 ) then

            x_estimate(1:n) = 0.0D+00

          else

            write ( *, '(a)' ) '  Second try, use random initial vector:'

            call r8vec_uniform_01 ( n, seed, x_estimate )

          end if

          x_error = sqrt ( sum ( ( x_exact(1:n) - x_estimate(1:n) )**2 ) )

          write ( *, '(a,g14.6)' ) '  Before solving, X_ERROR = ', x_error

          itr_max = 20
          mr = n - 1
          tol_abs = 1.0D-08
          tol_rel = 1.0D-08
                    
         call create_preconditioner_ILU_CR_mgmres()            
         !write(*,*) ' call mgmres'
         call mgmres ( &     
            n, matvec_ILU_CR_mgmres, apply_ILU_CR_mgmres, &
            x_estimate, rhs, itr_max, mr, tol_abs, tol_rel )
         call ax_cr ( n, nz_num, ia, ja, a, x_estimate, Ax_for_soln )

          soln_error = sqrt ( sum ( ( rhs(1:n) - Ax_for_soln(1:n) )**2 ) )
          x_error = sqrt ( sum ( ( x_exact(1:n) - x_estimate(1:n) )**2 ) )

          write ( *, '(a,g14.6)' ) '  x=soln error ', x_error
          write ( *, '(a)' ) '  x:'
          do i = 1, n
            write ( *, '(2x,i8,2x,g14.6)' ) i, x_estimate(i)
          end do

        end do

        contains
        
        subroutine create_preconditioner_ILU_CR_mgmres()
           !call rearrange_cr ( n, nz_num, ia, ja, a )
           call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )
           call ilu_cr ( n, nz_num, ia, ja, a, ua, l ) 
        end subroutine create_preconditioner_ILU_CR_mgmres
     
        subroutine apply_ILU_CR_mgmres(x) ! apply preconditioner to x
           real(dp), intent(inout) :: x(:) ! (neq)
           call lus_cr ( n, nz_num, ia, ja, l, ua, x, x ) 
        end subroutine apply_ILU_CR_mgmres

        subroutine matvec_ILU_CR_mgmres(x, r) ! set r = Jacobian*x
           real(dp), intent(in) :: x(:) ! (neq)
           real(dp), intent(out) :: r(:) ! (neq)
           call ax_cr ( n, nz_num, ia, ja, a, x, r ) 
        end subroutine matvec_ILU_CR_mgmres
        
      end subroutine test_ILU_CR_mgmres_nvar2
      

!*****************************************************************************
!*****************************************************************************
!
!  support routines for sparse arrays
!
!*****************************************************************************
!*****************************************************************************

      
      subroutine atx_cr ( n, nz_num, ia, ja, a, x, w )

      !*****************************************************************************
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
      !    Input, integer :: N, the order of the system.
      !
      !    Input, integer :: NZ_NUM, the number of nonzeros.
      !
      !    Input, integer :: IA(N+1), JA(NZ_NUM), the row and column
      !    indices of the matrix values.  The row vector has been compressed.
      !
      !    Input, real(dp) :: A(NZ_NUM), the matrix values.
      !
      !    Input, real(dp) :: X(N), the vector to be multiplied by A'.
      !
      !    Output, real(dp) :: W(N), the value of A'*X.
      !
        implicit none

        integer :: n
        integer :: nz_num

        real(dp) :: a(nz_num)
        integer :: i
        integer :: ia(n+1)
        integer :: ja(nz_num)
        integer :: k1
        integer :: k2
        real(dp) :: w(n)
        real(dp) :: x(n)

        w(1:n) = 0.0D+00

        do i = 1, n
          k1 = ia(i)
          k2 = ia(i+1) - 1
          w(ja(k1:k2)) = w(ja(k1:k2)) + a(k1:k2) * x(i)
        end do

        return
      end subroutine atx_cr
      
      subroutine atx_st ( n, nz_num, ia, ja, a, x, w )

      !*****************************************************************************
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
      !    Input, integer :: N, the order of the system.
      !
      !    Input, integer :: NZ_NUM, the number of nonzeros.
      !
      !    Input, integer :: IA(NZ_NUM), JA(NZ_NUM), the row and column
      !    indices of the matrix values.
      !
      !    Input, real(dp) :: A(NZ_NUM), the matrix values.
      !
      !    Input, real(dp) :: X(N), the vector to be multiplied by A'.
      !
      !    Output, real(dp) :: W(N), the value of A'*X.
      !
        implicit none

        integer :: n
        integer :: nz_num

        real(dp) :: a(nz_num)
        integer :: i
        integer :: ia(nz_num)
        integer :: j
        integer :: ja(nz_num)
        integer :: k
        real(dp) :: w(n)
        real(dp) :: x(n)

        w(1:n) = 0.0D+00

        do k = 1, nz_num
          i = ia(k)
          j = ja(k)
          w(j) = w(j) + a(k) * x(i)
        end do

        return
      end subroutine atx_st
      
      subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )

      !*****************************************************************************
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
      !    Input, integer :: N, the order of the system.
      !
      !    Input, integer :: NZ_NUM, the number of nonzeros.
      !
      !    Input, integer :: IA(N+1), JA(NZ_NUM), the row and column
      !    indices of the matrix values.  The row vector has been compressed.
      !
      !    Input, real(dp) :: A(NZ_NUM), the matrix values.
      !
      !    Input, real(dp) :: X(N), the vector to be multiplied by A.
      !
      !    Output, real(dp) :: W(N), the value of A*X.
      !
        implicit none

        integer :: n
        integer :: nz_num

        real(dp) :: a(nz_num)
        integer :: i
        integer :: ia(n+1)
        integer :: ja(nz_num)
        integer :: k1
        integer :: k2
        real(dp) :: w(n)
        real(dp) :: x(n)

        w(1:n) = 0.0D+00

        do i = 1, n
          k1 = ia(i)
          k2 = ia(i+1) - 1
          w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
        end do

        return
      end subroutine ax_cr
      
      subroutine ax_st ( n, nz_num, ia, ja, a, x, w )

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
        implicit none

        integer :: n
        integer :: nz_num

        real(dp) :: a(nz_num)
        integer :: i
        integer :: ia(nz_num)
        integer :: j
        integer :: ja(nz_num)
        integer :: k
        real(dp) :: w(n)
        real(dp) :: x(n)

        w(1:n) = 0.0D+00

        do k = 1, nz_num
          i = ia(k)
          j = ja(k)
          w(i) = w(i) + a(k) * x(j)
        end do

        return
      end subroutine ax_st
      
      subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

      !*****************************************************************************
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
      !    Input, integer :: N, the order of the system.
      !
      !    Input, integer :: NZ_NUM, the number of nonzeros.
      !
      !    Input, integer :: IA(N+1), JA(NZ_NUM), the row and column
      !    indices of the matrix values.  The row vector has been compressed.
      !    On output, the order of the entries of JA may have changed because of
      !    the sorting.
      !
      !    Output, integer :: UA(N), the index of the diagonal element
      !    of each row.
      !
        implicit none

        integer :: n
        integer :: nz_num

        integer :: i
        integer :: ia(n+1)
        integer :: k
        integer :: ja(nz_num)
        integer :: ua(n)

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

      !*****************************************************************************
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
      !    Input, integer :: N, the order of the system.
      !
      !    Input, integer :: NZ_NUM, the number of nonzeros.
      !
      !    Input, integer :: IA(N+1), JA(NZ_NUM), the row and column
      !    indices of the matrix values.  The row vector has been compressed.
      !
      !    Input, real(dp) :: A(NZ_NUM), the matrix values.
      !
      !    Input, integer :: UA(N), the index of the diagonal element
      !    of each row.
      !
      !    Output, real(dp) :: L(NZ_NUM), the ILU factorization of A.
      !
        implicit none

        integer :: n
        integer :: nz_num

        real(dp) :: a(nz_num)
        integer :: i
        integer :: ia(n+1)
        integer :: iw(n)
        integer :: j
        integer :: ja(nz_num)
        integer :: jj
        integer :: jrow
        integer :: jw
        integer :: k
        real(dp) :: l(nz_num)
        real(dp) :: tl
        integer :: ua(n)
        include 'formats'
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
        
        write(*,1) 'sparse LU Dhat', l(ua(1:n))

        return
      end subroutine ilu_cr
      
      subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )

      !*****************************************************************************
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
      !    Input, integer :: N, the order of the system.
      !
      !    Input, integer :: NZ_NUM, the number of nonzeros.
      !
      !    Input, integer :: IA(N+1), JA(NZ_NUM), the row and column
      !    indices of the matrix values.  The row vector has been compressed.
      !
      !    Input, real(dp) :: L(NZ_NUM), the matrix values.
      !
      !    Input, integer :: UA(N), the index of the diagonal element
      !    of each row.
      !
      !    Input, real(dp) :: R(N), the right hand side.
      !
      !    Output, real(dp) :: Z(N), the solution of the system M * Z = R.
      !
        implicit none

        integer :: n
        integer :: nz_num

        integer :: i
        integer :: ia(n+1)
        integer :: j
        integer :: ja(nz_num)
        real(dp) :: l(nz_num)
        real(dp) :: r(n)
        integer :: ua(n)
        real(dp) :: w(n)
        real(dp) :: z(n)
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
      
      subroutine r8vec_uniform_01 ( n, seed, r )

      !*****************************************************************************
      !
      !! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
      !
      !  Discussion:
      !
      !    An R8VEC is a vector of real(dp) :: values.
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
      !    Input, integer :: N, the number of entries in the vector.
      !
      !    Input/output, integer :: SEED, the "seed" value, which
      !    should NOT be 0.  On output, SEED has been updated.
      !
      !    Output, real(dp) :: R(N), the vector of pseudorandom values.
      !
        implicit none

        integer :: n

        integer :: i
        integer :: k
        integer :: seed
        real(dp) :: r(n)

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

      !*****************************************************************************
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
      !    Input, integer :: N, the order of the system.
      !
      !    Input, integer :: NZ_NUM, the number of nonzeros.
      !
      !    Input, integer :: IA(N+1), the compressed row indices.
      !
      !    Input/output, integer :: JA(NZ_NUM), the column indices.
      !    On output, these may have been rearranged by the sorting.
      !
      !    Input/output, real(dp) :: A(NZ_NUM), the matrix values.  On output,
      !    the matrix values may have been moved somewhat because of the sorting.
      !
        implicit none

        integer :: n
        integer :: nz_num

        real(dp) :: a(nz_num)
        integer :: i
        integer :: ia(n+1)
        integer :: i4temp
        integer :: ja(nz_num)
        integer :: k
        integer :: l
        real(dp) :: r8temp

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

      
!*****************************************************************************
!*****************************************************************************
!
!  matrix support routines
!
!*****************************************************************************
!*****************************************************************************
      
      
      subroutine BTD_mv1(nvar, nz, lblk, dblk, ublk, b1, r1)
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
      end subroutine BTD_mv1                  
      
      subroutine m_inverse(k, nvar, blk, blk_inv, ipiv, ierr)
         use star_bcyclic, only: my_getf2
         integer, intent(in) :: k, nvar
         real(dp), intent(in) :: blk(:,:) ! (nvar,nvar)
         real(dp), intent(out) :: blk_inv(:,:) ! (nvar,nvar)
         integer, intent(out) :: ipiv(:), ierr
         include 'formats'
         real(dp) :: work(nvar), binv(nvar,nvar)
         integer :: i, j, ip(nvar)
         blk_inv = blk
         call my_getf2(nvar, blk_inv, nvar, ipiv, ierr)
         !call DGETRF(nvar, nvar, blk_inv, nvar, ipiv, ierr)
         if (ierr /= 0) then
            write(*,3) 'DGETRF failed', k, ierr
            stop 'create_pcond'
         end if
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
         if (ierr /= 0) then
            write(*,3) 'DGETRI failed', k, ierr
            stop 'create_pcond'
         end if         
         do j=1,nvar
            !$omp simd
            do i=1,nvar
               blk_inv(i,j) = binv(i,j)
            end do
         end do
      end subroutine m_inverse
         
      subroutine mm_0(a, b, c) ! c := a*b, c different than b
         real(dp), dimension(:,:) :: a, b, c ! (nvar,nvar)
         integer :: j, i, nvar
         nvar=size(a,dim=1)
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

      subroutine mv_self(a,x) ! x = a*x
         real(dp) :: a(:,:) ! (nvar,nvar)
         real(dp) :: x(:) ! (nvar)
         real(dp) :: tmp, temp(size(x,dim=1))
         integer :: j, nvar
         nvar=size(a,dim=1)
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
      subroutine mgmres ( n, matvec, psolve, x, rhs, itr_max, mr, &
              tol_abs, tol_rel )
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
        integer :: n
        real(dp) :: x(n)
        real(dp) :: rhs(n)
        integer :: itr_max
        integer :: mr
        real(dp) :: tol_abs
        real(dp) :: tol_rel
         
        ! locals
        real(dp) :: r(n)
        !integer :: ua(n)
        real(dp) :: v(n,mr+1);

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
        logical, parameter :: verbose = .true.
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
      
      
      end module abtilu
