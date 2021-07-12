! ***********************************************************************
! Copyright (C) 2012  The MESA Team
! This file is part of MESA.
! MESA is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Library Public License as published
! by the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
! MESA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General Public License
! along with this software; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
! ***********************************************************************

! derived from BCYCLIC written hirshman et. al.
! S.P.Hirshman, K.S.Perumalla, V.E.Lynch, & R.Sanchez,
! BCYCLIC: A parallel block tridiagonal matrix cyclic solver,
! J. Computational Physics, 229 (2010) 6392-6404.


      module star_bcyclic

      use star_private_def
      use const_def, only: dp, ln10
      use utils_lib, only: set_nan

      implicit none

      private
      public :: bcyclic_factor, bcyclic_solve, clear_storage, &
         my_getf2_n4, my_getf2_n5, my_getf2, my_gemm0_p1, &
         my_getrs1_n4, my_getrs1_n5, my_getrs1, my_gemv_p1

      logical, parameter :: dbg = .false.
      logical, parameter :: do_set_nan = .false.

      logical, parameter :: blas3_dgemm3 = .false.

      contains

      subroutine bcyclic_factor ( &
            s, nvar, nz, lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipivot1, &
            B1, row_scale_factors1, col_scale_factors1, &
            equed1, iter, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar ! linear size of each block
         integer, intent(in) :: nz ! number of block rows
         real(dp), pointer, dimension(:) :: &
            lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, &
            B1, row_scale_factors1, col_scale_factors1
         integer, pointer :: ipivot1(:)
         character (len=nz) :: equed1
         integer, intent(in) :: iter ! solver iteration number for debugging output
         integer, intent(out) :: ierr

         integer, pointer :: iptr(:,:), nslevel(:), ipivot(:)
         integer :: neq, ncycle, nstemp, maxlevels, nlevel, i, j, k
         logical :: have_odd_storage
         real(dp), pointer, dimension(:,:) :: dmat, dmatF
         real(dp), pointer, dimension(:) :: row_scale_factors, col_scale_factors
         character (len=1) :: equed
         real(dp) :: min_rcond_from_DGESVX, rpgfac
         integer :: k_min_rcond_from_DGESVX
         
         integer, allocatable :: factored(:)

         include 'formats'
         
         if (s% use_DGESVX_in_bcyclic .and. s% report_min_rcond_from_DGESXV) &
            min_rcond_from_DGESVX = 1d99
         
         allocate(factored(nz))
         do k=1,nz
            factored(k) = 0
         end do

         ierr = 0
         neq = nvar*nz
         !$omp simd
         do i = 1,nvar*neq
            lblkF1(i) = lblk1(i)
            dblkF1(i) = dblk1(i)
            ublkF1(i) = ublk1(i)
         end do

         if (dbg) write(*,*) 'start bcyclic_factor'

         ! compute number of cyclic reduction levels
         ncycle = 1
         maxlevels = 0
         do while (ncycle < nz)
            ncycle = 2*ncycle
            maxlevels = maxlevels+1
         end do
         maxlevels = max(1, maxlevels)

         have_odd_storage = associated(s% bcyclic_odd_storage)
         if (have_odd_storage) then
            if (size(s% bcyclic_odd_storage) < maxlevels) then
               call clear_storage(s)
               have_odd_storage = .false.
            end if
         end if

         if (.not. have_odd_storage) then
            allocate (s% bcyclic_odd_storage(maxlevels+3), stat=ierr)
            if (ierr /= 0) then
               write(*,*) 'alloc failed for odd_storage in bcyclic'
               return
            end if
            do nlevel = 1, size(s% bcyclic_odd_storage)
               s% bcyclic_odd_storage(nlevel)% ul_size = 0
            end do
         end if

         allocate (nslevel(maxlevels), stat=ierr)
         if (ierr /= 0) return

         ncycle = 1
         nstemp = nz
         nlevel = 1

         if (dbg) write(*,*) 'start factor_cycle'

         factor_cycle: do ! perform cyclic-reduction factorization

            nslevel(nlevel) = nstemp

            if (dbg) write(*,2) 'call cycle_onestep', nstemp

            call cycle_onestep( &
               s, nvar, nz, nstemp, ncycle, nlevel, iter, &
               lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipivot1, &
               B1, row_scale_factors1, col_scale_factors1, equed1, factored, &
               min_rcond_from_DGESVX, k_min_rcond_from_DGESVX, rpgfac, &
               ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if

            if (nstemp == 1) exit

            nstemp = (nstemp+1)/2
            nlevel = nlevel+1
            ncycle = 2*ncycle

            if (nlevel > maxlevels) exit

         end do factor_cycle

         if (dbg) write(*,*) 'done factor_cycle'

         ! factor row 1
         dmat(1:nvar,1:nvar) => dblk1(1:nvar*nvar)
         dmatF(1:nvar,1:nvar) => dblkF1(1:nvar*nvar)
         ipivot(1:nvar) => ipivot1(1:nvar)
         row_scale_factors(1:nvar) => row_scale_factors1(1:nvar)
         col_scale_factors(1:nvar) => col_scale_factors1(1:nvar)
         factored(1) = factored(1) + 1
         call dense_factor(s, 1, nvar, dmat, dmatF, ipivot, &
            row_scale_factors, col_scale_factors, equed, &
            min_rcond_from_DGESVX, k_min_rcond_from_DGESVX, rpgfac, &
            ierr)
         equed1(1:1) = equed(1:1)
         if (ierr /= 0) then
            write(*,*) 'dense_factor failed'
            call dealloc
            return
         end if
         
         do k=1,nz ! check that every cell factored exactly once
            if (factored(k) /= 1) then
               write(*,3) 'factored /= 1', k, factored(k)
               stop 'bcyclic_factor'
            end if
         end do

         call dealloc
            
         if (s% use_DGESVX_in_bcyclic .and. s% report_min_rcond_from_DGESXV) then
            write(*,4) 'DGESVX: k_min, iter, model, min rcond, rpgfac', &
               k_min_rcond_from_DGESVX, iter, s% model_number, min_rcond_from_DGESVX, rpgfac
         end if

         if (dbg) write(*,*) 'done bcyclic_factor'

         contains

         subroutine dealloc
            deallocate (nslevel)
         end subroutine dealloc


      end subroutine bcyclic_factor


      subroutine cycle_onestep( &
            s, nvar, nz, nblk, ncycle, nlevel, iter, &
            lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipivot1, &
            B1, row_scale_factors1, col_scale_factors1, equed1, factored, &
            min_rcond_from_DGESVX, k_min_rcond_from_DGESVX, rpgfac, &
            ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar, nz, nblk, ncycle, nlevel, iter
         real(dp), pointer, dimension(:), intent(inout) :: &
            lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, &
            B1, row_scale_factors1, col_scale_factors1
         character (len=nz) :: equed1
         integer, pointer, intent(inout) :: ipivot1(:)
         integer :: factored(:)
         real(dp) :: min_rcond_from_DGESVX, rpgfac
         integer :: k_min_rcond_from_DGESVX
         integer, intent(out) :: ierr

         integer, pointer :: ipivot(:)
         real(dp), pointer, dimension(:,:) :: dmat, umat, lmat, umat0, lmat0, dmatF
         real(dp), pointer, dimension(:,:) :: lnext, unext, lprev, uprev
         real(dp), pointer, dimension(:) :: mat1
         integer :: i, j, shift, min_sz, new_sz, shift1, shift2, nvar2, &
            ns, op_err, nmin, kcount, k, ii, jj, kk
         real(dp), pointer, dimension(:) :: row_scale_factors, col_scale_factors
         character (len=1) :: equed

         include 'formats'

         ierr = 0
         nvar2 = nvar*nvar
         nmin = 1
         kcount = 1+(nblk-nmin)/2
         min_sz = nvar2*kcount
         if (s% bcyclic_odd_storage(nlevel)% ul_size < min_sz) then
            if (s% bcyclic_odd_storage(nlevel)% ul_size > 0) &
               deallocate( &
                  s% bcyclic_odd_storage(nlevel)% umat1, &
                  s% bcyclic_odd_storage(nlevel)% lmat1)
            new_sz = min_sz*1.1 + 100
            s% bcyclic_odd_storage(nlevel)% ul_size = new_sz
            allocate (s% bcyclic_odd_storage(nlevel)% umat1(new_sz), &
                      s% bcyclic_odd_storage(nlevel)% lmat1(new_sz), stat=ierr)
            if (ierr /= 0) then
               write(*,*) 'allocation error in cycle_onestep'
               return
            end if
         end if

!$OMP PARALLEL DO private(ns,kcount,shift,shift2,i) SCHEDULE(static,3)
         do ns = nmin, nblk, 2  ! copy umat and lmat
            kcount = (ns-nmin)/2 + 1
            shift = nvar2*(kcount-1)
            shift2 = nvar2*ncycle*(ns-1)
            do i=1,nvar2
               s% bcyclic_odd_storage(nlevel)% umat1(shift+i) = ublkF1(shift2+i)
               s% bcyclic_odd_storage(nlevel)% lmat1(shift+i) = lblkF1(shift2+i)
            end do
         end do
!$OMP END PARALLEL DO

         if (nvar2*kcount > s% bcyclic_odd_storage(nlevel)% ul_size) then
            write(*,*) 'nvar2*kcount > ul_size in cycle_onestep'
            ierr = -1
            return
         end if

         if (dbg) write(*,*) 'start lu factorization'
         ! compute lu factorization of even diagonal blocks
         nmin = 2
!$OMP PARALLEL DO SCHEDULE(static,3) &
!$OMP PRIVATE(ipivot,dmat,dmatF,ns,op_err,shift1,shift2,i,j,k,row_scale_factors,col_scale_factors,equed)
         do ns = nmin, nblk, 2

            k = ncycle*(ns-1) + 1
            shift1 = nvar*(k-1)
            shift2 = nvar*shift1
            dmat(1:nvar,1:nvar) => dblk1(shift2+1:shift2+nvar2)
            dmatF(1:nvar,1:nvar) => dblkF1(shift2+1:shift2+nvar2)
            op_err = 0
            ipivot(1:nvar) => ipivot1(shift1+1:shift1+nvar)
            row_scale_factors(1:nvar) => row_scale_factors1(shift1+1:shift1+nvar)
            col_scale_factors(1:nvar) => col_scale_factors1(shift1+1:shift1+nvar)
            factored(k) = factored(k) + 1
            call dense_factor(s, k, nvar, dmat, dmatF, ipivot, &
               row_scale_factors, col_scale_factors, equed, &
               min_rcond_from_DGESVX, k_min_rcond_from_DGESVX, rpgfac, &
               op_err)
            equed1(k:k) = equed(1:1)
            if (op_err /= 0) then
               ierr = op_err
            end if

         end do
!$OMP END PARALLEL DO
         if (ierr /= 0) then
            !write(*,*) 'factorization failed in bcyclic'
            return
         end if

         if (dbg) write(*,*) 'done lu factorization; start solve'

!$OMP PARALLEL DO SCHEDULE(static,3) &
!$OMP PRIVATE(ns,k,shift1,shift2,ipivot,dmat,dmatF,umat,lmat,mat1,i,j,row_scale_factors,col_scale_factors,equed,op_err)
         do ns = nmin, nblk, 2
            ! compute new l=-d[-1]l, u=-d[-1]u for even blocks
            k = ncycle*(ns-1) + 1
            shift1 = nvar*(k-1)
            shift2 = nvar*shift1

            lmat(1:nvar,1:nvar) => lblkF1(shift2+1:shift2+nvar2)

            dmat(1:nvar,1:nvar) => dblk1(shift2+1:shift2+nvar2)
            dmatF(1:nvar,1:nvar) => dblkF1(shift2+1:shift2+nvar2)
            ipivot(1:nvar) => ipivot1(shift1+1:shift1+nvar)
            row_scale_factors(1:nvar) => row_scale_factors1(shift1+1:shift1+nvar)
            col_scale_factors(1:nvar) => col_scale_factors1(shift1+1:shift1+nvar)
            equed(1:1) = equed1(k:k)
            call dense_solve(s, k, nvar, dmat, dmatF, ipivot, lmat, &
               row_scale_factors, col_scale_factors, equed, op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if

            do j=1,nvar
               do i=1,nvar
                  lmat(i,j) = -lmat(i,j)
               end do
            end do

            umat(1:nvar,1:nvar) => ublkF1(shift2+1:shift2+nvar2)

            call dense_solve(s, k, nvar, dmat, dmatF, ipivot, umat, &
               row_scale_factors, col_scale_factors, equed, op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if

            do j=1,nvar
               do i=1,nvar
                  umat(i,j) = -umat(i,j)
               end do
            end do

         end do
!$OMP END PARALLEL DO
         if (dbg) write(*,*) 'done solve'

         if (ierr /= 0) return

         ! compute new odd blocks in terms of even block factors
         ! compute odd hatted matrix elements except at boundaries
         nmin = 1
!$OMP PARALLEL DO SCHEDULE(static,3) &
!$OMP PRIVATE(i,ns,shift2,dmat,umat,lmat,lnext,unext,lprev,uprev,kcount,shift,umat0,lmat0,k)
         do i= 1, 3*(1+(nblk-nmin)/2)

            ns = 2*((i-1)/3) + nmin
            k = ncycle*(ns-1) + 1
            if (factored(k) > 0) then
               write(*,2) 'compute new dmat after already factored', k
               stop 'cycle_onestep'
            end if
            shift2 = nvar2*(k-1)
            dmat(1:nvar,1:nvar) => dblkF1(shift2+1:shift2+nvar2)
            umat(1:nvar,1:nvar) => ublkF1(shift2+1:shift2+nvar2)
            lmat(1:nvar,1:nvar) => lblkF1(shift2+1:shift2+nvar2)

            if (ns < nblk) then
               shift2 = nvar2*ncycle*ns
               lnext(1:nvar,1:nvar) => lblkF1(shift2+1:shift2+nvar2)
               unext(1:nvar,1:nvar) => ublkF1(shift2+1:shift2+nvar2)
            end if

            if (ns > 1) then
               shift2 = nvar2*ncycle*(ns-2)
               lprev(1:nvar,1:nvar) => lblkF1(shift2+1:shift2+nvar2)
               uprev(1:nvar,1:nvar) => ublkF1(shift2+1:shift2+nvar2)
            end if

            kcount = 1+(ns-nmin)/2
            shift = nvar2*(kcount-1)
            lmat0(1:nvar,1:nvar) => &
               s% bcyclic_odd_storage(nlevel)% lmat1(shift+1:shift+nvar2)
            umat0(1:nvar,1:nvar) => &
               s% bcyclic_odd_storage(nlevel)% umat1(shift+1:shift+nvar2)

            select case(mod(i-1,3))
            case (0)
               if (ns > 1) then
                  ! lmat = matmul(lmat0, lprev)
                  if (blas3_dgemm3) then
                     call my_gemm3_0_p1(nvar,nvar,nvar,lmat0,lprev,lmat)
                  else
                     call my_gemm0_p1(nvar,nvar,nvar,lmat0,nvar,lprev,nvar,lmat,nvar)
                  end if
               end if
            case (1)
               if (ns < nblk) then
                  ! umat = matmul(umat0, unext)
                  if (blas3_dgemm3) then
                     call my_gemm3_0_p1(nvar,nvar,nvar,umat0,unext,umat)
                  else
                     call my_gemm0_p1(nvar,nvar,nvar,umat0,nvar,unext,nvar,umat,nvar)
                  end if
               end if
            case (2)
               if (ns < nblk) then
                  if (ns > 1) then
                     ! dmat = dmat + matmul(umat0, lnext) + matmul(lmat0,uprev)
                     if (blas3_dgemm3) then ! my_gemm3_plus_mm BAD
                        call my_gemm3_plus_mm(nvar,nvar,nvar,umat0,lnext,lmat0,uprev,dmat)
                     else
                        call my_gemm_plus_mm(nvar,nvar,nvar,umat0,lnext,lmat0,uprev,dmat)
                     end if
                  else
                     ! dmat = dmat + matmul(umat0, lnext)
                     if (blas3_dgemm3) then
                        call my_gemm3_p1(nvar,nvar,nvar,umat0,lnext,dmat)
                     else
                        call my_gemm_p1(nvar,nvar,nvar,umat0,nvar,lnext,nvar,dmat,nvar)
                     end if
                  end if
               else if (ns > 1) then
                  ! dmat = dmat + matmul(lmat0,uprev)
                  if (blas3_dgemm3) then
                     call my_gemm3_p1(nvar,nvar,nvar,lmat0,uprev,dmat)
                  else
                     call my_gemm_p1(nvar,nvar,nvar,lmat0,nvar,uprev,nvar,dmat,nvar)
                  end if
               end if
            end select

         end do
!$OMP END PARALLEL DO
         if (dbg) write(*,*) 'done cycle_onestep'

      end subroutine cycle_onestep


      subroutine cycle_rhs( &
            s, nz, nblk, nvar, ncycle, nlevel, &
            dblk1, dblkF1, soln1, ipivot1, &
            row_scale_factors1, col_scale_factors1, equed1, ierr)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         integer, intent(in) :: nz, nblk, nvar, ncycle, nlevel
         real(dp), pointer, intent(in), dimension(:) :: &
            dblk1, dblkF1, row_scale_factors1, col_scale_factors1
         real(dp), pointer, intent(inout) :: soln1(:)
         integer, pointer, intent(in) :: ipivot1(:)
         character (len=nz) :: equed1
         integer, intent(out) :: ierr

         integer :: i, k, ns, op_err, nmin, kcount, shift, shift1, shift2, nvar2
         integer, pointer :: ipivot(:)
         real(dp), pointer, dimension(:,:) :: dmatF, dmat, umat, lmat
         real(dp), pointer, dimension(:) :: X, Xprev, Xnext
         real(dp), pointer, dimension(:) :: row_scale_factors, col_scale_factors
         character (len=1) :: equed
         logical :: okay

         include 'formats'

         ierr = 0
         nvar2 = nvar*nvar
         ! compute dblk[-1]*brhs for even indices and store in brhs(even)
         nmin = 2
         op_err = 0
!$OMP PARALLEL DO SCHEDULE(static,3) &
!$OMP PRIVATE(ns,shift1,ipivot,shift2,k,dmat,dmatF,X,row_scale_factors,col_scale_factors,equed,i,okay,op_err)
         do ns = nmin, nblk, 2
            k = ncycle*(ns-1) + 1
            shift1 = nvar*(k-1)
            shift2 = nvar*shift1
            dmat(1:nvar,1:nvar) => dblk1(shift2+1:shift2+nvar2)
            dmatF(1:nvar,1:nvar) => dblkF1(shift2+1:shift2+nvar2)
            ipivot(1:nvar) => ipivot1(shift1+1:shift1+nvar)
            row_scale_factors(1:nvar) => row_scale_factors1(shift1+1:shift1+nvar)
            col_scale_factors(1:nvar) => col_scale_factors1(shift1+1:shift1+nvar)
            equed(1:1) = equed1(k:k)
            X(1:nvar) => soln1(shift1+1:shift1+nvar)
            call dense_solve1(s, k, nvar, X, dmat, dmatF, ipivot, .true., &
               row_scale_factors, col_scale_factors, equed, op_err)
            if (op_err /= 0) then
               ierr = op_err
               cycle
            end if

         end do
!$OMP END PARALLEL DO

        if (ierr /= 0) return

        ! compute odd (hatted) sources (b-hats) for interior rows
         nmin = 1
         kcount = 0
!$OMP PARALLEL DO SCHEDULE(static,3) &
!$OMP PRIVATE(ns,shift1,X,kcount,shift,umat,lmat,Xnext,Xprev)
         do ns = nmin, nblk, 2
            shift1 = nvar*ncycle*(ns-1)
            X(1:nvar) => soln1(shift1+1:shift1+nvar)
            kcount = 1+(ns-nmin)/2
            shift = nvar2*(kcount-1)
            umat(1:nvar,1:nvar) => &
               s% bcyclic_odd_storage(nlevel)% umat1(shift+1:shift+nvar2)
            lmat(1:nvar,1:nvar) => &
               s% bcyclic_odd_storage(nlevel)% lmat1(shift+1:shift+nvar2)
            if (ns > 1) then
               shift1 = nvar*ncycle*(ns-2)
               Xprev => soln1(shift1+1:shift1+nvar)
            end if
            if (ns < nblk) then
               shift1 = nvar*ncycle*ns
               Xnext => soln1(shift1+1:shift1+nvar)
               if (ns > 1) then
                  ! bptr = bptr - matmul(umat,bnext) - matmul(lmat,bprev)
                  call my_gemv_mv(nvar,nvar,umat,Xnext,lmat,Xprev,X)
               else
                  ! bptr = bptr - matmul(umat,bnext)
                  call my_gemv(nvar,nvar,umat,nvar,Xnext,X)
               end if
            else if (ns > 1) then
               ! bptr = bptr - matmul(lmat,bprev)
               call my_gemv(nvar,nvar,lmat,nvar,Xprev,X)
            end if
         end do
!$OMP END PARALLEL DO

         if (nvar2*kcount > s% bcyclic_odd_storage(nlevel)% ul_size) then
            write(*,*) 'nvar2*kcount > ul_size in cycle_rhs'
            ierr = -1
            return
         end if

      end subroutine cycle_rhs


      ! computes even index solution from the computed (at previous,higher level)
      ! odd index solutions at this level.
      ! note at this point, the odd brhs values have been replaced (at the highest cycle)
      ! with the solution values (x), at subsequent (lower) cycles, the
      ! odd values are replaced by the even solutions at the next highest cycle. the even
      ! brhs values were multiplied by d[-1] and stored in cycle_rhs
      ! solve for even index values in terms of (computed at this point) odd index values
      subroutine cycle_solve( &
            s, nvar, nz, ncycle, nblk, nlevel, lblk1, ublk1, lblkF1, ublkF1, soln1)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar, nz, ncycle, nblk, nlevel
         real(dp), pointer, intent(in), dimension(:) :: lblk1, ublk1, lblkF1, ublkF1
         real(dp), pointer, intent(inout) :: soln1(:)

         real(dp), pointer :: umat(:,:), lmat(:,:), bprev(:), bnext(:), bptr(:)
         real(dp), pointer, dimension(:) :: bprevr, bnextr
         integer :: shift1, shift2, nvar2, ns, ierr, nmin, i, j

         include 'formats'

         nvar2 = nvar*nvar
         nmin = 2
!$OMP PARALLEL DO SCHEDULE(static,3) &
!$OMP PRIVATE(ns,shift1,bptr,shift2,lmat,bprev,umat,bnext)
         do ns = nmin, nblk, 2
            shift1 = ncycle*nvar*(ns-1)
            bptr(1:nvar) => soln1(shift1+1:shift1+nvar)
            shift2 = nvar*shift1
            lmat(1:nvar,1:nvar) => lblkF1(shift2+1:shift2+nvar2)
            if (ns > 1) then
               shift1 = ncycle*nvar*(ns-2)
               bprev(1:nvar) => soln1(shift1+1:shift1+nvar)
            end if
            if (ns < nblk) then
               umat(1:nvar,1:nvar) => ublkF1(shift2+1:shift2+nvar2)
               shift1 = ncycle*nvar*ns
               bnext(1:nvar) => soln1(shift1+1:shift1+nvar)
               if (ns > 1) then
                  ! bptr = bptr + matmul(umat,bnext) + matmul(lmat,bprev)
                  call my_gemv_p_mv(nvar,nvar,umat,bnext,lmat,bprev,bptr)
               else
                  ! bptr = bptr + matmul(umat,bnext)
                  call my_gemv_p1(nvar,nvar,umat,nvar,bnext,bptr)
               end if
            else if (ns > 1) then
               ! bptr = bptr + matmul(lmat,bprev)
               call my_gemv_p1(nvar,nvar,lmat,nvar,bprev,bptr)
            end if
         end do
!$OMP END PARALLEL DO

      end subroutine cycle_solve


      subroutine dense_factor(s, k, nvar, mtx, mtxF, ipivot, &
            row_scale_factors, col_scale_factors, equed, &
            min_rcond_from_DGESVX, k_min_rcond_from_DGESVX, rpgfac, &
            ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         real(dp), pointer :: mtx(:,:), mtxF(:,:)
         integer, pointer :: ipivot(:)
         real(dp), pointer :: row_scale_factors(:), col_scale_factors(:)
         character (len=1) :: equed
         real(dp) :: min_rcond_from_DGESVX, rpgfac
         integer :: k_min_rcond_from_DGESVX
         integer, intent(out) :: ierr
         logical :: singular
         integer :: i, j
         real(dp), pointer :: work(:)
         integer, pointer :: iwork(:)
         real(dp) :: anorm, rcond
         include 'formats'
         ierr = 0
         
         if (s% use_DGESVX_in_bcyclic) then
            call factor_with_DGESVX
            return
         end if
         
         if (nvar == 4) then
            call my_getf2_n4(mtxF, ipivot, ierr)
         else if (nvar == 5) then
            call my_getf2_n5(mtxF, ipivot, ierr)
         else
            call my_getf2(nvar, mtxF, nvar, ipivot, ierr)
         end if
         
         contains
         
         subroutine factor_with_DGESVX
            character (len=1) :: fact, trans
            integer, parameter :: nrhs = 0
            real(dp) :: rcond
            real(dp) :: a(nvar,nvar), af(nvar,nvar), b(nvar,nrhs), x(nvar,nrhs), &
               r(nvar), c(nvar), ferr(nrhs), berr(nrhs), work(4*nvar)
            integer :: ipiv(nvar), iwork(nvar)
            integer :: i, j
            include 'formats'

            do i=1,nvar
               do j=1,nvar
                  a(i,j) = mtxF(i,j)
               end do
            end do
            
            if (s% use_equilibration_in_DGESVX) then
               fact = 'E' ! matrix A will be equilibrated, then copied to AF and factored
            else
               fact = 'N' ! matrix A will be copied to AF and factored
            end if
            trans = 'N' ! no transpose

!      SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
!     $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
!     $                   WORK, IWORK, INFO )

            call DGESVX(fact, trans, nvar, nrhs, a, nvar, af, nvar, ipiv, &
                        equed, r, c, b, nvar, x, nvar, rcond, ferr, berr, &
                        work, iwork, ierr)
               
            if (ierr > 0 .and. ierr <= nvar) then ! singular
               write(*,3) 'singular matrix for DGESVX', k, ierr
               stop 'factor_with_DGESVX'
            end if
            if (ierr == nvar+1) then ! means bad rcond, but may not be fatal
               write(*,2) 'DGESVX reports bad matrix conditioning: k, rcond', k, rcond
               ierr = 0
            end if
            
            do i=1,nvar
               do j=1,nvar
                  mtx(i,j) = a(i,j)
                  mtxF(i,j) = af(i,j)
               end do
               row_scale_factors(i) = r(i)
               col_scale_factors(i) = c(i)
               ipivot(i) = ipiv(i)
            end do
            
            if (s% report_min_rcond_from_DGESXV .and. rcond < min_rcond_from_DGESVX) then
               !$OMP critical (bcyclic_dense_factor_crit)
               min_rcond_from_DGESVX = rcond
               k_min_rcond_from_DGESVX = k
               rpgfac = work(1)
               !$OMP end critical (bcyclic_dense_factor_crit)
            end if

         end subroutine factor_with_DGESVX
         
      end subroutine dense_factor
   

      subroutine bcyclic_solve ( &
            s, nvar, nz, lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, ipivot1, &
            B1, soln1, row_scale_factors1, col_scale_factors1, equed1, &
            iter, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar, nz, iter
         real(dp), pointer, dimension(:) :: &
            lblk1, dblk1, ublk1, lblkF1, dblkF1, ublkF1, &
            B1, soln1, row_scale_factors1, col_scale_factors1
         integer, pointer :: ipivot1(:)
         character (len=nz) :: equed1
         integer, intent(out) :: ierr

         integer, pointer :: iptr(:,:), nslevel(:), ipivot(:)
         integer :: ncycle, nstemp, maxlevels, nlevel, nvar2, i
         real(dp), pointer, dimension(:,:) :: dmat, dmatF
         real(dp), pointer, dimension(:) :: row_scale_factors, col_scale_factors
         character (len=1) :: equed
         logical :: okay

         include 'formats'


         if (dbg) write(*,*) 'start bcyclic_solve'
         
         ! copy B to soln
         do i=1,nvar*nz
            soln1(i) = B1(i)
         end do

         ierr = 0

         nvar2 = nvar*nvar
         ncycle = 1
         maxlevels = 0
         do while (ncycle < nz)
            ncycle = 2*ncycle
            maxlevels = maxlevels+1
         end do
         maxlevels = max(1, maxlevels)

         allocate (nslevel(maxlevels), stat=ierr)
         if (ierr /= 0) return

         ncycle = 1
         nstemp = nz
         nlevel = 1

         if (dbg) write(*,*) 'start forward_cycle'

         forward_cycle: do

            nslevel(nlevel) = nstemp
            if (dbg) write(*,2) 'call cycle_rhs', nstemp
            call cycle_rhs( &
               s, nz, nstemp, nvar, ncycle, nlevel, &
               dblk1, dblkF1, soln1, ipivot1, &
               row_scale_factors1, col_scale_factors1, equed1, ierr)
            if (ierr /= 0) then
               call dealloc
               return
            end if

            if (nstemp == 1) exit

            nstemp = (nstemp+1)/2
            nlevel = nlevel+1
            ncycle = 2*ncycle

            if (nlevel > maxlevels) exit

         end do forward_cycle

         if (dbg) write(*,*) 'done forward_cycle'

         dmat(1:nvar,1:nvar) => dblk1(1:nvar2)
         dmatF(1:nvar,1:nvar) => dblkF1(1:nvar2)
         ipivot(1:nvar) => ipivot1(1:nvar)
         row_scale_factors(1:nvar) => row_scale_factors1(1:nvar)
         col_scale_factors(1:nvar) => col_scale_factors1(1:nvar)
         equed(1:1) = equed1(1:1)
         call dense_solve1(s, 1, nvar, soln1, dmat, dmatF, ipivot, .false., &
            row_scale_factors, col_scale_factors, equed, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in my_getrs1'
            call dealloc
            return
         end if

         ! back solve for even x's
         back_cycle: do while (ncycle > 1)
            ncycle = ncycle/2
            nlevel = nlevel-1
            if (nlevel < 1) then
               ierr = -1
               exit
            end if
            nstemp = nslevel(nlevel)
            call cycle_solve( &
               s, nvar, nz, ncycle, nstemp, nlevel, &
               lblk1, ublk1, lblkF1, ublkF1, soln1)
         end do back_cycle

         call dealloc

         if (dbg) write(*,*) 'done bcyclic_solve'


         contains


         subroutine dealloc
            deallocate (nslevel)
         end subroutine dealloc


      end subroutine bcyclic_solve


      subroutine clear_storage(s)
         type (star_info), pointer :: s
         integer :: nlevel
         nlevel = size(s% bcyclic_odd_storage)
         do while (nlevel > 0)
            if (s% bcyclic_odd_storage(nlevel)% ul_size > 0) then
               deallocate(s% bcyclic_odd_storage(nlevel)% umat1)
               deallocate(s% bcyclic_odd_storage(nlevel)% lmat1)
            end if
            nlevel = nlevel-1
         end do
         deallocate(s% bcyclic_odd_storage)
         nullify(s% bcyclic_odd_storage)
      end subroutine clear_storage


      subroutine dense_solve(s, k, nvar, mtx, mtxF, ipivot, X_mtx, &
            row_scale_factors, col_scale_factors, equed, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         real(dp), pointer, dimension(:,:) :: mtx, mtxF, X_mtx
         integer, pointer :: ipivot(:)
         real(dp), pointer, dimension(:) :: row_scale_factors, col_scale_factors
         character (len=1) :: equed
         integer, intent(out) :: ierr
         integer :: i
         real(dp), pointer :: X(:)
         ierr = 0
         
         if (s% use_DGESVX_in_bcyclic) then
            call solve_with_DGESVX
            return
         end if

         do i=1,nvar
            X(1:nvar) => X_mtx(1:nvar,i)
            call dense_solve1(s, k, nvar, X, mtx, mtxF, ipivot, .false., &
               row_scale_factors, col_scale_factors, equed, ierr)
            if (ierr /= 0) return
         end do
         
         contains
         
         subroutine solve_with_DGESVX
            character (len=1) :: fact, trans
            real(dp) :: rcond
            real(dp) :: a(nvar,nvar), af(nvar,nvar), b(nvar,nvar), x(nvar,nvar), &
               r(nvar), c(nvar), ferr(nvar), berr(nvar), work(4*nvar)
            integer :: ipiv(nvar), iwork(nvar)
            integer :: i, j, nrhs
            include 'formats'

            nrhs = nvar

            do i=1,nvar
               do j=1,nvar
                  a(i,j) = mtx(i,j)
                  af(i,j) = mtxF(i,j)
                  b(i,j) = X_mtx(i,j)
                  x(i,j) = 0d0
               end do
               r(i) = row_scale_factors(i)
               c(i) = col_scale_factors(i)
               ipiv(i) = ipivot(i)
            end do
            
            fact = 'F' ! factored
            trans = 'N' ! no transpose

!      SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
!     $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
!     $                   WORK, IWORK, INFO )

            call DGESVX(fact, trans, nvar, nrhs, a, nvar, af, nvar, ipiv, &
                        equed, r, c, b, nvar, x, nvar, rcond, ferr, berr, &
                        work, iwork, ierr)
            if (ierr /= 0) then
               write(*,2) 'solve_with_DGESVX failed', k
            end if
            
            do i=1,nvar
               do j=1,nvar
                  X_mtx(i,j) = x(i,j)
               end do
            end do

         end subroutine solve_with_DGESVX

      end subroutine dense_solve


      subroutine dense_solve1(s, k, nvar, X_vec, mtx, mtxF, ipivot, dbg, &
            row_scale_factors, col_scale_factors, equed, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         real(dp), pointer :: X_vec(:), mtx(:,:), mtxF(:,:)
         integer, pointer :: ipivot(:)
         real(dp), pointer, dimension(:) :: row_scale_factors, col_scale_factors
         character (len=1) :: equed
         logical, intent(in) :: dbg
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         
         if (s% use_DGESVX_in_bcyclic) then
            call solve1_with_DGESVX
            return
         end if
         
         if (nvar == 4) then
            call my_getrs1_n4(mtxF, ipivot, X_vec, ierr)
         else if (nvar == 5) then
            call my_getrs1_n5(mtxF, ipivot, X_vec, ierr)
         else
            call my_getrs1(nvar, mtxF, nvar, ipivot, X_vec, nvar, ierr)
         end if
         
         contains
         
         subroutine solve1_with_DGESVX
            character (len=1) :: fact, trans
            real(dp) :: rcond
            integer, parameter :: nrhs = 1
            real(dp) :: a(nvar,nvar), af(nvar,nvar), b(nvar,nrhs), x(nvar,nrhs), &
               r(nvar), c(nvar), ferr(nrhs), berr(nrhs), work(4*nvar)
            integer :: ipiv(nvar), iwork(nvar)
            integer :: i, j

            include 'formats'

            do i=1,nvar
               do j=1,nvar
                  a(i,j) = mtx(i,j)
                  af(i,j) = mtxF(i,j)
               end do
               b(i,1) = X_vec(i)
               x(i,1) = 0d0
               r(i) = row_scale_factors(i)
               c(i) = col_scale_factors(i)
               ipiv(i) = ipivot(i)
            end do
            
            fact = 'F' ! factored
            trans = 'N' ! no transpose

!      SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
!     $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
!     $                   WORK, IWORK, INFO )

            call DGESVX(fact, trans, nvar, nrhs, a, nvar, af, nvar, ipiv, &
                        equed, r, c, b, nvar, x, nvar, rcond, ferr, berr, &
                        work, iwork, ierr)
            
            do i=1,nvar
               X_vec(i) = x(i,1)
            end do

         end subroutine solve1_with_DGESVX

      end subroutine dense_solve1



      subroutine bcyclic_deallocate (s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine bcyclic_deallocate


      include 'mtx_solve_routines.inc'


      subroutine my_gemm3_0_p1(m,n,k,a,b,c) ! c := a*b
         integer, intent(in) :: k,m,n
         real(dp), dimension(:,:) :: a, b, c
         call dgemm3('N', 'N', m, n, k, 1d0, a, n, b, n, 0d0, c, n)   
      end subroutine my_gemm3_0_p1

      subroutine my_gemm3_plus_mm(m,n,k,a,b,d,e,c) ! c := c + a*b + d*e
         integer, intent(in) :: k,m,n
         real(dp), dimension(:,:) :: a, b, c, d, e
         call my_gemm3_p1(m,n,k,a,b,c)   ! c = c + a*b
         call my_gemm3_p1(m,n,k,d,e,c)   ! c = c + d*e
      end subroutine my_gemm3_plus_mm

      subroutine my_gemm3_p1(m,n,k,a,b,c) ! c := c + a*b
         integer, intent(in) :: k,m,n
         real(dp), dimension(:,:) :: a, b, c 
         call dgemm3('N', 'N', m, n, k, 1d0, a, n, b, n, 1d0, c, n)
      end subroutine my_gemm3_p1
   
      subroutine my_gemm3(m,n,k,a,b,c) ! c := c - a*b
         integer, intent(in) :: k,m,n
         real(dp), dimension(:,:) :: a, b, c
         call dgemm3('N', 'N', m, n, k, -1d0, a, n, b, n, 1d0, c, n)
      end subroutine my_gemm3
   
      logical function lsame(ca,cb)
         character, intent(in) :: ca, cb
         integer :: inta, intb, zcode
         zcode = ichar('Z')
         inta = ichar(ca)
         intb = ichar(cb)
         if (zcode == 90 .or. zcode == 122) then
            ! ASCII is assumed - ZCODE is the ASCII code of either lower or
            ! upper case 'Z'.
            if (inta >= 97 .and. inta <= 122) inta = inta - 32
            if (intb >= 97 .and. intb <= 122) intb = intb - 32
         else 
            stop 'my lsame assumes ASCII'
         end if
         lsame = (inta == intb)
      end function lsame
   
      subroutine xerbla(srname, info)
         character (len=*), intent(in) :: srname
         integer, intent(in) :: info
   9999   format( ' ** On entry to ', A, ' parameter number ', I2, ' had ', &
            'an illegal value' )
         write( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info
      end subroutine xerbla
      
   ! Level 3 BLAS tuned for single processors with caches
   ! http://www.netlib.org/blas/gemm_based/  

      SUBROUTINE dgemm3 ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
   !     .. Scalar Arguments ..
         CHARACTER        TRANSA, TRANSB
         INTEGER            M, N, K, LDA, LDB, LDC
         DOUBLE PRECISION   ALPHA, BETA
   !     .. Array Arguments ..
         DOUBLE PRECISION   A( :, : ), B( :, : ), C( :, : )
         !DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
   !     ..
   !
   !  Purpose
   !  =======
   !
   !  DGEMM  performs one of the matrix-matrix operations
   !
   !     C := alpha*op( A )*op( B ) + beta*C,
   !
   !  where  op( X ) is one of
   !
   !     op( X ) = X   or   op( X ) = X',
   !
   !  alpha and beta are scalars, and A, B and C are matrices, with op( A )
   !  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
   !
   !  Parameters
   !  ==========
   !
   !  TRANSA - CHARACTER*1.
   !           On entry, TRANSA specifies the form of op( A ) to be used in
   !           the matrix multiplication as follows:
   !
   !              TRANSA = 'N' or 'n',  op( A ) = A.
   !
   !              TRANSA = 'T' or 't',  op( A ) = A'.
   !
   !              TRANSA = 'C' or 'c',  op( A ) = A'.
   !
   !           Unchanged on exit.
   !
   !  TRANSB - CHARACTER*1.
   !           On entry, TRANSB specifies the form of op( B ) to be used in
   !           the matrix multiplication as follows:
   !
   !              TRANSB = 'N' or 'n',  op( B ) = B.
   !
   !              TRANSB = 'T' or 't',  op( B ) = B'.
   !
   !              TRANSB = 'C' or 'c',  op( B ) = B'.
   !
   !           Unchanged on exit.
   !
   !  M      - INTEGER.
   !           On entry,  M  specifies  the number  of rows  of the  matrix
   !           op( A )  and of the  matrix  C.  M  must  be at least  zero.
   !           Unchanged on exit.
   !
   !  N      - INTEGER.
   !           On entry,  N  specifies the number  of columns of the matrix
   !           op( B ) and the number of columns of the matrix C. N must be
   !           at least zero.
   !           Unchanged on exit.
   !
   !  K      - INTEGER.
   !           On entry,  K  specifies  the number of columns of the matrix
   !           op( A ) and the number of rows of the matrix op( B ). K must
   !           be at least  zero.
   !           Unchanged on exit.
   !
   !  ALPHA  - DOUBLE PRECISION.
   !           On entry, ALPHA specifies the scalar alpha.
   !           Unchanged on exit.
   !
   !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
   !           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
   !           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
   !           part of the array  A  must contain the matrix  A,  otherwise
   !           the leading  k by m  part of the array  A  must contain  the
   !           matrix A.
   !           Unchanged on exit.
   !
   !  LDA    - INTEGER.
   !           On entry, LDA specifies the first dimension of A as declared
   !           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
   !           LDA must be at least  max( 1, m ), otherwise  LDA must be at
   !           least  max( 1, k ).
   !           Unchanged on exit.
   !
   !  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
   !           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
   !           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
   !           part of the array  B  must contain the matrix  B,  otherwise
   !           the leading  n by k  part of the array  B  must contain  the
   !           matrix B.
   !           Unchanged on exit.
   !
   !  LDB    - INTEGER.
   !           On entry, LDB specifies the first dimension of B as declared
   !           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
   !           LDB must be at least  max( 1, k ), otherwise  LDB must be at
   !           least  max( 1, n ).
   !           Unchanged on exit.
   !
   !  BETA   - DOUBLE PRECISION.
   !           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
   !           supplied as zero then C need not be set on input.
   !           Unchanged on exit.
   !
   !  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
   !           Before entry, the leading  m by n  part of the array  C must
   !           contain the matrix  C,  except when  beta  is zero, in which
   !           case C need not be set on entry.
   !           On exit, the array  C  is overwritten by the  m by n  matrix
   !           ( alpha*op( A )*op( B ) + beta*C ).
   !
   !  LDC    - INTEGER.
   !           On entry, LDC specifies the first dimension of C as declared
   !           in  the  calling  (sub)  program.   LDC  must  be  at  least
   !           max( 1, m ).
   !           Unchanged on exit.
   !
   !
   !  Level 3 Blas routine.
   !
   !  -- Written on 8-February-1989.
   !     Jack Dongarra, Argonne National Laboratory.
   !     Iain Duff, AERE Harwell.
   !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
   !     Sven Hammarling, Numerical Algorithms Group Ltd.
   !
   !  -- Modified in October-1997.
   !     Superscalar GEMM-Based Level 3 BLAS (Version 0.1).
   !     Per Ling, Department of Computing Science,
   !     Umea University, Sweden.
   !
   !
   !     .. Local Scalars ..
         INTEGER            I, II, ISEC, UISEC, J, JJ, JSEC, UJSEC, &
                           L, LL, LSEC, ULSEC, INFO, NROWA, NROWB
         LOGICAL            NOTA, NOTB
         DOUBLE PRECISION   DELTA
         DOUBLE PRECISION   F11, F12, F21, F22, F31, F32, F41, F42
         DOUBLE PRECISION   F13, F14, F23, F24, F33, F34, F43, F44
   !     .. Intrinsic Functions ..
         INTRINSIC          MAX, MIN, MOD
   !     .. External Functions ..
         LOGICAL            LSAME
         EXTERNAL           LSAME
   !     .. External Subroutines ..
         EXTERNAL           XERBLA
   !     .. Parameters ..
         DOUBLE PRECISION   ZERO, ONE
         PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
   !     .. User specified parameters for DGEMM ..
         INTEGER            MB, NB, NBT, KB
         PARAMETER        ( MB = 32, NB = 1024, NBT = 96, KB = 32 )
         DOUBLE PRECISION   T1( KB, MB ), T2( KB, NBT )
   !     ..
   !     .. Executable Statements ..
   !
   !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
   !     transposed and set NROWA and NROWB as the number of rows of A  and
   !     the number of rows of B respectively.
   !
         NOTA = LSAME( TRANSA, 'N' )
         NOTB = LSAME( TRANSB, 'N' )
         IF ( NOTA ) THEN
            NROWA = M
         ELSE
            NROWA = K
         END IF
         IF ( NOTB ) THEN
            NROWB = K
         ELSE
            NROWB = N
         END IF
   !
   !     Test the input parameters.
   !
         INFO = 0
         IF( ( .NOT.NOTA ).AND.( .NOT. LSAME( TRANSA,'C' ) ) &
                                  .AND.( .NOT.LSAME( TRANSA, 'T' ) ) )THEN
            INFO = 1
         ELSE IF( ( .NOT.NOTB ).AND.( .NOT.LSAME( TRANSB, 'C' ) ) &
                                  .AND.( .NOT.LSAME( TRANSB, 'T' ) ) )THEN
            INFO = 2
         ELSE IF( M.LT.0 )THEN
            INFO = 3
         ELSE IF( N.LT.0 )THEN
            INFO = 4
         ELSE IF( K.LT.0 )THEN
            INFO = 5
         ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
            INFO = 8
         ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
            INFO = 10
         ELSE IF( LDC.LT.MAX( 1, M ) )THEN
            INFO = 13
         END IF
         IF( INFO.NE.0 )THEN
            CALL XERBLA( 'DGEMM ', INFO )
            RETURN
         END IF
   !
   !     Quick return if possible.
   !
         IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) then
            !write(*,*) 'quick return since ALPHA == 0 and BETA == 1'
            RETURN
         end if
   !
   !     And when alpha.eq.zero.
   !
         IF( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) )THEN
            IF( BETA.EQ.ZERO )THEN
               UISEC = M-MOD( M, 4 )
               DO 30, J = 1, N
                  DO 10, I = 1, UISEC, 4
                     C( I, J ) = ZERO
                     C( I+1, J ) = ZERO
                     C( I+2, J ) = ZERO
                     C( I+3, J ) = ZERO
      10          CONTINUE
                  DO 20, I = UISEC+1, M
                     C( I, J ) = ZERO
      20          CONTINUE
      30       CONTINUE
            ELSE
               UISEC = M-MOD( M, 4 )
               DO 60, J = 1, N
                  DO 40, I = 1, UISEC, 4
                     C( I, J ) = BETA*C( I, J )
                     C( I+1, J ) = BETA*C( I+1, J )
                     C( I+2, J ) = BETA*C( I+2, J )
                     C( I+3, J ) = BETA*C( I+3, J )
      40          CONTINUE
                  DO 50, I = UISEC+1, M
                     C( I, J ) = BETA*C( I, J )
      50          CONTINUE
      60       CONTINUE
            END IF
            RETURN
         END IF
   !
   !     Start the operations.
   
         if (nota .and. notb .and. alpha == 1d0 &
               .and. (beta == 0d0 .or. beta == 1d0)) then
            ! special case
            ! C := A*B  or  C := A*B + C



            DO 1250 JJ = 1, N, NB
               JSEC = MIN( NB, N-JJ+1 )
               UJSEC = JSEC-MOD( JSEC, 4 )
               DO 1240 LL = 1, K, KB
                  LSEC = MIN( KB, K-LL+1 )
                  ULSEC = LSEC-MOD( LSEC, 2 )
   !
   !              Determine if the block of C should be updated with
   !              beta or not.
   !
                  DELTA = ONE
                  IF( LL.EQ.1 ) DELTA = BETA
   !
                  DO 1230 II = 1, M, MB
                     ISEC = MIN( MB, M-II+1 )
   !
                     UISEC = ISEC-MOD( ISEC, 2 )

                        DO 1080, L = LL, LL+ULSEC-1, 2
                           !$omp simd
                           DO 1070, I = II, II+UISEC-1, 2
                              T1( L-LL+1, I-II+1 ) = A( I, L )
                              T1( L-LL+2, I-II+1 ) = A( I, L+1 )
                              T1( L-LL+1, I-II+2 ) = A( I+1, L )
                              T1( L-LL+2, I-II+2 ) = A( I+1, L+1 )
    1070                   CONTINUE
                           IF( UISEC.LT.ISEC )THEN
                              T1( L-LL+1, ISEC ) = A( II+ISEC-1, L )
                              T1( L-LL+2, ISEC ) = A( II+ISEC-1, L+1 )
                           END IF
    1080                CONTINUE
                        IF( ULSEC.LT.LSEC )THEN
                           !$omp simd
                           DO 1090, I = II, II+ISEC-1
                              T1( LSEC, I-II+1 ) = A( I, LL+LSEC-1 )
    1090                   CONTINUE
                        END IF
   !
   !                 C := T1'*B + beta*C, update a rectangular block
   !                 of C using 4 by 4 unrolling.
   !
                     UISEC = ISEC-MOD( ISEC, 4 )
                     DO 1170 J = JJ, JJ+UJSEC-1, 4
                        DO 1140 I = II, II+UISEC-1, 4
                           if (delta == 1d0) then
                              F11 = C( I,J )
                              F21 = C( I+1,J )
                              F12 = C( I,J+1 )
                              F22 = C( I+1,J+1 )
                              F13 = C( I,J+2 )
                              F23 = C( I+1,J+2 )
                              F14 = C( I,J+3 )
                              F24 = C( I+1,J+3 )
                              F31 = C( I+2,J )
                              F41 = C( I+3,J )
                              F32 = C( I+2,J+1 )
                              F42 = C( I+3,J+1 )
                              F33 = C( I+2,J+2 )
                              F43 = C( I+3,J+2 )
                              F34 = C( I+2,J+3 )
                              F44 = C( I+3,J+3 )
                           else if (delta == 0d0) then
                              F11 = 0d0
                              F21 = 0d0
                              F12 = 0d0
                              F22 = 0d0
                              F13 = 0d0
                              F23 = 0d0
                              F14 = 0d0
                              F24 = 0d0
                              F31 = 0d0
                              F41 = 0d0
                              F32 = 0d0
                              F42 = 0d0
                              F33 = 0d0
                              F43 = 0d0
                              F34 = 0d0
                              F44 = 0d0
                           else
                              stop 'dgemm special case requires beta == 0 or 1'
                           end if
                           !$omp simd
                           DO 1130 L = LL, LL+LSEC-1
                              F11 = F11 + T1( L-LL+1, I-II+1 )*B( L, J )
                              F21 = F21 + T1( L-LL+1, I-II+2 )*B( L, J )
                              F12 = F12 + T1( L-LL+1, I-II+1 )*B( L, J+1 )
                              F22 = F22 + T1( L-LL+1, I-II+2 )*B( L, J+1 )
                              F13 = F13 + T1( L-LL+1, I-II+1 )*B( L, J+2 )
                              F23 = F23 + T1( L-LL+1, I-II+2 )*B( L, J+2 )
                              F14 = F14 + T1( L-LL+1, I-II+1 )*B( L, J+3 )
                              F24 = F24 + T1( L-LL+1, I-II+2 )*B( L, J+3 )
                              F31 = F31 + T1( L-LL+1, I-II+3 )*B( L, J )
                              F41 = F41 + T1( L-LL+1, I-II+4 )*B( L, J )
                              F32 = F32 + T1( L-LL+1, I-II+3 )*B( L, J+1 )
                              F42 = F42 + T1( L-LL+1, I-II+4 )*B( L, J+1 )
                              F33 = F33 + T1( L-LL+1, I-II+3 )*B( L, J+2 )
                              F43 = F43 + T1( L-LL+1, I-II+4 )*B( L, J+2 )
                              F34 = F34 + T1( L-LL+1, I-II+3 )*B( L, J+3 )
                              F44 = F44 + T1( L-LL+1, I-II+4 )*B( L, J+3 )
    1130                   CONTINUE
                           C( I,J ) = F11
                           C( I+1, J ) = F21
                           C( I, J+1 ) = F12
                           C( I+1, J+1 ) = F22
                           C( I, J+2 ) = F13
                           C( I+1, J+2 ) = F23
                           C( I, J+3 ) = F14
                           C( I+1, J+3 ) = F24
                           C( I+2, J ) = F31
                           C( I+3, J ) = F41
                           C( I+2, J+1 ) = F32
                           C( I+3, J+1 ) = F42
                           C( I+2, J+2 ) = F33
                           C( I+3, J+2 ) = F43
                           C( I+2, J+3 ) = F34
                           C( I+3, J+3 ) = F44
    1140                CONTINUE
                        IF( UISEC.LT.ISEC )THEN
                           DO 1160 I = II+UISEC, II+ISEC-1
                              F11 = C( I, J )
                              F12 = C( I, J+1 )
                              F13 = C( I, J+2 )
                              F14 = C( I, J+3 )
                              !$omp simd
                              DO 1150 L = LL, LL+LSEC-1
                                 F11 = F11 + T1( L-LL+1, I-II+1 )*B( L, J )
                                 F12 = F12 + T1( L-LL+1, I-II+1 )* &
                                                               B( L, J+1 )
                                 F13 = F13 + T1( L-LL+1, I-II+1 )* &
                                                               B( L, J+2 )
                                 F14 = F14 + T1( L-LL+1, I-II+1 )* &
                                                               B( L, J+3 )
    1150                      CONTINUE
                              C( I, J ) = F11
                              C( I, J+1 ) = F12
                              C( I, J+2 ) = F13
                              C( I, J+3 ) = F14
    1160                   CONTINUE
                        END IF
    1170             CONTINUE
                     IF( UJSEC.LT.JSEC )THEN
                        DO 1220 J = JJ+UJSEC, JJ+JSEC-1
                           DO 1190 I = II, II+UISEC-1, 4
                              F11 = C( I,J )
                              F21 = C( I+1, J )
                              F31 = C( I+2, J )
                              F41 = C( I+3, J )
                              !$omp simd
                              DO 1180 L = LL, LL+LSEC-1
                                 F11 = F11 + T1( L-LL+1, I-II+1 )*B( L, J )
                                 F21 = F21 + T1( L-LL+1, I-II+2 )*B( L, J )
                                 F31 = F31 + T1( L-LL+1, I-II+3 )*B( L, J )
                                 F41 = F41 + T1( L-LL+1, I-II+4 )*B( L, J )
    1180                      CONTINUE
                              C( I,J ) = F11
                              C( I+1, J ) = F21
                              C( I+2, J ) = F31
                              C( I+3, J ) = F41
    1190                   CONTINUE
                           DO 1210 I = II+UISEC, II+ISEC-1
                              F11 = C( I, J )
                              !$omp simd
                              DO 1200 L = LL, LL+LSEC-1
                                 F11 = F11 + T1( L-LL+1, I-II+1 )*B( L, J )
    1200                      CONTINUE
                              C( I, J ) = F11
    1210                   CONTINUE
    1220                CONTINUE
                     END IF
    1230          CONTINUE
    1240       CONTINUE
    1250    CONTINUE
          
          return ! end of special case

         end if   
   !
         IF (NOTB) THEN
   !
   !        Form  C := alpha*A*B + beta*C or C := alpha*A'*B + beta*C.
   !
            DO 250 JJ = 1, N, NB
               JSEC = MIN( NB, N-JJ+1 )
               UJSEC = JSEC-MOD( JSEC, 4 )
               DO 240 LL = 1, K, KB
                  LSEC = MIN( KB, K-LL+1 )
                  ULSEC = LSEC-MOD( LSEC, 2 )
   !
   !              Determine if the block of C should be updated with
   !              beta or not.
   !
                  DELTA = ONE
                  IF( LL.EQ.1 ) DELTA = BETA
   !
                  DO 230 II = 1, M, MB
                     ISEC = MIN( MB, M-II+1 )
   !
   !                 T1 := alpha*A' or T1 := alpha*A, copy the transpose
   !                 or the non-transpose of a rectangular block of
   !                 alpha*A to T1.
   !
                     UISEC = ISEC-MOD( ISEC, 2 )
                     IF( NOTA )THEN
                        DO 80, L = LL, LL+ULSEC-1, 2
                           DO 70, I = II, II+UISEC-1, 2
                              T1( L-LL+1, I-II+1 ) = ALPHA*A( I, L )
                              T1( L-LL+2, I-II+1 ) = ALPHA*A( I, L+1 )
                              T1( L-LL+1, I-II+2 ) = ALPHA*A( I+1, L )
                              T1( L-LL+2, I-II+2 ) = ALPHA*A( I+1, L+1 )
      70                   CONTINUE
                           IF( UISEC.LT.ISEC )THEN
                              T1( L-LL+1, ISEC ) = ALPHA*A( II+ISEC-1, L )
                              T1( L-LL+2, ISEC ) = ALPHA*A( II+ISEC-1, L+1 )
                           END IF
      80                CONTINUE
                        IF( ULSEC.LT.LSEC )THEN
                           DO 90, I = II, II+ISEC-1
                              T1( LSEC, I-II+1 ) = ALPHA*A( I, LL+LSEC-1 )
      90                   CONTINUE
                        END IF
                     ELSE
                        DO 110, I = II, II+UISEC-1, 2
                           DO 100, L = LL, LL+ULSEC-1, 2
                              T1( L-LL+1, I-II+1 ) = ALPHA*A( L, I )
                              T1( L-LL+1, I-II+2 ) = ALPHA*A( L, I+1 )
                              T1( L-LL+2, I-II+1 ) = ALPHA*A( L+1, I )
                              T1( L-LL+2, I-II+2 ) = ALPHA*A( L+1, I+1 )
     100                   CONTINUE
                           IF( ULSEC.LT.LSEC )THEN
                              T1( LSEC, I-II+1 ) = ALPHA*A( LL+LSEC-1, I )
                              T1( LSEC, I-II+2 ) = ALPHA*A( LL+LSEC-1, I+1 )
                           END IF
     110                CONTINUE
                        IF( UISEC.LT.ISEC )THEN
                           DO 120, L = LL, LL+LSEC-1
                              T1( L-LL+1, ISEC ) = ALPHA*A( L, II+ISEC-1 )
     120                   CONTINUE
                        END IF
                     END IF
   !
   !                 C := T1'*B + beta*C, update a rectangular block
   !                 of C using 4 by 4 unrolling.
   !
                     UISEC = ISEC-MOD( ISEC, 4 )
                     DO 170 J = JJ, JJ+UJSEC-1, 4
                        DO 140 I = II, II+UISEC-1, 4
                           F11 = DELTA*C( I,J )
                           F21 = DELTA*C( I+1,J )
                           F12 = DELTA*C( I,J+1 )
                           F22 = DELTA*C( I+1,J+1 )
                           F13 = DELTA*C( I,J+2 )
                           F23 = DELTA*C( I+1,J+2 )
                           F14 = DELTA*C( I,J+3 )
                           F24 = DELTA*C( I+1,J+3 )
                           F31 = DELTA*C( I+2,J )
                           F41 = DELTA*C( I+3,J )
                           F32 = DELTA*C( I+2,J+1 )
                           F42 = DELTA*C( I+3,J+1 )
                           F33 = DELTA*C( I+2,J+2 )
                           F43 = DELTA*C( I+3,J+2 )
                           F34 = DELTA*C( I+2,J+3 )
                           F44 = DELTA*C( I+3,J+3 )
                           DO 130 L = LL, LL+LSEC-1
                              F11 = F11 + T1( L-LL+1, I-II+1 )*B( L, J )
                              F21 = F21 + T1( L-LL+1, I-II+2 )*B( L, J )
                              F12 = F12 + T1( L-LL+1, I-II+1 )*B( L, J+1 )
                              F22 = F22 + T1( L-LL+1, I-II+2 )*B( L, J+1 )
                              F13 = F13 + T1( L-LL+1, I-II+1 )*B( L, J+2 )
                              F23 = F23 + T1( L-LL+1, I-II+2 )*B( L, J+2 )
                              F14 = F14 + T1( L-LL+1, I-II+1 )*B( L, J+3 )
                              F24 = F24 + T1( L-LL+1, I-II+2 )*B( L, J+3 )
                              F31 = F31 + T1( L-LL+1, I-II+3 )*B( L, J )
                              F41 = F41 + T1( L-LL+1, I-II+4 )*B( L, J )
                              F32 = F32 + T1( L-LL+1, I-II+3 )*B( L, J+1 )
                              F42 = F42 + T1( L-LL+1, I-II+4 )*B( L, J+1 )
                              F33 = F33 + T1( L-LL+1, I-II+3 )*B( L, J+2 )
                              F43 = F43 + T1( L-LL+1, I-II+4 )*B( L, J+2 )
                              F34 = F34 + T1( L-LL+1, I-II+3 )*B( L, J+3 )
                              F44 = F44 + T1( L-LL+1, I-II+4 )*B( L, J+3 )
     130                   CONTINUE
                           C( I,J ) = F11
                           C( I+1, J ) = F21
                           C( I, J+1 ) = F12
                           C( I+1, J+1 ) = F22
                           C( I, J+2 ) = F13
                           C( I+1, J+2 ) = F23
                           C( I, J+3 ) = F14
                           C( I+1, J+3 ) = F24
                           C( I+2, J ) = F31
                           C( I+3, J ) = F41
                           C( I+2, J+1 ) = F32
                           C( I+3, J+1 ) = F42
                           C( I+2, J+2 ) = F33
                           C( I+3, J+2 ) = F43
                           C( I+2, J+3 ) = F34
                           C( I+3, J+3 ) = F44
     140                CONTINUE
                        IF( UISEC.LT.ISEC )THEN
                           DO 160 I = II+UISEC, II+ISEC-1
                              F11 = DELTA*C( I, J )
                              F12 = DELTA*C( I, J+1 )
                              F13 = DELTA*C( I, J+2 )
                              F14 = DELTA*C( I, J+3 )
                              DO 150 L = LL, LL+LSEC-1
                                 F11 = F11 + T1( L-LL+1, I-II+1 )*B( L, J )
                                 F12 = F12 + T1( L-LL+1, I-II+1 )* &
                                                               B( L, J+1 )
                                 F13 = F13 + T1( L-LL+1, I-II+1 )* &
                                                               B( L, J+2 )
                                 F14 = F14 + T1( L-LL+1, I-II+1 )* &
                                                               B( L, J+3 )
     150                      CONTINUE
                              C( I, J ) = F11
                              C( I, J+1 ) = F12
                              C( I, J+2 ) = F13
                              C( I, J+3 ) = F14
     160                   CONTINUE
                        END IF
     170             CONTINUE
                     IF( UJSEC.LT.JSEC )THEN
                        DO 220 J = JJ+UJSEC, JJ+JSEC-1
                           DO 190 I = II, II+UISEC-1, 4
                              F11 = DELTA*C( I,J )
                              F21 = DELTA*C( I+1, J )
                              F31 = DELTA*C( I+2, J )
                              F41 = DELTA*C( I+3, J )
                              DO 180 L = LL, LL+LSEC-1
                                 F11 = F11 + T1( L-LL+1, I-II+1 )*B( L, J )
                                 F21 = F21 + T1( L-LL+1, I-II+2 )*B( L, J )
                                 F31 = F31 + T1( L-LL+1, I-II+3 )*B( L, J )
                                 F41 = F41 + T1( L-LL+1, I-II+4 )*B( L, J )
     180                      CONTINUE
                              C( I,J ) = F11
                              C( I+1, J ) = F21
                              C( I+2, J ) = F31
                              C( I+3, J ) = F41
     190                   CONTINUE
                           DO 210 I = II+UISEC, II+ISEC-1
                              F11 = DELTA*C( I, J )
                              DO 200 L = LL, LL+LSEC-1
                                 F11 = F11 + T1( L-LL+1, I-II+1 )*B( L, J )
     200                      CONTINUE
                              C( I, J ) = F11
     210                   CONTINUE
     220                CONTINUE
                     END IF
     230          CONTINUE
     240       CONTINUE
     250    CONTINUE
         ELSE
   !
   !        Form  C := alpha*A*B' + beta*C or C := alpha*A'*B' + beta*C.
   !
            DO 470 JJ = 1, N, NBT
               JSEC = MIN( NBT, N-JJ+1 )
               DO 460 LL = 1, K, KB
                  LSEC = MIN( KB, K-LL+1 )
   !
   !              Determine if the block of C should be updated with
   !              beta or not.
   !
                  DELTA = ONE
                  IF( LL.EQ.1 ) DELTA = BETA
   !
   !              T2 := alpha*B', copy the transpose of a rectangular
   !              block of alpha*A to T2.
   !
                  ULSEC = LSEC-MOD( LSEC, 2 )
                  UJSEC = JSEC-MOD( JSEC, 2 )
                  DO 270, L = LL, LL+ULSEC-1, 2
                     DO 260, J = JJ, JJ+UJSEC-1, 2
                        T2( L-LL+1, J-JJ+1 ) = ALPHA*B( J, L )
                        T2( L-LL+2, J-JJ+1 ) = ALPHA*B( J, L+1 )
                        T2( L-LL+1, J-JJ+2 ) = ALPHA*B( J+1, L )
                        T2( L-LL+2, J-JJ+2 ) = ALPHA*B( J+1, L+1 )
     260             CONTINUE
                     IF( UJSEC.LT.JSEC )THEN
                        T2( L-LL+1, JSEC ) = ALPHA*B( JJ+JSEC-1, L )
                        T2( L-LL+2, JSEC ) = ALPHA*B( JJ+JSEC-1, L+1 )
                     END IF
     270          CONTINUE
                  IF( ULSEC.LT.LSEC )THEN
                     DO 280, J = JJ, JJ+JSEC-1
                        T2( LSEC, J-JJ+1 ) = ALPHA*B( J, LL+LSEC-1 )
     280             CONTINUE
                  END IF
   !
                  UJSEC = JSEC-MOD( JSEC, 4 )
                  DO 450 II = 1, M, MB
                     ISEC = MIN( MB, M-II+1 )
   !
   !                 T1 := alpha*A' or T1 := alpha*A, copy the transpose
   !                 or the non-transpose of a rectangular block of
   !                 alpha*A to T1.
   !
                     UISEC = ISEC-MOD( ISEC, 2 )
                     IF( NOTA )THEN
                        DO 300, L = LL, LL+ULSEC-1, 2
                           DO 290, I = II, II+UISEC-1, 2
                              T1( L-LL+1, I-II+1 ) = A( I, L )
                              T1( L-LL+2, I-II+1 ) = A( I, L+1 )
                              T1( L-LL+1, I-II+2 ) = A( I+1, L )
                              T1( L-LL+2, I-II+2 ) = A( I+1, L+1 )
     290                   CONTINUE
                           IF( UISEC.LT.ISEC )THEN
                              T1( L-LL+1, ISEC ) = A( II+ISEC-1, L )
                              T1( L-LL+2, ISEC ) = A( II+ISEC-1, L+1 )
                           END IF
     300                CONTINUE
                        IF( ULSEC.LT.LSEC )THEN
                           DO 310, I = II, II+ISEC-1
                              T1( LSEC, I-II+1 ) = A( I, LL+LSEC-1 )
     310                   CONTINUE
                        END IF
                     ELSE
                        DO 330, I = II, II+UISEC-1, 2
                           DO 320, L = LL, LL+ULSEC-1, 2
                              T1( L-LL+1, I-II+1 ) = A( L, I )
                              T1( L-LL+1, I-II+2 ) = A( L, I+1 )
                              T1( L-LL+2, I-II+1 ) = A( L+1, I )
                              T1( L-LL+2, I-II+2 ) = A( L+1, I+1 )
     320                   CONTINUE
                           IF( ULSEC.LT.LSEC )THEN
                              T1( LSEC, I-II+1 ) = A( LL+LSEC-1, I )
                              T1( LSEC, I-II+2 ) = A( LL+LSEC-1, I+1 )
                           END IF
     330                CONTINUE
                        IF( UISEC.LT.ISEC )THEN
                           DO 340, L = LL, LL+LSEC-1
                              T1( L-LL+1, ISEC ) = A( L, II+ISEC-1 )
     340                   CONTINUE
                        END IF
                     END IF
   !
   !                 C := T1'*B + beta*C, update a rectangular block
   !                 of C using 4 by 4 unrolling.
   !
                     UISEC = ISEC-MOD( ISEC, 4 )
                     DO 390 J = JJ, JJ+UJSEC-1, 4
                        DO 360 I = II, II+UISEC-1, 4
                           F11 = DELTA*C( I,J )
                           F21 = DELTA*C( I+1,J )
                           F12 = DELTA*C( I,J+1 )
                           F22 = DELTA*C( I+1,J+1 )
                           F13 = DELTA*C( I,J+2 )
                           F23 = DELTA*C( I+1,J+2 )
                           F14 = DELTA*C( I,J+3 )
                           F24 = DELTA*C( I+1,J+3 )
                           F31 = DELTA*C( I+2,J )
                           F41 = DELTA*C( I+3,J )
                           F32 = DELTA*C( I+2,J+1 )
                           F42 = DELTA*C( I+3,J+1 )
                           F33 = DELTA*C( I+2,J+2 )
                           F43 = DELTA*C( I+3,J+2 )
                           F34 = DELTA*C( I+2,J+3 )
                           F44 = DELTA*C( I+3,J+3 )
                           DO 350 L = LL, LL+LSEC-1
                              F11 = F11 + T1( L-LL+1, I-II+1 )* &
                                                      T2( L-LL+1, J-JJ+1 )
                              F21 = F21 + T1( L-LL+1, I-II+2 )* &
                                                      T2( L-LL+1, J-JJ+1 )
                              F12 = F12 + T1( L-LL+1, I-II+1 )* &
                                                      T2( L-LL+1, J-JJ+2 )
                              F22 = F22 + T1( L-LL+1, I-II+2 )* &
                                                      T2( L-LL+1, J-JJ+2 )
                              F13 = F13 + T1( L-LL+1, I-II+1 )* &
                                                      T2( L-LL+1, J-JJ+3 )
                              F23 = F23 + T1( L-LL+1, I-II+2 )* &
                                                      T2( L-LL+1, J-JJ+3 )
                              F14 = F14 + T1( L-LL+1, I-II+1 )* &
                                                      T2( L-LL+1, J-JJ+4 )
                              F24 = F24 + T1( L-LL+1, I-II+2 )* &
                                                      T2( L-LL+1, J-JJ+4 )
                              F31 = F31 + T1( L-LL+1, I-II+3 )* &
                                                      T2( L-LL+1, J-JJ+1 )
                              F41 = F41 + T1( L-LL+1, I-II+4 )* &
                                                      T2( L-LL+1, J-JJ+1 )
                              F32 = F32 + T1( L-LL+1, I-II+3 )* &
                                                      T2( L-LL+1, J-JJ+2 )
                              F42 = F42 + T1( L-LL+1, I-II+4 )* &
                                                      T2( L-LL+1, J-JJ+2 )
                              F33 = F33 + T1( L-LL+1, I-II+3 )* &
                                                      T2( L-LL+1, J-JJ+3 )
                              F43 = F43 + T1( L-LL+1, I-II+4 )* &
                                                      T2( L-LL+1, J-JJ+3 )
                              F34 = F34 + T1( L-LL+1, I-II+3 )* &
                                                      T2( L-LL+1, J-JJ+4 )
                              F44 = F44 + T1( L-LL+1, I-II+4 )* &
                                                      T2( L-LL+1, J-JJ+4 )
     350                   CONTINUE
                           C( I,J ) = F11
                           C( I+1, J ) = F21
                           C( I, J+1 ) = F12
                           C( I+1, J+1 ) = F22
                           C( I, J+2 ) = F13
                           C( I+1, J+2 ) = F23
                           C( I, J+3 ) = F14
                           C( I+1, J+3 ) = F24
                           C( I+2, J ) = F31
                           C( I+3, J ) = F41
                           C( I+2, J+1 ) = F32
                           C( I+3, J+1 ) = F42
                           C( I+2, J+2 ) = F33
                           C( I+3, J+2 ) = F43
                           C( I+2, J+3 ) = F34
                           C( I+3, J+3 ) = F44
     360                CONTINUE
                        IF( UISEC.LT.ISEC )THEN
                           DO 380 I = II+UISEC, II+ISEC-1
                              F11 = DELTA*C( I, J )
                              F12 = DELTA*C( I, J+1 )
                              F13 = DELTA*C( I, J+2 )
                              F14 = DELTA*C( I, J+3 )
                              DO 370 L = LL, LL+LSEC-1
                                 F11 = F11 + T1( L-LL+1, I-II+1 )* &
                                                      T2( L-LL+1, J-JJ+1 )
                                 F12 = F12 + T1( L-LL+1, I-II+1 )* &
                                                       T2( L-LL+1, J-JJ+2 )
                                 F13 = F13 + T1( L-LL+1, I-II+1 )* &
                                                      T2( L-LL+1, J-JJ+3 )
                                 F14 = F14 + T1( L-LL+1, I-II+1 )* &
                                                      T2( L-LL+1, J-JJ+4 )
     370                      CONTINUE
                              C( I, J ) = F11
                              C( I, J+1 ) = F12
                              C( I, J+2 ) = F13
                              C( I, J+3 ) = F14
     380                   CONTINUE
                        END IF
     390             CONTINUE
                     IF( UJSEC.LT.JSEC )THEN
                        DO 440 J = JJ+UJSEC, JJ+JSEC-1
                           DO 410 I = II, II+UISEC-1, 4
                              F11 = DELTA*C( I, J )
                              F21 = DELTA*C( I+1, J )
                              F31 = DELTA*C( I+2, J )
                              F41 = DELTA*C( I+3, J )
                              DO 400 L = LL, LL+LSEC-1
                                 F11 = F11 + T1( L-LL+1, I-II+1 )* &
                                                      T2( L-LL+1, J-JJ+1 )
                                F21 = F21 + T1( L-LL+1, I-II+2 )* &
                                                      T2( L-LL+1, J-JJ+1 )
                                 F31 = F31 + T1( L-LL+1, I-II+3 )* &
                                                      T2( L-LL+1, J-JJ+1 )
                                 F41 = F41 + T1( L-LL+1, I-II+4 )* &
                                                      T2( L-LL+1, J-JJ+1 )
     400                      CONTINUE
                              C( I,J ) = F11
                              C( I+1, J ) = F21
                              C( I+2, J ) = F31
                              C( I+3, J ) = F41
     410                   CONTINUE
                           DO 430 I = II+UISEC, II+ISEC-1
                              F11 = DELTA*C( I, J )
                              DO 420 L = LL, LL+LSEC-1
                                 F11 = F11 + T1( L-LL+1, I-II+1 )* &
                                                      T2( L-LL+1, J-JJ+1 )
     420                      CONTINUE
                              C( I, J ) = F11
     430                   CONTINUE
     440                CONTINUE
                     END IF
     450          CONTINUE
     460       CONTINUE
     470    CONTINUE
         END IF
   !
         RETURN
   !
   !     End of DGEMM.
   !
      END SUBROUTINE dgemm3

      end module star_bcyclic
