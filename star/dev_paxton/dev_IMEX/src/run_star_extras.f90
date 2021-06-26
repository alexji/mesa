! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      
      implicit none
      
      include 'test_suite_extras_def.inc'
      
      !real(dp) :: example_of_data_to_save_in_photos
      
      contains

      include 'test_suite_extras.inc'
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         include 'read_imex_controls.inc'
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         s% other_photo_read => extras_photo_read
         s% other_photo_write => extras_photo_write
      end subroutine extras_controls


      include 'imex.inc'


      subroutine extras_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         ierr = 0
         !read(iounit,iostat=ierr) example_of_data_to_save_in_photos
      end subroutine extras_photo_read


      subroutine extras_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         !write(iounit) example_of_data_to_save_in_photos
      end subroutine extras_photo_write
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 9
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         use imex_work
         use imex_output, only: prim, T
         use imex, only: get_imex_total_energies
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: i, k
         real(dp) :: &
            sum_EKIN_actual, sum_EKIN_for_ETOT, sum_EGAS, sum_ERAD, sum_ETOT
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call get_imex_total_energies(prim, T, &
            sum_EKIN_actual, sum_EKIN_for_ETOT, sum_EGAS, sum_ERAD, sum_ETOT)
         i = 1
         names(i) = 'dt_advec'; vals(i) = dt_advection; i=i+1
         names(i) = 'dt_grid'; vals(i) = dt_grid; i=i+1
         names(i) = 'dt_max'; vals(i) = dt_max_new; i=i+1
         names(i) = 'dt_rel_dE'; vals(i) = dt_rel_dE; i=i+1
         names(i) = 'dt_front'; vals(i) = dt_front_v; i=i+1
         names(i) = 'sum_ETOT'; vals(i) = sum_ETOT; i=i+1
         names(i) = 'sum_EKIN'; vals(i) = sum_EKIN_actual; i=i+1
         names(i) = 'sum_EGAS'; vals(i) = sum_EGAS; i=i+1
         names(i) = 'sum_ERAD'; vals(i) = sum_ERAD; i=i+1
         !names(i) = ''; vals(i) = ; i=i+1
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 16
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         use imex, only: set_imex_plot_data
         use imex_plot_data
         use imex_work, only: imex1, imex2, imex3, imex4, imex5, imex6
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, i
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call set_imex_plot_data()
         i = 1
         names(i) = 'plot_r'; i=i+1
         names(i) = 'plot_L'; i=i+1
         names(i) = 'plot_rho'; i=i+1
         names(i) = 'plot_T'; i=i+1
         names(i) = 'plot_v'; i=i+1
         names(i) = 'plot_csound'; i=i+1
         names(i) = 'plot_etot'; i=i+1
         names(i) = 'plot_ekin'; i=i+1
         names(i) = 'plot_eeos'; i=i+1
         names(i) = 'plot_Peos'; i=i+1
         names(i) = 'imex1'; i=i+1
         names(i) = 'imex2'; i=i+1
         names(i) = 'imex3'; i=i+1
         names(i) = 'imex4'; i=i+1
         names(i) = 'imex5'; i=i+1
         names(i) = 'imex6'; i=i+1
         
         do k=1,nz
            i = 1
            vals(k,i) = plot_r(k); i=i+1
            vals(k,i) = plot_L(k); i=i+1
            vals(k,i) = plot_rho(k); i=i+1
            vals(k,i) = plot_T(k); i=i+1
            vals(k,i) = plot_v(k); i=i+1
            vals(k,i) = plot_csound(k); i=i+1
            vals(k,i) = plot_etot(k); i=i+1
            vals(k,i) = plot_ekin(k); i=i+1
            vals(k,i) = plot_eeos(k); i=i+1
            vals(k,i) = plot_Peos(k); i=i+1
            vals(k,i) = imex1(k); i=i+1
            vals(k,i) = imex2(k); i=i+1
            vals(k,i) = imex3(k); i=i+1
            vals(k,i) = imex4(k); i=i+1
            vals(k,i) = imex5(k); i=i+1
            vals(k,i) = imex6(k); i=i+1
         end do
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         extras_finish_step = terminate
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call do_imex(s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step
      
      
      subroutine do_imex(s, ierr)
         use imex, only: &
            start_imex, steps_imex, finish_imex, &
            max_calls, nsteps_per_call, &
            model_number_for_nsteps_per_call2, nsteps_per_call2
         use imex_work, only: total_newton_iters, total_gmres_matvecs
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: nsteps_taken, i, nz, nsteps
         real(dp) :: time, dt, total_energy_initial
         integer(8) :: time0, time1, clock_rate
         real(dp) :: runtime
         logical :: final_step, must_write_files
         include 'formats'
         if (max_calls <= 0) then
            ierr = -1
            write(*,*) 'imex max_calls <= 0', max_calls
            return
         end if
         ierr = 0
         call system_clock(time0,clock_rate)
         call start_imex(total_energy_initial, ierr)
         if (ierr /= 0) stop 'start_imex error'
         do i = 1, max_calls 
            if (s% model_number >= model_number_for_nsteps_per_call2) then
               nsteps = nsteps_per_call2
            else
               nsteps = nsteps_per_call
            end if
            call steps_imex( &
               nsteps, time, dt, &
               final_step, nsteps_taken, nz, ierr)
            if (ierr /= 0) stop 'steps_imex error'
            s% nz = nz
            s% model_number = nsteps_taken
            s% time = time
            s% star_age = time/secyer
            s% dt = dt
            if (s% job% pgstar_flag) then
               s% RTI_flag = .true. ! for profile_getval
               s% v_flag = .true.; s% u_flag = .false. ! for profile_getval
               call read_pgstar_controls(s, ierr) 
               if (ierr /= 0) stop 'do_imex read_pgstar_controls error'
               if (final_step) write(*,*) 'done'
               must_write_files = final_step .and. s% job% save_pgstar_files_when_terminate
               call update_pgstar_plots(s, must_write_files, ierr)
               if (ierr /= 0) stop 'do_imex update_pgstar_plots error'
            end if
            if (final_step) exit
         end do
         call finish_imex(ierr)
         if (ierr /= 0) stop 'finish_imex error'
         call system_clock(time1,clock_rate)
         runtime = real(time1 - time0, dp) / clock_rate / 60
         write(*,'(/,a,f12.2,99i14/)') 'runtime (minutes), steps newton_iters gmres_matvecs ', &
              runtime, s% model_number, total_newton_iters, total_gmres_matvecs
         ierr = -1 ! to terminate the star run         
      end subroutine do_imex


      end module run_star_extras
      
