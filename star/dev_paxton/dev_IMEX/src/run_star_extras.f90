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
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
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
         use imex, only: start_imex, steps_imex, finish_imex, initialize_problem_A
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         integer :: nsteps_per_call, nsteps_taken, max_calls, i, nz
         real(dp) :: time, dt
         logical :: final_step, must_write_files
         include 'formats'
         nz = s% nz
         call start_imex(initialize_problem_A, ierr)
         if (ierr /= 0) stop 'start_imex error'
         nsteps_taken = 0
         time = 0d0
         max_calls = 10000
         nsteps_per_call = 500
         do i = 1, max_calls
            call steps_imex( &
               nsteps_per_call, time, dt, &
               final_step, nsteps_taken, nz, ierr)
            if (ierr /= 0) stop 'steps_imex error'
            s% model_number = nsteps_taken
            s% time = time
            s% star_age = time/secyer
            call update_star_data
            if (ierr /= 0) stop 'do_imex update_star_data error'
            if (s% job% pgstar_flag) then
               s% RTI_flag = .true. ! for profile_getval
               s% u_flag = .true. ! for profile_getval
               call read_pgstar_controls(s, ierr) 
               if (ierr /= 0) stop 'do_imex read_pgstar_controls error'
               must_write_files = final_step .and. s% job% save_pgstar_files_when_terminate
               call update_pgstar_plots(s, must_write_files, ierr)
               if (ierr /= 0) stop 'do_imex update_pgstar_plots error'
            end if
            if (final_step) exit
         end do
         call finish_imex(ierr)
         if (ierr /= 0) stop 'finish_imex error'
         !stop 'do_imex'
         ierr = -1 ! to terminate the star run
         
         contains
         
         subroutine update_star_data
            integer k
            call get1_vector(0, s% r, ierr); if (ierr /= 0) return
            call get1_vector(1, s% m, ierr); if (ierr /= 0) return
            call get1_vector(2, s% v, ierr); if (ierr /= 0) return
            do k=1,nz
               s% u_face_ad(k)%val = s% v(k)
            end do
            !call get1_vector(3, s% grav_potential, ierr); if (ierr /= 0) return
            call get1_vector(4, s% dm, ierr); if (ierr /= 0) return
            !call get1_vector(5, s% dr, ierr); if (ierr /= 0) return
            call get1_vector(6, s% rho, ierr); if (ierr /= 0) return
            call get1_vector(7, s% Peos, ierr); if (ierr /= 0) return
            call get1_vector(8, s% u, ierr); if (ierr /= 0) return
            !call get1_vector(9, s% X, ierr); if (ierr /= 0) return
            !call get1_vector(10, s% alpha_RTI, ierr); if (ierr /= 0) return
            !call get1_vector(11, s% cell_mass, ierr); if (ierr /= 0) return
            !call get1_vector(12, s% cell_momentum, ierr); if (ierr /= 0) return
            !call get1_vector(13, s% cell_total_energy, ierr); if (ierr /= 0) return
            !call get1_vector(14, s% cell_mass_X, ierr); if (ierr /= 0) return
            call get1_vector(15, s% csound, ierr); if (ierr /= 0) return
            do k=1,nz
               s% v_div_csound(k) = s% u(k)/s% csound(k)
            end do
         end subroutine update_star_data
         
         subroutine get1_vector(selector, d, ierr)
            integer, intent(in) :: selector
            real(dp) :: d(:)
            integer, intent(out) :: ierr
            ierr = 0
            !ipar(1) = selector
            !ierr = f_dispatch_c(3, &
            !   lrpar, rpar(1), lipar, ipar(1), s% nz, d(1), lidat, idat(1))
         end subroutine get1_vector
         
      end subroutine do_imex


      end module run_star_extras
      
