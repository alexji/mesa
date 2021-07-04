      program run
      use run_star_support, only: do_read_star_job
      use run_star, only: do_run_star
      use star_lib, only: star_test_abtilu
      
      implicit none
      
      integer :: ierr
      
      call star_test_abtilu()
      
      ierr = 0
      call do_read_star_job('inlist', ierr)
      if (ierr /= 0) stop 1
      
      call do_run_star
      
      end program
