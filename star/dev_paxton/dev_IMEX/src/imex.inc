
      subroutine read_inlist_imex_controls(ierr)
         use imex_input
         integer, intent(out) :: ierr
         character (len=256) :: filename, message
         integer :: unit         
         filename = 'inlist_imex_problem'
         ! set defaults
         include "imex_controls_defaults"
         write(*,*) 'read ' // trim(filename)
         open(newunit=unit, file=trim(filename), action='read', delim='quote', iostat=ierr)
         if (ierr /= 0) then
            write(*, *) 'Failed to open control namelist file ', trim(filename)
         else
            read(unit, nml=imex_controls, iostat=ierr)  
            close(unit)
            if (ierr /= 0) then
               write(*, *) 'Failed while trying to read control namelist file ', trim(filename)
               write(*, '(a)') &
                  'The following runtime error message might help you find the problem'
               write(*, *) 
               open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=imex_controls)
               close(unit)
            end if  
         end if

      end subroutine read_inlist_imex_controls
