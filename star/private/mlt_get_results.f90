! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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


      module mlt_get_results

      use star_private_def
      use const_def
      use num_lib
      use utils_lib
      use auto_diff_support
      use star_utils

      implicit none
      
      private
      public :: get_gradT, do1_mlt_eval, Get_results

      contains
      
      
      subroutine get_gradT(s, MLT_option, & ! used to create models
            r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
            iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
            mixing_type, gradT, Y_face, conv_vel, D, Gamma, ierr)
         type (star_info), pointer :: s
         character (len=*), intent(in) :: MLT_option
         real(dp), intent(in) :: &
            r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
            XH1, cgrav, m, gradL_composition_term, mixing_length_alpha
         integer, intent(in) :: iso
         real(dp), intent(out) :: gradT, Y_face, conv_vel, D, Gamma
         integer, intent(out) :: mixing_type, ierr 
         type(auto_diff_real_star_order1) :: &
            gradr_ad, grada_ad, scale_height_ad, gradT_ad, Y_face_ad, mlt_vc_ad, D_ad, &
            Gamma_ad, r_ad, L_ad, T_ad, P_ad, opacity_ad, rho_ad, dV_ad, chiRho_ad, chiT_ad, Cp_ad
         ierr = 0
         r_ad = r
         L_ad = L
         T_ad = T
         P_ad = P
         opacity_ad = opacity
         rho_ad = rho
         dV_ad = 0d0
         chiRho_ad = chiRho
         chiT_ad = chiT
         Cp_ad = Cp
         gradr_ad = gradr
         grada_ad = grada
         scale_height_ad = scale_height
         call Get_results(s, 0, MLT_option, &
            r_ad, L_ad, T_ad, P_ad, opacity_ad, rho_ad, dV_ad, chiRho_ad, chiT_ad, Cp_ad, &
            gradr_ad, grada_ad, scale_height_ad, &
            iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
            s% alpha_semiconvection, s% thermohaline_coeff, &
            mixing_type, gradT_ad, Y_face_ad, mlt_vc_ad, D_ad, Gamma_ad, ierr)
         gradT = gradT_ad%val
         Y_face = Y_face_ad%val
         conv_vel = mlt_vc_ad%val
         D = D_ad%val
         Gamma = Gamma_ad%val
      end subroutine get_gradT
      
         
      subroutine do1_mlt_eval( &
            s, k, MLT_option, gradL_composition_term, &
            gradr, grada, scale_height, mixing_length_alpha, &
            mixing_type, gradT, Y_face, mlt_vc, D, Gamma, ierr)
         use chem_def, only: ih1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         character (len=*), intent(in) :: MLT_option
         type(auto_diff_real_star_order1), intent(in) :: gradr, grada, scale_height
         real(dp), intent(in) :: gradL_composition_term, mixing_length_alpha
         integer, intent(out) :: mixing_type
         type(auto_diff_real_star_order1), intent(out) :: &
            gradT, Y_face, mlt_vc, D, Gamma
         integer, intent(out) :: ierr 
                 
         real(dp) :: cgrav, m, XH1, gradL_old, grada_face_old
         integer :: iso, old_mix_type
         type(auto_diff_real_star_order1) :: r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp
         include 'formats'
         ierr = 0
         
         cgrav = s% cgrav(k)
         m = s% m_grav(k)
         L = wrap_L_00(s,k)
         T = get_T_face(s,k)
         P = get_Peos_face(s,k)
         r = wrap_r_00(s,k)
         opacity = get_kap_face(s,k)
         rho = get_Rho_face(s,k)
         dV = 1d0/rho - 1d0/s% rho_start(k)
         chiRho = get_ChiRho_face(s,k)
         chiT = get_ChiT_face(s,k)
         Cp = get_Cp_face(s,k)
         iso = s% dominant_iso_for_thermohaline(k)
         XH1 = s% xa(s% net_iso(ih1),k)
         
         if (s% use_other_mlt_results) then
            call s% other_mlt_results(s% id, k, MLT_option, &
               r, L, T, P, opacity, rho, chiRho, chiT, Cp, gradr, grada, scale_height, &
               iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
               s% alpha_semiconvection, s% thermohaline_coeff, &
               mixing_type, gradT, Y_face, mlt_vc, D, Gamma, ierr)
         else         
            call Get_results(s, k, MLT_option, &
               r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, scale_height, &
               iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
               s% alpha_semiconvection, s% thermohaline_coeff, &
               mixing_type, gradT, Y_face, mlt_vc, D, Gamma, ierr)
         end if

      end subroutine do1_mlt_eval


      subroutine Get_results(s, k, MLT_option, &  ! NOTE: k=0 is a valid arg
            r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, scale_height, &
            iso, XH1, cgrav, m, gradL_composition_term, mixing_length_alpha, &
            alpha_semiconvection, thermohaline_coeff, &
            mixing_type, gradT, Y_face, conv_vel, D, Gamma, ierr)
         use star_utils
         type (star_info), pointer :: s
         integer, intent(in) :: k
         character (len=*), intent(in) :: MLT_option
         type(auto_diff_real_star_order1), intent(in) :: &
            r, L, T, P, opacity, rho, dV, chiRho, chiT, Cp, gradr, grada, scale_height
         integer, intent(in) :: iso
         real(dp), intent(in) :: &
            XH1, cgrav, m, gradL_composition_term, &
            mixing_length_alpha, alpha_semiconvection, thermohaline_coeff
         integer, intent(out) :: mixing_type
         type(auto_diff_real_star_order1), intent(out) :: &
            gradT, Y_face, conv_vel, D, Gamma
         integer, intent(out) :: ierr
         
         type(auto_diff_real_star_order1) :: &
            Pr, Pg, grav, scale_height2, Lambda, gradL, beta, Y_guess
         character (len=256) :: message        
         logical ::  test_partials, using_TDC, compare_TDC_to_MLT
         logical, parameter :: report = .false.
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         ierr = 0          
                  
         Pr = crad*pow4(T)/3d0
         Pg = P - Pr
         beta = Pg / P
         gradL = grada + gradL_composition_term ! Ledoux temperature gradient
         Lambda = mixing_length_alpha*scale_height
         grav = cgrav*m/pow2(r)
         conv_vel = 0d0
         if (k > 0) then
            s% SOURCE(k) = 0d0
            s% DAMP(k) = 0d0
            s% DAMPR(k) = 0d0
            s% COUPL(k) = 0d0
            s% tdc_num_iters(k) = 0
         end if

         if (report) then
            write(*,*)
            write(*,4) 'enter Get_results k slvr_itr model gradr grada scale_height ' // trim(MLT_option), &
               k, s% solver_iter, s% model_number, gradr%val, grada%val, scale_height%val
         end if

         ! Initialize with no mixing
         call set_no_mixing('')

         ! Bail if we asked for no mixing, or if parameters are bad.
         if (MLT_option == 'none' .or. beta < 1d-10 .or. mixing_length_alpha <= 0d0) return
         
         ! Bail if more parameters are bad
         if (opacity%val < 1d-10 .or. P%val < 1d-20 .or. T%val < 1d-10 .or. Rho%val < 1d-20 &
               .or. m < 1d-10 .or. r%val < 1d-10 .or. cgrav < 1d-10) return

         ! At this stage we're either using MLT or TDC.
         ! TDC requires the MLT answer as a starting point, so
         ! call MLT.
         call set_MLT

         ! check if this particular k needs to be done with TDC
         using_TDC = s% using_TDC
         if (.not. s% have_mlt_vc) using_TDC = .false.
         if (k <= 0 .or. s%dt <= 0d0) using_TDC = .false.
         if (using_TDC) then
            Y_guess = gradT - gradL
            using_TDC = .not. check_if_can_fall_back_to_MLT(s, k, mixing_length_alpha, Y_guess, &
                                                                  T, rho, Cp, dV, opacity, scale_height, grada, conv_vel)
         end if
         ! Run through assuming no TDC.         
         if (gradr > gradL) then ! convective
            if (report) write(*,3) 'call set_MLT', k, s% solver_iter
            call set_MLT
         else if (gradL_composition_term < 0) then
            if (report) write(*,3) 'call set_thermohaline', k, s% solver_iter
            call set_thermohaline
         else if (gradr > grada) then
            if (report) write(*,3) 'call set_semiconvection', k, s% solver_iter
            call set_semiconvection
         end if         

         if (k > 0) then ! save non-TDC values for debugging
            s% xtra1_array(k) = safe_log10(abs(gradT%val - grada%val))
            s% xtra2_array(k) = gradT%val
            s% xtra3_array(k) = conv_vel%val
         end if
         
         ! need to make use of gradL instead of grada consistent - at least for TDC
         if (using_TDC) then
            Y_guess = gradT - gradL
            compare_TDC_to_MLT = s% compare_TDC_to_MLT
            if (compare_TDC_to_MLT) then
               if (report) write(*,3) 'call do_compare_TDC_to_MLT', k, s% solver_iter
               call do_compare_TDC_to_MLT     
            else
               if (report) write(*,3) 'call set_TDC', k, s% solver_iter
               call set_TDC
            end if
         else if (report) then
            write(*,4) 'not okay_to_use_TDC mxtyp conv_vel', k, s% solver_iter, &
               mixing_type, conv_vel%val
         end if
         
         ! If there's too-little mixing to bother, or we hit a bad value, fall back on no mixing.
         if (D%val < s% remove_small_D_limit .or. is_bad(D%val)) then
            if (report) write(*,2) 'D < s% remove_small_D_limit', k, D%val, s% remove_small_D_limit
            mixing_type = no_mixing
         end if

         ! If we made it all that way and are still not mixing, call set_no_mixing.
         ! This catches places above where we might have thought we'd have mixing but
         ! ended up falling back on no mixing. It also reports the correct message.
         if (mixing_type == no_mixing) call set_no_mixing('final mixing_type == no_mixing')
         
         contains

         subroutine set_no_mixing(str)
            character (len=*) :: str
            include 'formats'            
            if (report .and. len_trim(str) > 0) &
               write(*,2) 'Get_results set_no_mixing ' // trim(str), k
            mixing_type = no_mixing
            gradT = gradr
            Y_face = gradT - gradL
            conv_vel = 0d0
            D = 0d0
            Gamma = 0d0
         end subroutine set_no_mixing

         subroutine do_compare_TDC_to_MLT
            ! for compare_TDC_to_MLT
            integer :: std_mixing_type
            type(auto_diff_real_star_order1) :: std_Y_face, c0, L0, A0, &
               std_gradT, std_gradr, std_conv_vel, std_D, std_Gamma, std_scale_height
            include 'formats'         
            std_mixing_type = mixing_type
            std_Y_face = gradT - gradL
            std_gradT = gradT
            std_gradr = gradr
            std_conv_vel = conv_vel
            std_D = D
            std_Gamma = Gamma
            std_scale_height = scale_height
            if (report) then
               write(*,*)
               write(*,1) 'do_compare_TDC_to_MLT'
            end if            
            call set_TDC
            if (D%val < s% remove_small_D_limit .or. is_bad(D%val)) then
               if (report) write(*,2) 'TDC D < s% remove_small_D_limit', k, D%val, s% remove_small_D_limit
               mixing_type = no_mixing
            end if
            if (report .or. s% x_integer_ctrl(19) <= 0) then
               if (std_mixing_type /= mixing_type .or. &
                   abs(std_gradT%val - gradT%val) > 1d-2) then
                  write(*,6) 'k solvr_iter model std_mxtyp tdc_mtyp', &
                     k, s% solver_iter, s% model_number, std_mixing_type, mixing_type
                  write(*,4) 'std_gradT tdc_gradT', k, s% solver_iter, s% model_number, std_gradT%val, gradT%val
                  write(*,4) 'std_vc tdc_vc', k, s% solver_iter, s% model_number, std_conv_vel%val, conv_vel%val
                  write(*,4) 'gradL', k, s% solver_iter, s% model_number, gradL%val
                  write(*,4) 'gradr', k, s% solver_iter, s% model_number, gradr%val
                  write(*,4) 'Y_guess', k, s% solver_iter, s% model_number, Y_guess%val
                  write(*,4) 'std Y', k, s% solver_iter, s% model_number, std_Y_face%val
                  write(*,4) 'tdc Y', k, s% solver_iter, s% model_number, Y_face%val
                  write(*,4) 'std gradT-gradr', k, s% solver_iter, s% model_number, std_gradT%val - std_gradr%val
                  write(*,4) 'tdc gradT-gradr', k, s% solver_iter, s% model_number, gradT%val - gradr%val
                  write(*,4) 'L', k, s% solver_iter, s% model_number, L%val
                  c0 = mixing_length_alpha*rho*T*Cp*4d0*pi*pow2(r)
                  L0 = (16d0*pi*crad*clight/3d0)*cgrav*m*pow4(T)/(P*opacity) ! assumes QHSE for dP/dm
                  if (s% okay_to_set_mlt_vc) then ! is also ok to use mlt_vc_old   
                     A0 = s% mlt_vc_old(k)/sqrt_2_div_3
                  else
                     A0 = s% mlt_vc(k)/sqrt_2_div_3
                  end if
                  write(*,4) 'c0', k, s% solver_iter, s% model_number, c0%val
                  write(*,4) 'L0', k, s% solver_iter, s% model_number, L0%val
                  write(*,4) 'A0', k, s% solver_iter, s% model_number, A0%val
                  write(*,4) 'Q(Y_tdc) for vc=0', k, s% solver_iter, s% model_number, &
                     (L%val - L0%val*grada%val) - L0%val*Y_face%val
                  write(*,4) 'Q(Y_guess) for vc=0', k, s% solver_iter, s% model_number, &
                     (L%val - L0%val*grada%val) - L0%val*Y_guess%val
                  write(*,4) 'gradL_composition_term', k, s% solver_iter, s% model_number, gradL_composition_term
                  write(*,4) 'COUPL', k, s% solver_iter, s% model_number, s% COUPL(k)
                  write(*,4) 'SOURCE', k, s% solver_iter, s% model_number, s% SOURCE(k)
                  write(*,4) 'DAMP', k, s% solver_iter, s% model_number, s% DAMP(k)
                  write(*,*)
                  stop 'do_compare_TDC_to_MLT'
               end if
            end if
         end subroutine do_compare_TDC_to_MLT

         subroutine set_TDC
            include 'formats'
            if (k == 1) then
               call set_no_mixing('set_TDC')
            else
               call get_TDC_solution(s, k, &
                  mixing_length_alpha, cgrav, m, Y_guess, report, &
                  mixing_type, L, r, P, T, rho, dV, Cp, opacity, &
                  scale_height, gradL, grada, conv_vel, Y_face, ierr)
               if (ierr /= 0) then
                  write(*,2) 'get_TDC_solution failed in set_TDC', k
                  write(*,*) 'Repeating call with reporting on.'
                  call get_TDC_solution(s, k, &
                     mixing_length_alpha, cgrav, m, Y_guess, .true., &
                     mixing_type, L, r, P, T, rho, dV, Cp, opacity, &
                     scale_height, gradL, grada, conv_vel, Y_face, ierr)
                  stop 'get_TDC_solution failed in set_TDC'
               end if
            end if
            gradT = Y_face + grada
         end subroutine set_TDC        
         
         subroutine set_MLT
            real(dp) :: ff1, ff2, ff3, ff4
            type(auto_diff_real_star_order1) :: &
               Q, omega, a0, ff4_omega2_plus_1, A_1, A_2, &
               A_numerator, A_denom, A, Bcubed, delta, Zeta, &
               f, f0, f1, f2, radiative_conductivity, convective_conductivity
            include 'formats' 
            Q = chiT/chiRho ! 'Q' param  C&G 14.24
            if (MLT_option == 'Cox' .or. MLT_option == 'TDC') then ! this assumes optically thick
               a0 = 9d0/4d0
               convective_conductivity = &
                  Cp*grav*pow2(Lambda)*rho*(sqrt(Q*rho/(2d0*P)))/9d0 ! erg / (K cm sec)
               radiative_conductivity = &
                  (4d0/3d0*crad*clight)*pow3(T)/(opacity*rho) ! erg / (K cm sec)
               A = convective_conductivity / radiative_conductivity !  unitless.
            else
               select case(trim(MLT_option))
               case ('Henyey')
                  ff1=1.0d0/s% Henyey_MLT_nu_param
                  ff2=0.5d0 
                  ff3=8.0d0/s% Henyey_MLT_y_param
                  ff4=1.0d0/s% Henyey_MLT_y_param
               case ('ML1')
                  ff1=0.125d0 
                  ff2=0.5d0 
                  ff3=24.0d0
                  ff4=0.0d0
               case ('ML2')
                  ff1=1.0d0
                  ff2=2.0d0
                  ff3=16.0d0
                  ff4=0.0d0
               case ('Mihalas')
                  ff1=0.125d0 
                  ff2=0.5d0 
                  ff3=16.0d0
                  ff4=2.0d0
               case default
                  write(*,'(3a)') 'Error: ', trim(MLT_option), ' is not an allowed MLT option'
                  call mesa_error(__FILE__,__LINE__)
               end select
               omega = Lambda*rho*opacity
               ff4_omega2_plus_1 = ff4/pow2(omega) + 1d0
               a0 = (3d0/16d0)*ff2*ff3/ff4_omega2_plus_1
               A_1 = 4d0*Cp*sqrt(ff1*P*Q*rho)
               A_2 = mixing_length_alpha*omega*ff4_omega2_plus_1
               A_numerator = A_1*A_2
               A_denom = ff3*crad*clight*pow3(T)
               A = A_numerator/A_denom   
            end if  
            ! 'B' param  C&G 14.81
            Bcubed = (pow2(A)/a0)*(gradr - gradL)   

            if (Bcubed <= 0d0) then
               ! Radiative zone, because this means that gradr < gradL
               Zeta = 0d0
            else
               ! Convection zone

               ! now solve cubic equation for convective efficiency, Gamma
               ! a0*Gamma^3 + Gamma^2 + Gamma - a0*Bcubed == 0   C&G 14.82, 
               ! leave it to Mathematica to find an expression for the root we want      
               delta = a0*Bcubed               
               f = -2d0 + 9d0*a0 + 27d0*a0*a0*delta
               if (f > 1d100) then
                  f0 = f
               else
                  f0 = pow2(f) + 4d0*(-1d0 + 3d0*a0)*(-1d0 + 3d0*a0)*(-1d0 + 3d0*a0)
                  if (f0 <= 0d0) then
                     f0 = f
                  else
                     f0 = sqrt(f0)         
                  end if
               end if
               f1 = -2d0 + 9d0*a0 + 27d0*a0*a0*delta + f0  
               if (f1 <= 0d0) return
               f1 = pow(f1,one_third)     
               f2 = 2d0*two_13*(1d0 - 3d0*a0) / f1       
               Gamma = (four_13*f1 + f2 - 2d0) / (6d0*a0)
               if (Gamma <= 0d0) return

               ! average convection velocity   C&G 14.86b
               conv_vel = mixing_length_alpha*sqrt(Q*P/(8d0*rho))*Gamma / A
               D = conv_vel*Lambda/3d0     ! diffusion coefficient [cm^2/sec]

               !Zeta = pow3(Gamma)/Bcubed  ! C&G 14.80
               Zeta = exp(3d0*log(Gamma) - log(Bcubed)) ! write it this way to avoid overflow problems
            end if
            
            ! Zeta must be >= 0 and <= 1
            if (is_bad(Zeta%val)) return
            if (Zeta < 0d0) then
               Zeta = 0d0
            else if (Zeta > 1d0) then
               Zeta = 1d0
            end if            
            
            
            if (.false.) then ! need this old form for a few test cases
               gradT = (1d0 - Zeta)*gradr + Zeta*grada ! C&G 14.79      
               Y_face = gradT - grada
            else ! switch to this when resolve the problems with those cases
               gradT = (1d0 - Zeta)*gradr + Zeta*gradL ! C&G 14.79      
               Y_face = gradT - gradL
            end if
            
            
            mixing_type = convective_mixing
            if (k > 0) then
               s% xtra1_array(k) = conv_vel%val
               s% xtra2_array(k) = gradT%val
            end if

            if (report) then
               write(*,2) 'set_MLT val for Zeta gradr grada gradT Y_face', k, &
                  Zeta%val, gradr%val, grada%val, gradT%val, Y_face%val
               write(*,2) 'set_MLT d_dlnd_00 for Zeta gradr grada gradT Y_face', k, &
                  Zeta%d1Array(i_lnd_00), gradr%d1Array(i_lnd_00), &
                  grada%d1Array(i_lnd_00), gradT%d1Array(i_lnd_00), &
                  Y_face%d1Array(i_lnd_00)
            end if

         end subroutine set_MLT   

!------------------------------ Semiconvection

         subroutine set_semiconvection ! Langer 1983 & 1985
            type(auto_diff_real_star_order1) :: bc, LG, &
               radiative_conductivity, a0, a1, a2, a3, a4, a5, a6, a, &
               b1, b2, b3, b4, b5, b6, b7, b, div, bsq    
            real(dp) :: alpha
            include 'formats'
            radiative_conductivity = &
               (4d0/3d0*crad*clight)*pow3(T)/(opacity*rho) ! erg / (K cm sec)
            D = alpha_semiconvection*radiative_conductivity/(6d0*Cp*rho) &
                  *(gradr - grada)/(gradL - gradr)
            if (D%val <= 0) return         
            if (s% semiconvection_option == 'Langer_85 mixing; gradT = gradr') then
               gradT = gradr
               Y_face = gradT - grada
               conv_vel = 3d0*D/Lambda             
               mixing_type = semiconvective_mixing
               return
            end if            
            if (s% semiconvection_option /= 'Langer_85') then
               write(*,*) 'MLT: unknown values for semiconvection_option ' // &
                  trim(s% semiconvection_option)
               ierr = -1
               return
            end if            
!            Solve[{
!                  L/Lrad - Lsc/Lrad - 1 == 0, 
!                  Lrad == grad LG, 
!                  gradMu == (4 - 3*beta)/beta*gradL_composition_term,
!                  Lsc/Lrad == alpha (grad - gradA)/(2 grad (gradL - grad))
!                              (grad - gradA - (beta (8 - 3 beta))/bc gradMu)}, 
!                  grad, {Lsc, Lrad, gradMu}] // Simplify
            alpha = min(1d0, alpha_semiconvection)
            bc = 32d0 - 24d0*beta - beta*beta            
            LG = (16d0*pi*clight*m*cgrav*Pr)/(P*opacity)            
            a0 = alpha*gradL_composition_term*LG            
            a1 = -2d0*bc*L            
            a2 = 2d0*alpha*bc*grada*LG            
            a3 = -2d0*bc*gradL*LG            
            a4 = 32d0*a0            
            a5 = -36d0*beta*a0            
            a6 = 9d0*beta*beta*a0            
            a = a1 + a2 + a3 + a4 + a5 + a6                           
            b1 = 32d0 - 36d0*beta + 9d0*beta*beta            
            b2 = b1*a0            
            b3 = -2d0*gradL*L + alpha*grada*grada*LG            
            b4 = (-alpha*gradA + gradL)*LG            
            b5 = -b2 + 2d0*bc*(L + b4)            
            b6 = b2*grada + bc*b3            
            b7 = -4d0*(alpha - 2d0)*bc*LG*b6            
            b = b7 + b5*b5            
            div = 2d0*(alpha - 2d0)*bc*LG
            bsq = sqrt(b)
            gradT = (a + bsq)/div
            Y_face = gradT - grada
            conv_vel = 3d0*D/Lambda             
            mixing_type = semiconvective_mixing
         end subroutine set_semiconvection
        
         subroutine set_thermohaline
            real(dp) :: D_thrm
            call get_D_thermohaline(s, &
               grada%val, gradr%val, T%val, opacity%val, rho%val, &
               Cp%val, gradL_composition_term, &
               iso, XH1, thermohaline_coeff, D_thrm, ierr)
            D = D_thrm
            gradT = gradr
            Y_face = gradT - grada
            conv_vel = 3d0*D/Lambda
            mixing_type = thermohaline_mixing 
         end subroutine set_thermohaline

      end subroutine Get_results
      
!------------------------------ Time-dependent convection (TDC)

      
      subroutine get_TDC_solution(s, k, &
            mixing_length_alpha, cgrav, m, Y_guess, report, &
            mixing_type, L, r, P, T, rho, dV, Cp, kap, Hp, gradL, grada, cv, Y_face, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: mixing_length_alpha, cgrav, m
         integer, intent(out) :: mixing_type
         type(auto_diff_real_star_order1), intent(in) :: &
            L, r, P, T, rho, dV, Cp, kap, Hp, gradL, grada, Y_guess
         logical, intent(in) :: report
         type(auto_diff_real_star_order1),intent(out) :: cv, Y_face
         integer, intent(out) :: ierr
         
         type(auto_diff_real_star_order1) :: A0, c0, L0
         type(auto_diff_real_tdc) :: Af, Y, Z, Q, Z_new, dQdZ, correction, lower_bound_Z, upper_bound_Z
         type(auto_diff_real_tdc) :: Q_start, Y_start
         real(dp) ::  gradT, Lr, Lc, scale
         integer :: iter
         logical :: converged, Y_is_positive, first_Q_is_positive
         real(dp), parameter :: tolerance = 1d-8
         real(dp), parameter :: alpha_c  = (1d0/2d0)*sqrt_2_div_3
         integer, parameter :: max_iter = 200
         include 'formats'

         ierr = 0
         if (mixing_length_alpha == 0d0 .or. k < 1 .or. s% dt <= 0d0) then
            stop 'bad call to TDC get_TDC_solution'
         end if         

         ! Set up for solve
         c0 = mixing_length_alpha*alpha_c*rho*T*Cp*4d0*pi*pow2(r)
         L0 = (16d0*pi*crad*clight/3d0)*cgrav*m*pow4(T)/(P*kap) ! assumes QHSE for dP/dm
         if (s% okay_to_set_mlt_vc) then
            A0 = s% mlt_vc_old(k)/sqrt_2_div_3
         else
            A0 = s% mlt_vc(k)/sqrt_2_div_3
         end if

         ! Set scale for judging the solution.
         ! Q has units of a luminosity, so the scale should be a luminosity.
         if (s% solver_iter == 0) then
            scale = max(abs(s% L(k)), 1d-3*maxval(s% L(1:s% nz)))
         else
            scale = max(abs(s% L_start(k)), 1d-3*maxval(s% L_start(1:s% nz)))
         end if

         ! First, find a guess for Y.
         !
         ! If Q(Y=0) is positive then the luminosity is too great to be carried radiatively, so
         ! we'll necessarily have Y > 0.
         !
         ! If Q(Y=0) is negative then the luminosity can be carried by radiation alone, so we'll
         ! necessarily have Y < 0.
         !
         ! This isn't just a physics argument: if unconvinced examine the mathematical form of Q and notice that it's
         ! monotonically decreasing in Y, so if Q > 0 when Y=0 the solution to Q=0 has Y > 0. Likewise if Q < 0 when Y=0
         ! the solution to Q=0 has Y < 0.
         converged = .false.
         Y = 0d0
         call compute_Q(s, k, mixing_length_alpha, &
            Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q, Af)
         if (abs(Q / scale) < tolerance) then
            converged = .true.
         else
            if (Q > 0d0) then
               Y_is_positive = .true.
               Y = convert(abs(Y_guess))
            else
               Y_is_positive = .false.
               Y = convert(-abs(Y_guess))
            end if

            ! Newton's method to find solution Y
            Y%d1val1 = Y%val ! Fill in starting dY/dZ. Using Y = \pm exp(Z) we find dY/dZ = Y.
            Z = log(abs(Y))

            ! Save starting values
            Q_start = Q
            Y_start = Y

            ! We use the fact that Q(Y) is monotonic for Y > 0 to produce iteratively refined bounds on Q.
            lower_bound_Z = -100d0
            upper_bound_Z = 50d0 
            
            if (report) write(*,2) 'initial Z', 0, Z%val
            do iter = 1, max_iter
               call compute_Q(s, k, mixing_length_alpha, &
                  Y, c0, L, L0, A0, T, rho, dV, Cp, kap, Hp, gradL, grada, Q, Af)
               if (is_bad(Q%val)) exit
               if (abs(Q%val)/scale <= tolerance) then
                  if (report) write(*,2) 'converged', iter, abs(Q%val)/scale, tolerance
                  converged = .true.
                  exit
               end if

               if (gradL > 0d0 .and. Y_is_positive) then ! Take advantage of monotonicity when we can...
                  if (Q > 0d0) then
                     ! Q(Y) is monotonic so this means Z is a lower-bound.
                     lower_bound_Z = Z
                  else
                     ! Q(Y) is monotonic so this means Z is an upper-bound.
                     upper_bound_Z = Z
                  end if
               end if

               dQdZ = differentiate_1(Q)

               if (is_bad(dQdZ%val) .or. abs(dQdZ%val) < 1d-99) then
                  if (report) write(*,2) 'dQdZ', iter, dQdZ%val
                  exit
               end if

               correction = -Q/dQdz
               if (correction > 2d0) then
                  correction = 2d0
               else if (correction < -2d0) then
                  correction = -2d0
               end if

               if (abs(correction) < 1d-13) then
                  ! Can't get much more precision than this.
                  converged = .true.
                  exit
               end if

               Z_new = Z + correction

               ! If the correction pushes the solution out of bounds then we know
               ! that was a bad step. Bad steps are still in the same direction, they just
               ! go too far, so we replace that result with one that's halfway to the relevant bound.
               if (Z_new > upper_bound_Z) then
                  Z_new = (Z + upper_bound_Z) / 2d0
               else if (Z_new < lower_bound_Z) then
                  Z_new = (Z + lower_bound_Z) / 2d0
               end if



               if (report) write(*,2) 'iter Z_new, Z, low_bnd, upr_bnc, Q/dQdZ, Q, dQdZ, corr', iter, &
                  Z_new%val, Z%val, lower_bound_Z%val, upper_bound_Z%val, Q%val/dQdZ%val, Q%val, dQdZ%val, correction%val
               Z_new%d1val1 = 1d0            
               Z = Z_new

               if (Y_is_positive) then
                  Y = exp(Z)
               else
                  Y = -exp(Z)
               end if

            end do
         end if

         if (.not. converged) then
            ierr = 1
            if (report .or. s% x_integer_ctrl(19) <= 0) then
            !$OMP critical (tdc_crit0)
               write(*,5) 'failed get_TDC_solution k slvr_iter model TDC_iter', &
                  k, s% solver_iter, s% model_number, iter
               write(*,*) 'Q_start', Q_start%val
               write(*,*) 'Y_start', Y_start%val
               write(*,2) 'Q', k, Q%val
               write(*,2) 'scale', k, scale
               write(*,2) 'Q/scale', k, Q%val/scale
               write(*,2) 'tolerance', k, tolerance
               write(*,2) 'dQdZ', k, dQdZ%val
               write(*,2) 'Y', k, Y%val
               write(*,2) 'dYdZ', k, Y%d1val1
               write(*,2) 'exp(Z)', k, exp(Z%val)
               write(*,2) 'Z', k, Z%val
               write(*,2) 'Af', k, Af%val
               write(*,2) 'dAfdZ', k, Af%d1val1
               write(*,2) 'A0', k, A0%val
               write(*,2) 'c0', k, c0%val
               write(*,2) 'L', k, L%val
               write(*,2) 'L0', k, L0%val
               write(*,2) 'grada', k, grada%val
               write(*,2) 'gradL', k, gradL%val
               write(*,*)
            !$OMP end critical (tdc_crit0)
            end if
            return
         end if

         cv = sqrt_2_div_3*unconvert(Af)   
         Y_face = unconvert(Y)
         gradT = Y_face%val + gradL%val
         Lr = L0%val*gradT
         Lc = L%val - Lr
         if (cv > 0d0) then
            mixing_type = convective_mixing
         else
            mixing_type = no_mixing
         end if
         if (k > 0) s% tdc_num_iters(k) = iter          
      end subroutine get_TDC_solution
            

      !> Q is the residual in the TDC equation, namely:
      !!
      !! Q = (L - L0 * gradL) - (L0 + c0 * Af) * Y
      !!
      !! @param s star pointer
      !! @param k face index
      !! @param Y superadiabaticity
      !! @param c0_in A proportionality factor for the convective luminosity
      !! @param L_in luminosity
      !! @param L0_in L0 = (Lrad / grad_rad) is the luminosity radiation would carry if dlnT/dlnP = 1.
      !! @param A0 Initial convection speed
      !! @param T Temperature
      !! @param rho Density (g/cm^3)
      !! @param dV 1/rho_face - 1/rho_start_face (change in specific volume at the face)
      !! @param Cp Heat capacity
      !! @param kap Opacity
      !! @param Hp Pressure scale height
      !! @param gradL_in gradL is the neutrally buoyant dlnT/dlnP (= grad_ad + grad_mu),
      !! @param grada_in grada is the adiabatic dlnT/dlnP,
      !! @param Q The residual of the above equaiton (an output).
      !! @param Af The final convection speed (an output).
      subroutine compute_Q(s, k, mixing_length_alpha, &
            Y, c0_in, L_in, L0_in, A0, T, rho, dV, Cp, kap, Hp, gradL_in, grada_in, Q, Af)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: mixing_length_alpha
         type(auto_diff_real_star_order1), intent(in) :: &
            c0_in, L_in, L0_in, A0, T, rho, dV, Cp, kap, Hp, gradL_in, grada_in
         type(auto_diff_real_tdc), intent(in) :: Y
         type(auto_diff_real_tdc), intent(out) :: Q, Af
         type(auto_diff_real_tdc) :: xi0, xi1, xi2, c0, L0, L, gradL

         call eval_xis(s, k, mixing_length_alpha, &
            Y, T, rho, Cp, dV, kap, Hp, grada_in, xi0, xi1, xi2)          

         Af = eval_Af(s, k, A0, xi0, xi1, xi2)
         L = convert(L_in)
         L0 = convert(L0_in)
         gradL = convert(gradL_in)
         c0 = convert(c0_in)
         Q = (L - L0*gradL) - (L0 + c0*Af)*Y

      end subroutine compute_Q

      !> Determines if it is safe (physically) to use MLT instead of TDC.
      !!
      !! The TDC velocity solution approaches the MLT solution asymptotically like exp(-Jt),
      !! so when Jt >> 1 we can fall back to MLT safely.
      !!
      !! Likewise in cases where the predicted difference between Af and A0 (i.e. change in convection speed)
      !! is small we can fall back safely to MLT.
      !!
      !! For both tests we need to know Y (the superadiabaticity), because that's what TDC calculates,
      !! so we use the Y value from MLT as a guess.
      !!
      !! @param s star pointer
      !! @param k face index
      !! @param mixing_length_alpha mixing length parameter
      !! @param Y_in guess of superadiabaticity
      !! @param T temperature (K)
      !! @param rho density (g/cm^3)
      !! @param Cp heat capacity (erg/g/K)
      !! @param kap opacity (cm^2/g)
      !! @param Hp pressure scale height (cm)
      !! @param grada is the adiabatic dlnT/dlnP.
      !! @param A0 convection speed from the start of the step (cm/s)
      !! @param fallback False if we need to use TDC, True if we can fall back to MLT.
      logical function check_if_can_fall_back_to_MLT(s, k, mixing_length_alpha, &
                                                      Y_in, T, rho, Cp, dV, kap, Hp, grada, A0) result(fallback)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: mixing_length_alpha
         type(auto_diff_real_star_order1), intent(in) :: Y_in, T, rho, Cp, dV, kap, Hp, grada, A0
         type(auto_diff_real_tdc) :: Y, Af, xi0, xi1, xi2, J2, Jt2

         if (abs(s%mstar_dot) > 1d-99 .and. k < s% k_const_mass) then
            fallback = .true.
            return
         end if

         Y = convert(Y_in)

         call eval_xis(s, k, mixing_length_alpha, &
            Y, T, rho, Cp, dV, kap, Hp, grada, xi0, xi1, xi2)

         Af = eval_Af(s, k, A0, xi0, xi1, xi2)

         J2 = pow2(xi1) - 4d0 * xi0 * xi2
         Jt2 = abs(J2) * pow2(s%dt)

         ! Note the '1d-50' here is in cm/s, and is just there to avoid division by zero.
         if (Jt2 > 3d3 .or. abs(Af%val - A0%val) / (1d-50 + abs(A0%val)) < 1d-8) then
            fallback = .true.
         else
            fallback = .false.
         end if

      end function check_if_can_fall_back_to_MLT

      !! Calculates the coefficients of the TDC velocity equation.
      !! The velocity equation is
      !!
      !! 2 dw/dt = xi0 + w * xi1 + w**2 * xi2 - Lambda
      !!
      !! where Lambda (currently set to zero) captures coupling between cells.
      !!
      !! The coefficients xi0/1/2 are given by
      !!
      !! xi0 = (epsilon_q + S * Y) / w
      !! xi1 = (-D_r + p_turb dV/dt) / w^2    [here V = 1/rho is the volume]
      !! xi2 = -D / w^3
      !!
      !! Note that these terms are evaluated on faces, because the convection speed
      !! is evaluated on faces. As a result all the inputs must either natively
      !! live on faces or be interpolated to faces.
      !!
      !! This function also has some side effects in terms of storing some of the terms
      !! it calculates for plotting purposes.
      !! TODO: These should be removed or changed from xtra arrays to named arrays.
      !!
      !! @param s star pointer
      !! @param k face index
      !! @param mixing_length_alpha mixing length parameter
      !! @param Y superadiabaticity
      !! @param T temperature (K)
      !! @param rho density (g/cm^3)
      !! @param Cp heat capacity (erg/g/K)
      !! @param dV 1/rho_face - 1/rho_start_face (change in specific volume at the face)
      !! @param kap opacity (cm^2/g)
      !! @param Hp pressure scale height (cm)
      !! @param grada is the adiabatic dlnT/dlnP.
      !! @param xi0 Output, the constant term in the convective velocity equation.
      !! @param xi1 Output, the prefactor of the linear term in the convective velocity equation.
      !! @param xi2 Output, the prefactor of the quadratic term in the convective velocity equation.
      subroutine eval_xis(s, k, mixing_length_alpha, &
            Y, T, rho, Cp, dV, kap, Hp, grada, xi0, xi1, xi2) 
         ! eval_xis sets up Y with partial wrt Z
         ! so results come back with partials wrt Z
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: mixing_length_alpha
         type(auto_diff_real_star_order1), intent(in) :: T, rho, dV, Cp, kap, Hp, grada
         type(auto_diff_real_tdc), intent(in) :: Y
         type(auto_diff_real_tdc), intent(out) :: xi0, xi1, xi2
         type(auto_diff_real_tdc) :: S0, D0, DR0
         type(auto_diff_real_star_order1) :: gammar_div_alfa, Pt0, dVdt
         real(dp), parameter :: x_ALFAS = (1.d0/2.d0)*sqrt_2_div_3
         real(dp), parameter :: x_CEDE  = (8.d0/3.d0)*sqrt_2_div_3
         real(dp), parameter :: x_ALFAP = 2.d0/3.d0
         real(dp), parameter :: x_GAMMAR = 2.d0*sqrt(3.d0)
         include 'formats'

         S0 = convert(x_ALFAS*mixing_length_alpha*Cp*T*grada/Hp)
         S0 = S0*Y
         D0 = convert(s% alpha_TDC_DAMP*x_CEDE/(mixing_length_alpha*Hp))
         if (s% alpha_TDC_DAMPR == 0d0) then
            DR0 = 0d0
         else
            gammar_div_alfa = s% alpha_TDC_DAMPR*x_GAMMAR/(mixing_length_alpha*Hp)
            DR0 = convert(4d0*boltz_sigma*pow2(gammar_div_alfa)*pow3(T)/(pow2(rho)*Cp*kap))
         end if
         Pt0 = s% alpha_TDC_PtdVdt*x_ALFAP*rho
         if (s% dt > 0) then
            dVdt = dV/s% dt
         else
            dVdt = 0d0
         end if
         xi0 = S0
         xi1 = -(DR0 + convert(Pt0*dVdt))
         xi2 = -D0

         if (k > 0) then   
            s% xtra4_array(k) = S0%val
            s% xtra5_array(k) = D0%val
            s% xtra6_array(k) = DR0%val
         end if
      end subroutine eval_xis
      
      !! Calculates the solution to the TDC velocity equation.
      !! The velocity equation is
      !!
      !! 2 dw/dt = xi0 + w * xi1 + w**2 * xi2 - Lambda
      !!
      !! where Lambda (currently set to zero) captures coupling between cells.
      !! The xi0/1/2 variables are constants for purposes of solving this equation.
      !!
      !! An important related parameter is J:
      !! 
      !! J^2 = xi1^2 - 4 * xi0 * xi2
      !!
      !! When J^2 > 0 the solution for w is hyperbolic in time.
      !! When J^2 < 0 the solution is trigonometric, behaving like tan(J * t / 4 + c).
      !!
      !! In the trigonometric branch note that once the solution passes through its first
      !! root the solution for all later times is w(t) = 0. This is not a solution to the
      !! convective velocity equation as written above, but *is* a solution to the convective
      !! energy equation (multiply the velocity equation by w), which is the one we actually
      !! want to solve (per Radek Smolec's thesis). We just solve the velocity form
      !! because it's more convenient.
      !!
      !! This function also has some side effects in terms of storing some of the terms
      !! it calculates for plotting purposes.
      !! TODO: These should be removed or changed from xtra arrays to named arrays.
      !!
      !! @param s star pointer
      !! @param k face index
      !! @param A0_in convection speed from the start of the step (cm/s)
      !! @param xi0 The constant term in the convective velocity equation.
      !! @param xi1 The prefactor of the linear term in the convective velocity equation.
      !! @param xi2 The prefactor of the quadratic term in the convective velocity equation.            
      !! @param Af Output, the convection speed at the end of the step (cm/s)
      function eval_Af(s, k, A0_in, xi0, xi1, xi2) result(Af)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(in) :: A0_in
         type(auto_diff_real_tdc), intent(in) :: xi0, xi1, xi2
         type(auto_diff_real_tdc) :: Af ! output
         type(auto_diff_real_tdc) :: J2, J, Jt, Jt4, num, den, y_for_atan, root, A0, lk 
         real(dp) :: dt    

         ! Debugging
         logical :: dbg = .false.
         integer :: dbg_k = 1448  

         include 'formats'

         J2 = pow2(xi1) - 4d0 * xi0 * xi2
         dt = s%dt
         A0 = convert(A0_in)

         if (J2 > 0d0) then ! Hyperbolic branch
            J = sqrt(J2)
            Jt = dt * J
            Jt4 = 0.25d0 * Jt
            num = safe_tanh(Jt4) * (2d0 * xi0 + A0 * xi1) + A0 * J
            den = safe_tanh(Jt4) * (xi1 + 2d0 * A0 * xi2) - J
            Af = num / den 
            if (Af < 0d0) then
               Af = -Af
            end if
         else if (J2 < 0d0) then ! Trigonometric branch
            J = sqrt(-J2)
            Jt = dt * J

            ! This branch contains decaying solutions that reach A = 0, at which point
            ! they switch onto the 'zero' branch. So we have to calculate the position of
            ! the first root to check it against dt.
            y_for_atan = xi1 + 2d0 * A0 * xi2
            root = safe_atan(J, xi1) - safe_atan(J, y_for_atan)

            ! The root enters into a tangent, so we can freely shift it by pi and
            ! get another root. We care about the first positive root, and the above prescription
            ! is guaranteed to give an answer between (-2*pi,2*pi) because atan produces an answer in [-pi,pi],
            ! so we add/subtract a multiple of pi to get the root into [0,pi).
            if (root > pi) then
               root = root - pi
            else if (root < -pi) then
               root = root + 2d0*pi
            else if (root < 0d0) then
               root = root + pi
            end if

            if (0.25d0 * Jt < root) then
               num = -xi1 + J * tan(0.25d0 * Jt + atan(y_for_atan / J)) 
               den = 2d0 * xi2
               Af = num / den
            else
               Af = 0d0
            end if
         else ! if (J2 == 0d0) then         
            Af = A0            
         end if         
         if (k > 0) then ! save for plots
            s% SOURCE(k) = xi0%val*Af%val
            s% DAMPR(k) = -xi1%val*pow2(Af%val)
            s% DAMP(k) = -xi2%val*pow3(Af%val)
            s% COUPL(k) = s% SOURCE(k) - s% DAMP(k) - s% DAMPR(k)
         end if


         if (k == dbg_k .and. dbg) then
            write(*,*) J2%val,J2%d1val1
            write(*,*) Jt%val,Jt%d1val1
            write(*,*) xi0%val,xi0%d1val1
            write(*,*) xi1%val,xi1%d1val1
            write(*,*) xi2%val,xi2%d1val1
            write(*,*) num%val
            write(*,*) den%val
            write(*,*) Af%val,Af%d1val1
            write(*,*) ''
         end if


      end function eval_Af
      
      !> Computes the hyperbolic tangent of x in a way that is numerically safe.
      !!
      !! @param x Input
      !! @param z Output
      type(auto_diff_real_tdc) function safe_tanh(x) result(z)
         type(auto_diff_real_tdc), intent(in) :: x

         if (x > 50d0) then
            z = 1d0
         else if (x < -50d0) then
            z = -1d0
         else
            z = tanh(x)
         end if
      end function safe_tanh

      !> Computes the arctangent of y/x in a way that is numerically safe near x=0.
      !!
      !! @param x x coordinate for the arctangent.
      !! @param y y coordinate for the arctangent.
      !! @param z Polar angle z such that tan(z) = y / x.
      type(auto_diff_real_tdc) function safe_atan(x,y) result(z)
         type(auto_diff_real_tdc), intent(in) :: x,y
         type(auto_diff_real_tdc) :: x1, y1
         if (abs(x) < 1d-50) then
            ! x is basically zero, so for ~any non-zero y the ratio y/x is ~infinity.
            ! That means that z = +- pi. We want z to be positive, so we return pi.
            z = pi
         else
            z = atan(y/x)
         end if
      end function safe_atan
      
      !> The TDC newton solver needs higher-order partial derivatives than
      !! the star newton solver, because the TDC one needs to pass back a result
      !! which itself contains the derivatives that the star solver needs.
      !! These additional derivatives are provided by the auto_diff_real_tdc type.
      !!
      !! This method converts a auto_diff_real_star_order1 variable into a auto_diff_real_tdc,
      !! setting the additional partial derivatives to zero. This 'upgrades' variables storing
      !! stellar structure to a form the TDC solver can use.
      !!
      !! @param K_in, input, an auto_diff_real_star_order1 variable
      !! @param K, output, an auto_diff_real_tdc variable.
      type(auto_diff_real_tdc) function convert(K_in) result(K)
         type(auto_diff_real_star_order1), intent(in) :: K_in
         K%val = K_in%val
         K%d1Array(1:auto_diff_star_num_vars) = K_in%d1Array(1:auto_diff_star_num_vars)
         K%d1val1 = 0d0
         K%d1val1_d1Array(1:auto_diff_star_num_vars) = 0d0
      end function convert
      
      !> The TDC newton solver needs higher-order partial derivatives than
      !! the star newton solver, because the TDC one needs to pass back a result
      !! which itself contains the derivatives that the star solver needs.
      !! These additional derivatives are provided by the auto_diff_real_tdc type.
      !!
      !! This method converts a auto_diff_real_tdc variable into a auto_diff_real_star_order1,
      !! dropping the additional partial derivatives which (after the TDC solver is done) are
      !! no longer needed. This allows the output of the TDC solver to be passed back to the star solver.
      !!
      !! @param K_in, input, an auto_diff_real_tdc variable
      !! @param K, output, an auto_diff_real_star_order1 variable.      
      type(auto_diff_real_star_order1) function unconvert(K_in) result(K)
         type(auto_diff_real_tdc), intent(in) :: K_in
         K%val = K_in%val
         K%d1Array(1:auto_diff_star_num_vars) = K_in%d1Array(1:auto_diff_star_num_vars)
      end function unconvert


!------------------------------ Thermohaline 

      !> Computes the diffusivity of thermohaline mixing when the
      !! thermal gradient is stable and the composition gradient is unstable.
      !!
      !! @param grada Adiabatic gradient dlnT/dlnP
      !! @param gradr Radiative temperature gradient dlnT/dlnP, equals the actual gradient because there's no convection
      !! @param T Temperature
      !! @param opacity opacity
      !! @param rho Density
      !! @param Cp Heat capacity at constant pressure
      !! @param gradL_composition_term dlnMu/dlnP where Mu is the mean molecular weight.
      !! @param iso The index of the species that drives thermohaline mixing.
      !! @param XH1 Mass fraction of H1.
      !! @param thermohaline_coeff Free parameter multiplying the thermohaline diffusivity.
      !! @param D_thrm Output, diffusivity.
      !! @param ierr Output, error index.
      subroutine get_D_thermohaline(s, &
            grada, gradr, T, opacity, rho, Cp, gradL_composition_term, &
            iso, XH1, thermohaline_coeff, D_thrm, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: &
            grada, gradr, T, opacity, rho, Cp, gradL_composition_term, XH1, &
            thermohaline_coeff
         integer, intent(in) :: iso      
         real(dp), intent(out) :: D_thrm
         integer, intent(out) :: ierr
         real(dp) :: dgrad, K_therm, K_T, K_mu, nu, R0, Pr, tau, r_th            
         include 'formats'     
         dgrad = max(1d-40, grada - gradr) ! positive since Schwarzschild stable               
         K_therm = 4d0*crad*clight*pow3(T)/(3d0*opacity*rho) ! thermal conductivity
         if (s% thermohaline_option == 'Kippenhahn') then
            ! Kippenhahn, R., Ruschenplatt, G., & Thomas, H.-C. 1980, A&A, 91, 175
            D_thrm = -3d0*K_therm/(2*rho*cp)*gradL_composition_term/dgrad
         else if (s% thermohaline_option == 'Traxler_Garaud_Stellmach_11' .or. &
                  s% thermohaline_option == 'Brown_Garaud_Stellmach_13') then
            call get_diff_coeffs(s, &
               K_therm, Cp, rho, T, opacity, iso, XH1, K_T, K_mu, nu)
            R0 = (gradr - grada)/gradL_composition_term
            Pr = nu/K_T
            tau = K_mu/K_T
            r_th = (R0 - 1d0)/(1d0/tau - 1d0)
            if (r_th >= 1d0) then ! stable if R0 >= 1/tau
               D_thrm = 0d0
            else if (Pr < 0d0) then
               ! Bad results from get_diff_coeffs will just result in NaNs from thermohaline options, so skip
               D_thrm = 0d0
            else if (s% thermohaline_option == 'Traxler_Garaud_Stellmach_11') then 
               ! Traxler, Garaud, & Stellmach, ApJ Letters, 728:L29 (2011).
               ! also see Denissenkov. ApJ 723:563–579, 2010.
               D_thrm = 101d0*sqrt(K_mu*nu)*exp(-3.6d0*r_th)*pow(1d0 - r_th,1.1d0) ! eqn 24
            else ! if (s% thermohaline_option == 'Brown_Garaud_Stellmach_13') then
               D_thrm = K_mu*(Numu(R0,r_th,pr,tau) - 1d0)
            endif
         else
            D_thrm = 0
            ierr = -1
            write(*,*) 'unknown for MLT thermohaline_option' // trim(s% thermohaline_option)
         end if
         D_thrm = thermohaline_coeff*D_thrm
      end subroutine get_D_thermohaline


      subroutine get_diff_coeffs(s, &
            K_therm, Cp, rho, T, opacity, iso, XH1, kt, kmu, vis)
         use chem_def, only: chem_isos
         type (star_info), pointer :: s
         real(dp), intent(in) :: K_therm, Cp, rho, T, opacity, XH1
         integer, intent(in) :: iso      
         real(dp), intent(out) :: kt, kmu, vis
         real(dp) :: loglambdah, loglambdacx, loglambdacy, ccx, ccy, qe4
         real(dp) :: Bcoeff, chemA, chemZ, acx, acy, nu_mol, nu_rad      
         real(dp), parameter :: sqrt5 = sqrt(5d0)           
         kt = K_therm/(Cp*rho)       ! thermal diffusivity (assumes radiatively dominated)
         qe4=pow4(qe)
      
         ! Log Lambda for pure H (equation 10 from Proffitt Michaud 93)
         loglambdah = -19.26d0 - 0.5d0*log(rho) + 1.5d0*log(T) - 0.5d0*log(1d0 + 0.5d0*(1+XH1)) 
         nu_rad = 4d0*crad*pow4(T)/(15d0*clight*opacity*pow2(rho)) ! radiative viscosity
         nu_mol = 0.406d0*sqrt(amu)*pow(boltzm*T,2.5d0)/(qe4*loglambdah*rho) 
         ! From Spitzer "Physics of Fully Ionized Gases equation 5-54
         ! Assumes pure H. Still trying to work out what it would be for a mixture. 
         vis = nu_mol + nu_rad   ! total viscosity

         ! The following is from Proffitt & Michaud, 1993.
         ! Their constant B (equation 15)
         Bcoeff = (15.d0/16.d0)*sqrt(2.d0*amu/(5*pi))*pow(boltzm,2.5d0)/qe4
         ! Extract what species drives the thermohaline concvection
         chemA = chem_isos%Z_plus_N(iso)
         chemZ = chem_isos%Z(iso)

         if(chemZ.gt.2) then
         ! This is if the driving chemical is NOT He.
            ! Log Lambda for H-dominant chem mixture (equation 10)
            loglambdacx = loglambdah - log(chemz)  
            ! Log Lambda for He-dominant chem mixture (equation 10)
            loglambdacy = loglambdah - log(2.d0*chemz)
            ! Calculation of C_ij coeffs (equation 12)
            ccx = log(exp(1.2d0*loglambdacx)+1.)/1.2d0
            ccy = log(exp(1.2d0*loglambdacy)+1.)/1.2d0
            ! Reduced masses (I had to guess, from Bahcall & Loeb 1990), with H and He
            acx = (1.d0*chemA)/(1.d0+chemA)
            acy = 4*chemA/(4.d0+chemA)
            ! My formula (see notes) based on Proffitt and Michaud 1993
            kmu = 2*Bcoeff*pow(T,2.5d0)/(sqrt5*rho*chemZ*chemZ)/ &
               (XH1*sqrt(acx)*ccx + (1-XH1)*sqrt(acy)*ccy)

         else
            ! Log Lambda for H-He mixture (equation 10)
            loglambdah = -19.26d0 - log(2d0) - 0.5d0*log(rho) + &
               1.5d0*log(T) - 0.5d0*log(1d0 + 0.5d0*(1+XH1)) 
            ! Calculation of C_ij coeffs (equation 12)
            ccy = log(exp(1.2d0*loglambdah)+1d0)/1.2d0
            ! My formula (see notes) based on Proffitt and Michaud 1993
            kmu = (Bcoeff*pow(T,2.5d0)/(rho*ccy))*(3+XH1)/((1+XH1)*(3+5*XH1)*(0.7d0+0.3d0*XH1))
         
         endif
         ! write(57,*) kt,kmu,vis,chemZ

      end subroutine get_diff_coeffs


      real(dp) function numu(R0,r_th,prandtl,diffratio)
         !Function calculates Nu_mu from input parameters, following Brown et al. 2013.
         !Written by P. Garaud (2013). Please email pgaraud@ucsc.edu for troubleshooting. 

         real(dp), intent(in) :: R0,r_th,prandtl,diffratio
         real(dp) :: maxl2,maxl,lambdamax
         real(dp) :: myvars(2)
         integer :: ierr, iter, max_iters

         ! Initialize guess using estimates from Brown et al. 2013
         call analytical_estimate_th(maxl,lambdamax,r_th,prandtl,diffratio)
               
         myvars(1) = maxl
         myvars(2) = lambdamax

        !Call Newton relaxation algorithm
        call NR(myvars,prandtl,diffratio,R0,ierr)
       
        !If the growth rate is negative, then try another set of parameters as first guess.  
        !Repeat as many times as necessary until convergence is obtained.
        iter = 1
        max_iters = 200
        do while(iter<=max_iters .and. ((myvars(2)<0).or.(ierr /= 0))) 
           !write(*,*) 'Alternative', r_th,prandtl,diffratio,iter
           !Reset guess values
           myvars(1) = maxl
           myvars(2) = lambdamax
           !Call relaxation for slightly different Pr, tau, R0.
           call NR(myvars,prandtl*(1d0+iter*1.d-2),diffratio,R0/(1d0+iter*1.d-2),ierr)
           !If it converged this time, call NR for the real parameters.
           if(ierr.eq.0) call NR(myvars,prandtl,diffratio,R0,ierr)
           !write(*,*) prandtl,diffratio,R0,myvars(1),myvars(2),ierr
           !Otherwise, increase counter and try again.
           iter = iter + 1            
        enddo
        
        if((myvars(2)<0).or.(ierr /= 0)) then
           write(*,*) "WARNING: thermohaline Newton relaxation failed to converge, falling back to estimate"
           maxl2 = maxl*maxl
        else ! NR succeeded, so use results in myvars
           !Plug solution into "l^2" and lambda.
           maxl2 = myvars(1)*myvars(1)
           lambdamax = myvars(2)
           !write(*,*) prandtl,diffratio,r_th,maxl2,lambdamax
        end if

        !Calculate Nu_mu using Formula (33) from Brown et al, with C = 7.
        numu = 1.d0 + 49.d0*lambdamax*lambdamax/(diffratio*maxl2*(lambdamax+diffratio*maxl2))

         return
      end function numu 


      subroutine thermohaline_rhs(myx,myf,myj,prandtl,diffratio,R0)
         ! This routine is needed for the NR solver.
         ! Inputs the two following equations for lambda and maxl2:
         ! lambda^3 + a_2 lambda^2 + a_1 lambda + a_0 = 0 (eq. 19 of Brown et al.)
         ! b_2 lambda^2 + b_1 lambda + b_0 = 0 (eq. 20 of Brown et al.)
         ! Inputs f, the equations, and j, their jacobian.
         ! Written by P. Garaud (2013). Please email pgaraud@ucsc.edu for troubleshooting. 

         real(dp), intent(in) :: myx(2),  prandtl, diffratio, R0
         real(dp), intent(out) :: myf(2), myj(2,2)
         real(dp) :: a_2,a_1,a_0,b_2,b_1,b_0,myterm,myx1_2,myx1_3,myx1_4

         !This inputs the coefficients.
         b_2 = 1d0+prandtl+diffratio
         myx1_2 = myx(1)*myx(1)
         myx1_3 = myx1_2*myx(1)
         myx1_4 = myx1_3*myx(1)
         a_2 = myx1_2*b_2
         myterm = diffratio*prandtl+prandtl+diffratio
         b_1 = 2*myx1_2*myterm
         a_1 = myx1_4*myterm + prandtl*(1. - (1d0/R0))
         b_0 = 3.d0*myx1_4*diffratio*prandtl + prandtl*(diffratio - (1d0/R0))
         a_0 = myx1_4*myx1_2*diffratio*prandtl + myx1_2*prandtl*(diffratio - (1d0/R0))

!         write(*,*) a_2,a_1,a_0,b_2,b_1,b_0

         !These are equations 19 and 20
         myf(1) = ((myx(2) + a_2)*myx(2) + a_1)*myx(2) + a_0
         myf(2) = b_2*myx(2)*myx(2) + b_1*myx(2) + b_0
      
         !These are their Jacobians for the NR relaxation.
         myj(1,1) = 2*myx(1)*b_2*myx(2)*myx(2) + &
            4*myx1_3*myterm*myx(2) + 6*myx1_4*myx(1)*diffratio*prandtl   &
              + 2*myx(1)*prandtl*(diffratio - (1d0/R0))
         myj(1,2) = 3*myx(2)*myx(2) + 2*a_2*myx(2) + a_1
         myj(2,1) = 4*myx(1)*myterm*myx(2) + 12.d0*myx1_3*diffratio*prandtl
         myj(2,2) = 2*b_2*myx(2) + b_1

         return
      end subroutine thermohaline_rhs               


      subroutine analytical_estimate_th(maxl,lambdamax,r_th,prandtl,diffratio)
         !Inputs analytical estimates for l and lambda from Brown et al. 2013.

         real(dp) :: prandtl, diffratio, maxl, lambdamax, r_th, phi, maxl4, maxl6
      
         phi = diffratio/prandtl

         if(r_th .lt. 0.5d0) then
            if(r_th .gt. prandtl) then
               maxl = pow((1.d0/(1.d0+phi)) - 2.d0*dsqrt(r_th*phi)/pow(1d0+phi,2.5d0),0.25d0)   
                  ! Equation (B14)
               maxl4 = maxl*maxl*maxl*maxl
               maxl6 = maxl4*maxl*maxl
               lambdamax = 2*prandtl*phi*maxl6/(1d0-(1d0+phi)*maxl4)    ! Equation (B11)
            else
               maxl = dsqrt(dsqrt(1d0/(1d0+phi)) - dsqrt(prandtl)*(1d0+phi/((1d0+phi)*(1d0+phi))))  
                  ! Equation (B5)
               lambdamax = dsqrt(prandtl) - prandtl*dsqrt(1d0+phi)   !Equation (B5)
            endif
         else
            maxl = pow((one_third)*(1d0-r_th) + (1d0-r_th)*(1d0-r_th)*(5d0-4d0*phi)/27d0,0.25d0)
               ! Equation (B19) carried to next order (doesn't work well otherwise)
            maxl4 = maxl*maxl*maxl*maxl
            maxl6 = maxl4*maxl*maxl
            lambdamax = 2d0*prandtl*phi*maxl6/(1d0-(1d0+phi)*maxl4) ! Equation (B11)
         endif
         if(lambdamax<0) then   ! shouldn't be needed, but just as precaution
            maxl = 0.5d0
            lambdamax = 0.5d0
         endif
      
         return
      end subroutine analytical_estimate_th


      subroutine NR(xrk,prandtl,diffratio,R0,ierr)
         ! Newton Relaxation routine used to solve cubic & quadratic in thermohaline case.
         ! Written by P. Garaud (2013). Please email pgaraud@ucsc.edu for troubleshooting. 

         real(dp), parameter :: acy = 1.d-13 ! accuracy of NR solution.
         integer, parameter :: niter = 20  ! max number of iterations allowed before giving up.
         integer, parameter :: &  !array dimension input parameters for dgesvx
               n = 2, &
               nrhs = 1, &
               lda = n, &
               ldaf = n, &
               ldb = n, &
               ldx = n

         integer :: iter,ierr
         real(dp) :: xrk(2), f(2) ! Functions f 
         real(dp) :: j(2,2) ! Jacobian
         real(dp) :: err,errold ! Error at each iteration
         real(dp) :: x1_sav,x2_sav
         real(dp) :: prandtl, diffratio, R0
         real(dp) :: A(lda,n), AF(ldaf,n), R(n), C(n), B(ldb,nrhs), X(ldx,nrhs), &
               rcond, ferr(nrhs), berr(nrhs), work(4*n)
         character :: fact, trans, equed
         integer :: ipiv(n), iwork(n)

         include 'formats'

         !Initialize flags and other counters.
         ierr = 0
         iter = 0
         err = 1d99
         errold = 1d99
         !Save input guess
         x1_sav = xrk(1)
         x2_sav = xrk(2)

         !While error is too large .and. decreasing, iterate.
         do while ((err.gt.acy).and.(ierr.eq.0).and.(iter.lt.niter))
            call thermohaline_rhs(xrk,f,j,prandtl,diffratio,R0)    
         
            fact = 'E'
            trans = 'N'
            equed = ''
            
            A  = j
            B(1,1) = f(1)
            B(2,1) = f(2)

            call dgesvx( fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, &
               equed, r, c, B, ldb, x, ldx, rcond, ferr, berr, &
               work, iwork, ierr )

            if (ierr /= 0) then
               !write(*,*) 'dgesvx failed in thermohaline routine', iter
               !write(*,2) j(1,1),j(1,2)
               !write(*,2) j(2,1),j(2,2)
            else
               iter = iter + 1
               f(1) = X(1,1)
               f(2) = X(2,1)
               err = dsqrt(f(1)*f(1)+f(2)*f(2)) ! Calculate the new error
               ! If, after a while, the error is still not decreasing, give up and exit NR.
               ! Otherwise, continue.
               if((iter.gt.5).and.(err.gt.errold)) then              
                  ! Write(*,2) 'Error not decreasing at iter', iter, err, errold
                  ierr = 1
                  ! Reset xs and exit loop.
                  xrk(1) = x1_sav
                  xrk(2) = x2_sav                   
               else
                  xrk = xrk - f ! The solution is now in f, so update x 
                  errold = err
               endif
            endif
         enddo
         
         if(err<=acy) then
            ierr = 0
         else
            ! write(*,2) 'NR failed to converge', err, iter
            ierr = 1
         end if

         return
      end subroutine NR


      end module mlt_get_results
