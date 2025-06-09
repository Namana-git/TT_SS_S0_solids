module sort
     
   use global_variables
   use compute_lambda_and_xi
   implicit none

  contains
     
  subroutine sort_energy()
     integer:: I_ex,J_ex
     integer :: c,x,y,c2,z,i
     integer :: I_max,I_temp_s,I_temp_t
     double precision,allocatable :: temp_energy_arr(:),energy_arr(:)
     integer,allocatable :: temp_spin_ind(:)
     integer,allocatable :: temp_I_arr(:),temp_J_arr(:)
     double precision :: E_IJ
     
     sys%biex = sys%nex*(sys%nex+1)      ! No of unique I,J pairs
     allocate(temp_energy_arr(sys%biex))  ! sum of energies of I,J excitons 
     allocate(temp_I_arr(sys%biex))     !stores the I index 
     allocate(temp_J_arr(sys%biex))     !stores the J index
     allocate(temp_spin_ind(sys%biex))  ! stores whether it is tt or ss 
     
   
     !determine the nex_max and ncutoff
     c = 0
     c2 = 0
      do i = 1,2
        do I_ex =1,sys%nex
           do J_ex= 1,I_ex
              c = c+1
              if (i==1)then
               temp_energy_arr(c) = exciton_sys%eigenvalues_t(I_ex) + exciton_sys%eigenvalues_t(J_ex)
              end if
              if(i==2)then
                  temp_energy_arr(c) = exciton_sys%eigenvalues_s(I_ex) + exciton_sys%eigenvalues_s(J_ex)
              end if
              temp_I_arr(c) = I_ex 
              temp_J_arr(c) = J_ex 
              temp_spin_ind(c) = i
              !print*, "temp_energy_arr", temp_energy_arr
              if(temp_energy_arr(c) .le. 0.36997) then
                 c2 = c2 + 1
                 if(i == 1) then
                    !c2 = c2 + 1
                    I_temp_t = I_ex
                  end if
                 if(i == 2) then
                    !c3 = c3 + 1
                    I_temp_s = I_ex
                 end if

                !print *,temp_energy_arr(c),temp_I_arr(c),temp_J_arr(c)
               !print*,energy_arr(c),exciton_sys%I_arr(c),exciton_sys%J_arr(c)
              end if
           end do
        end do
      end do
  
     sys%ncutoff = c2
     !sys%ncutoff_s = c3
     sys%nex_max_t = I_temp_t
     sys%nex_max_s = I_temp_s
     !print*,sys%ncutoff,sys%nex_max
     print *,sys%ncutoff,sys%nex_max_s,sys%nex_max_t
      
     !do i = 1,2
      do x = 1,sys%biex
        do y = x+1,sys%biex
          if (temp_energy_arr(x) > temp_energy_arr(y)) then
                call swap(temp_energy_arr(x), temp_energy_arr(y))
                call swap_idx(temp_I_arr(x), temp_I_arr(y))
                call swap_idx(temp_J_arr(x), temp_J_arr(y))
                call swap_idx(temp_spin_ind(x), temp_spin_ind(y))
          end if 
        end do
      end do
      print*,"hello"
     allocate(exciton_sys%I_arr(sys%ncutoff))
     allocate(exciton_sys%J_arr(sys%ncutoff))
     allocate(exciton_sys%spin_ind(sys%ncutoff))

     exciton_sys%I_arr = temp_I_arr(1:sys%ncutoff)
     exciton_sys%J_arr = temp_J_arr(1:sys%ncutoff)
     exciton_sys%spin_ind = temp_spin_ind(1:sys%ncutoff)

     print*,"temp_arr:",temp_energy_arr
     print*,"exciton_sys%I_arr =", exciton_sys%I_arr
     print*,"exciton_sys%J_arr =", exciton_sys%J_arr
     !print*,"exciton_sys%spin_ind =", exciton_sys%spin_ind
     

     !print*,(exciton_sys%eigenvalues(exciton_sys%I_arr(sys%ncutoff)) + exciton_sys%eigenvalues(exciton_sys%J_arr(sys%ncutoff))),exciton_sys%I_arr(sys%ncutoff),exciton_sys%J_arr(sys%ncutoff)
    ! print*
     !print*,"I",exciton_sys%I_arr
     !print*,"J",exciton_sys%J_arr
    !print*,"spin_ind",exciton_sys%spin_ind
    print*,"exciton_sys%I_arr",exciton_sys%I_arr
    print*,"exciton_sys%J_arr",exciton_sys%J_arr
    deallocate(temp_energy_arr)
    deallocate(temp_I_arr)
    deallocate(temp_J_arr)


   end subroutine sort_energy


   subroutine swap(a, b)
      double precision, intent(inout) :: a, b
      double precision :: temp
      temp = a
      a = b
      b = temp
  end subroutine swap

  subroutine swap_idx(a, b)
   integer, intent(inout) :: a, b
   integer :: temp
   temp = a
   a = b
   b = temp
  end subroutine swap_idx






end module sort

