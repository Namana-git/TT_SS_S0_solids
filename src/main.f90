program mol_shiva
    use compute_AB
    use global_variables
    use input_read
    use kpoints
    !use sort
    use solve_generalized_eigenvalue_problem
    !use solve_generalized_eigenvalue_problem
    !use sort
    implicit none

    double precision :: start_time,end_time
    double precision :: kpt(3), kpt1(3), kpt2(3), kpt3(3), kpt4(3), kpt5(3), kpt6(3), kpt7(3)
    double precision :: Qkpt(3), Qpt1(3), Qpt2(3), Qpt3(3)
    integer :: index 
    integer :: error
    sys%nQ = 4
    sys%nc = 1
    sys%nv = 1
    sys%neb = 4*1*1
    call read_kpoints()
    call read_Qpoints()
    call read_triplet_exciton_eigenvalues()
   
    call read_triplet_exciton_eigenvectors()

   

    !Qkpt = [0.00,0.00,0.0]
    !Qpt1 = [0.00,0.00,-0.25]
    !Qpt2 = [0.00,0.00,-0.5]
    !Qpt3 = [0.00,0.00,0.25]

   ! call Qpoint_to_index(Qkpt,index)
   ! print*, "Index of Qkpt ", index,sys%Qpts(1,index),sys%Qpts(2,index),sys%Qpts(3,index)
   ! call Qpoint_to_index(Qpt1,index)
   ! print*, "Index of Qpt1 ", index,sys%Qpts(1,index),sys%Qpts(2,index),sys%Qpts(3,index)
    !call Qpoint_to_index(Qpt2,index)    
    !print*, "Index of Qpt2 ", index,sys%Qpts(1,index),sys%Qpts(2,index),sys%Qpts(3,index)
    !call Qpoint_to_index(Qpt3,index)
    !print*, "Index of Qpt3 ", index,sys%Qpts(1,index),sys%Qpts(2,index),sys%Qpts(3,index)
    !kpt4 = [0.00,0.00,0.5]
    !kpt5 = [0.00,0.00,-0.5]
    !kpt6 = [0.00,0.00,-0.25]
    !kpt7 = [0.00,0.00,-0.75]
  

 !   print*,"1",anint(1.0)
 !   print*,"0.5",anint(0.5)
 !   print*,"0.25",anint(0.25)
 !   print*,"0.75",anint(0.75)
 !   print*,"0.0",anint(0.0)
 !    print*,"-0.5",anint(-0.5)
 !   print*,"-0.25",anint(-0.25)
 !   print*,"-0.75",anint(-0.75)
 !   print*,"0",anint(0.0)


  !  call kpoint_to_index(kpt,index)
  !  print*, "Index of kpt  ", index,sys%kpts(1,index),sys%kpts(2,index),sys%kpts(3,index)
  !  call kpoint_to_index(kpt1,index)
  !  print*, "Index of kpt1 ", index,sys%kpts(1,index),sys%kpts(2,index),sys%kpts(3,index)
  !  call kpoint_to_index(kpt2,index)
  !  print*, "Index of kpt2 ", index,sys%kpts(1,index),sys%kpts(2,index),sys%kpts(3,index)
  !  call kpoint_to_index(kpt3,index)        
  ! print*, "Index of kpt3 ", index,sys%kpts(1,index),sys%kpts(2,index),sys%kpts(3,index)
    !call kpoint_to_index(kpt4,index)
    !print*, "Index of kpt4 ", index,sys%kpts(1,index),sys%kpts(2,index),sys%kpts(3,index)
    !call kpoint_to_index(kpt5,index)
  !  print*, "Index of kpt5", index,sys%kpts(1,index),sys%kpts(2,index),sys%kpts(3,index)
  !  call kpoint_to_index(kpt6,index)
  !  print*, "Index of kpt6 ", index,sys%kpts(1,index),sys%kpts(2,index),sys%kpts(3,index)
  !  call kpoint_to_index(kpt7,index)
  !  print*, "Index of kpt7 ", index,sys%kpts(1,index),sys%kpts(2,index),sys%kpts(3,index)
  !!  !start time
  !!  call cpu_time(start_time)
  !!  !initialize the no of excitons included
  !!  sys%nex = 4
  !  sys%nex_max = 4 
!
  !  !initialize the no of conduction anv valence states
  !  call read_input()
  !  ! read eigenvalues and eigenvect
  !  !print*,"hello"
  !  call read_singlet_exciton_eigenvalues()
  !  !print*,"hello1"
  !  call read_singlet_exciton_eigenvectors()
  !  call read_triplet_exciton_eigenvalues()
  !  call read_triplet_exciton_eigenvectors()
!
    !print*,"sys%nex",sys%nex
    !print*,"exciton_sys%eigenvalues_s",exciton_sys%eigenvalues_s
    !print*,"exciton_sys%eigenvalues_t",exciton_sys%eigenvalues_t
    !print*,"exciton_sys%eigenvectors_s",exciton_sys%eigenvectors_s
    !print*,"exciton_sys%eigenvectors_t",exciton_sys%eigenvectors_t

   ! print*,OMP_GET_NUM_THREADS()
    !sys%nex = 2
    !ham%l = sys%nex*(sys%nex + 1)/2
    !ham%l = sys%nex*sys%nex
   ! sys%northo = 600  
    !ham%l = 4
    !sys%nv = 20
    !sys%nc = 40
    !call sort_energy()
    
   ! call sort_energy()
     !print*,"hey"
    !compute A and B Hamiltonian from lambda and xi matrices
    !call compute_all_lambda()
    call compute_lambda()
  !  call compute_xi()
  !  call construct_AB_SS()
  !  call construct_AB_TT()q
  !  call construct_AB_TTSS()
  !  call construct_AB_SSTT()
    !call construct_A_SS()
    !call construct_A_TT()
    !call construct_A_TTSS()
    !call construct_A_SSTT()
    !call compute_xi()
    !call compute_xi_in()
    !call compute_TT_matrix()
   call solve_eigenvalue_overlap()
  ! call solve_eigenvalue_O()
    !call compute_H_ortho1()
   !call solve_eigenvalue_H()
   !print*,end_time-start_time
end program  

