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
    complex(kind=8), dimension(:,:,:,:,:,:,:),allocatable :: bsemat,bsemat_ee
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

    !sys%nv = 2  !debug statement

 

   ! call sort_energy()
     !print*,"hey"
    !compute A and B Hamiltonian from lambda and xi matrices
    !call compute_all_lambda()
    call compute_lambda()
    call compute_xi()
    call debug_xi()
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
    !call solve_eigenvalue_overlap()
   !call solve_eigenvalue_O()
    !call compute_H_ortho1()
   !call solve_eigenvalue_H()
   !print*,end_time-start_time
end program  

