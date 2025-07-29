module global_variables
    implicit none
    

    type :: input
    
    integer :: nc,nv,nex,nb,northo,ncutoff,nex_max,biex,num_ortho,nex_max1,nex_max_s,nex_max_t,lin_ind  !no of conduction, valence, exciton states and c+v total no of states 
    integer :: nk,nQ,neb
    double precision :: vol
    double precision, allocatable :: kpts(:,:),Qpts(:,:)
    end type input
    type(input) :: sys

    

    type :: exciton
        double precision, allocatable :: eigenvalues_s(:,:) !I,Q
        complex(kind=8), allocatable :: eigenvectors_s(:,:,:,:,:,:) !s,v,c,k,I,Q
        double precision, allocatable :: eigenvalues_t(:,:)!I,Q
        complex(kind=8), allocatable :: eigenvectors_t(:,:,:,:,:,:) !s,v,c,k,I,Q
        integer , allocatable :: I_arr(:),J_arr(:),spin_ind(:)
    end type exciton
    type(exciton) :: exciton_sys



    
    type :: Hamiltonian
          integer :: l
          double precision, allocatable :: A_tttt(:,:,:,:),A_ssss(:,:,:,:),A_sstt(:,:,:,:),A_ttss(:,:,:,:)
          double precision, allocatable :: B_tttt(:,:,:,:),B_ssss(:,:,:,:),B_sstt(:,:,:,:),B_ttss(:,:,:,:)
          complex(kind=8), allocatable :: lambda_tttt(:,:,:,:,:,:),lambda_ssss(:,:,:,:,:,:)
          complex(kind=8), allocatable :: lambda_sstt(:,:,:,:,:,:),lambda_ttss(:,:,:,:,:,:)  
          complex(kind=8), allocatable :: xi_d_tttt(:,:,:,:,:,:),xi_d_ssss(:,:,:,:,:,:)
          complex(kind=8), allocatable :: xi_x_tttt(:,:,:,:,:,:),xi_x_ssss(:,:,:,:,:,:)
          complex(kind=8), allocatable :: xi_d_in_tttt(:,:,:,:,:,:),xi_d_in_ssss(:,:,:,:,:,:)
          complex(kind=8), allocatable :: xi_x_in_tttt(:,:,:,:,:,:),xi_x_in_ssss(:,:,:,:,:,:)
          complex(kind=8), allocatable :: xi_d_in_ttss(:,:,:,:,:,:),xi_d_in_sstt(:,:,:,:,:,:)
          complex(kind=8), allocatable :: xi_x_in_ttss(:,:,:,:,:,:),xi_x_in_sstt(:,:,:,:,:,:)
          complex(kind=8), allocatable :: xi_x_ttss(:,:,:,:,:,:),xi_x_sstt(:,:,:,:,:,:)
          complex(kind=8),allocatable :: H(:,:),O(:,:),c(:,:),A(:,:)
    end  type Hamiltonian

    


    type(Hamiltonian) :: ham



          
end module global_variables
