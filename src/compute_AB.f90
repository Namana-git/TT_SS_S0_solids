module compute_AB
    !use stdlib_linalg, only: outer_product
    use global_variables
    use compute_lambda_and_xi
    use input_read
    use, intrinsic :: iso_fortran_env
    !use blas_interface
    implicit none



    ! calculates A and B from \lambda and \xi matrices
    !(E_{M} + E_{N} -\epsilon_{m})c_{MN} + \sum_{IJ} \xi^{tot} (M,N,I,J) = 0
    !\xi^{tot}(M,N,I,J) = \xi(M,N,I,J) + \xi^{in}(M,N,I,J) + (E_{I}+E_{J} - \epsilon_{m})\lambda (M,N,I,J)
    !A(x(MN),y(IJ)) = (E_{M} + E_{N})\delta_{MN,IJ) +  \xi(M,N,I,J) + \xi^{in}(M,N,I,J) + (E_{I}+E_{J})\lambda (M,N,I,J) 
    !B(x(MN,Y(IJ))) = \delta_{MN,IJ) + \lambda (M,N,I,J)     
contains

   subroutine compute_lambda()

       implicit none
       integer :: I,J,M,N
       integer :: iQ_r,iQ_l
       integer :: c1,c2
       complex(kind=8) :: lambda1,lambda2,lambda3,lambda4
       !double precision,allocatable,dimension(:,:,:,:) :: bsemat
        print*,"hey4"    
    
       lambda1 = cmplx(0.0,0.0) 
       lambda2 = cmplx(0.0,0.0)
       lambda3 = cmplx(0.0,0.0)
       lambda4 = cmplx(0.0,0.0)
       sys%nb = sys%nc+sys%nv
      ! allocate(ham%lambda_ee_ssss(sys%nex_max,sys%nex_max,sys%nex_max,sys%nex_max))
       allocate(ham%lambda_tttt(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
       !allocate(ham%lambda_ssss(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
       !allocate(ham%lambda_ttss(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
       !allocate(ham%lambda_sstt(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
      ! allocate(ham%lambda_ttss(sys%nex_max,sys%nex_max,sys%nex_max,sys%nex_max))
       !allocate(ham%lambda_hh_ssss(sys%nex_max,sys%nex_max,sys%nex_max,sys%nex_max))
      
       !allocate(ham%xi(sys%nex,sys%nex,sys%nex,sys%nex))
       ham%lambda_tttt = cmplx(0.0,0.0)
      ! ham%lambda_ssss = cmplx(0.0,0.0)
      ! ham%lambda_ttss =  cmplx(0.0,0.0)
      ! ham%lambda_sstt = cmplx(0.0,0.0)
       !ham%lambda_ttss = 0.0
       !ham%xi = 0.0
        print*,"hey5"
       !allocate(bsemat(sys%nb,sys%nb,sys%nb,sys%nb))
      ! call read_bsemat(bsemat)
       print*,"a."
      
        c1 = 0 
        c2 = 0
      ! open(UNIT=10, FILE="lambda.dat", STATUS="REPLACE")
      do iQ_r = 1,sys%nQ
         do iQ_l = 1,sys%nQ
            do M = 1,sys%neb
              do N = 1,sys%neb
                 do I = 1,sys%neb
                    do J = 1,sys%neb
                    
                     call compute_lambda_tttt(I,J,M,N,iQ_r,iQ_l,lambda1)
                      !print*,"lambda",I,J,M,N,iQ_r,iQ_l,lambda1
                      ham%lambda_tttt(I,J,M,N,iQ_r,iQ_l) = lambda1
                      if(I==1 .and. J==2 .and. M==1 .and. N==4 .and. iQ_r==1 .and. iQ_l== 2)then 
                         print*,"yellow",lambda1
                        print*,"lambda_tttt",I,J,M,N,iQ_r,iQ_l,ham%lambda_tttt(I,J,M,N,iQ_r,iQ_l)
                      end if 
                  !  call compute_lambda_ssss(I,J,M,N,Q_r,Q_l,lambda2)
                  !   !print*,"lambda",I,J,M,N,lambda2
                  !    ham%lambda_ssss(I,J,M,N) = lambda2
                  !  call compute_lambda_ttss(I,J,M,N,Q_r,Q_l,lambda3)
                  !   !print*,"lambda",I,J,M,N,lambda3
                  !    ham%lambda_ttss(I,J,M,N) = lambda3
!
                  !   call compute_lambda_sstt(I,J,M,N,Q_r,Q_l,lambda4)
                  !   !print*,"lambda",I,J,M,N,lambda4
                  !    ham%lambda_sstt(I,J,M,N) = lambda4
!
                  !    if(I==1 .and. J==1 .and. M==1 .and. N==1)then 
                  !       print*,"yellow",lambda1
                  !    end if 
                   
                      !print*,M,N,I,J,ham%lambda(I,J,M,N)
                  end do
              end do
            end do
          end do
        end do
      end do
        !ham%lambda_tttt = 0.0
        !ham%lambda_ssss = 0.0
      
       !close(10)
    
    end subroutine compute_lambda
    subroutine compute_xi()

        implicit none
        integer :: I,J,M,N
        integer :: iQ_r,iQ_l
        integer :: c1,c2
        complex(kind=8) :: xi_ee1_tttt,xi_ee2_tttt,xi_ee1_ssss,xi_ee2_ssss,xi_ee1_ttss,xi_ee2_ttss,xi_ee1_sstt,xi_ee2_sstt, &
                            xi_hh1_tttt,xi_hh2_tttt,xi_hh1_ssss,xi_hh2_ssss,xi_hh1_ttss,xi_hh2_ttss,xi_hh1_sstt,xi_hh2_sstt, &
                            xi_in_ehd1_tttt,xi_in_ehd2_tttt,xi_in_ehd1_ssss,xi_in_ehd2_ssss,xi_in_ehd1_sstt,xi_in_ehd2_sstt,xi_in_ehd1_ttss,xi_in_ehd2_ttss, &
                            xi_out_ehd1_tttt,xi_out_ehd2_tttt,xi_out_ehd1_ssss,xi_out_ehd2_ssss,xi_out_ehd1_sstt,xi_out_ehd2_sstt,xi_out_ehd1_ttss,xi_out_ehd2_ttss, &
                            xi_in_ehx1_tttt,xi_in_ehx2_tttt,xi_in_ehx1_ssss,xi_in_ehx2_ssss,xi_in_ehx1_sstt,xi_in_ehx2_sstt,xi_in_ehx1_ttss,xi_in_ehx2_ttss, &
                            xi_out_ehx1_tttt,xi_out_ehx2_tttt,xi_out_ehx1_ssss,xi_out_ehx2_ssss,xi_out_ehx1_sstt,xi_out_ehx2_sstt,xi_out_ehx1_ttss,xi_out_ehx2_ttss
        complex(kind=8),allocatable,dimension(:,:,:,:,:,:,:) :: bsemat_ee,bsemat_hh,bsemat_d,bsemat_x
     
       ! ham%lambda = 0.0 
        !lambda1 = 0.0    
        !allocate(ham%lambda(sys%nex,sys%nex,sys%nex,sys%nex))
        !allocate(ham%lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
        allocate(ham%xi_d_tttt(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        allocate(ham%xi_x_tttt(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        allocate(ham%xi_d_in_tttt(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        !allocate(ham%xi_d_ssss(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        !allocate(ham%xi_x_ssss(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        !allocate(ham%xi_d_in_ssss(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        !allocate(ham%xi_x_in_ssss(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        !allocate(ham%xi_d_in_sstt(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        !allocate(ham%xi_x_in_sstt(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        !allocate(ham%xi_x_sstt(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        !allocate(ham%xi_d_in_ttss(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        !allocate(ham%xi_x_in_ttss(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        !allocate(ham%xi_x_ttss(sys%neb,sys%neb,sys%neb,sys%neb,sys%nQ,sys%nQ))
        !ham%lambda = 0.0
        !ham%lambda_e = 0.0
        ham%xi_d_tttt = cmplx(0.0,0.0)
        ham%xi_x_tttt = cmplx(0.0,0.0)
        ham%xi_d_in_tttt = cmplx(0.0,0.0)
       ! ham%xi_d_ssss = 0.0
       ! ham%xi_x_ssss = 0.0
       ! ham%xi_d_in_ssss = 0.0
       ! ham%xi_x_in_ssss = 0.0
       ! ham%xi_x_in_sstt = 0.0
       ! ham%xi_x_sstt = 0.0
       ! ham%xi_d_in_sstt = 0.0
       ! ham%xi_d_in_ttss = 0.0
       ! ham%xi_x_in_ttss = 0.0
       ! ham%xi_x_ttss = 0.0
        !print*,"hello1"
        !
       ! print*,"hey6"

     
        !allocate(bsemat_d(sys%nb,sys%nb,sys%nb,sys%nb))
        !allocate(bsemat_x(sys%nb,sys%nb,sys%nb,sys%nb))
        !bsemat_d = 0.0
        !bsemat_x = 0.0
        !print*,"hello4"
         allocate(bsemat_ee(sys%nc,sys%nc,sys%nc,sys%nc,sys%nk,sys%nk,sys%nQ))
         bsemat_ee = cmplx(0.0,0.0)
         call read_bsemat_ee(bsemat_ee)
         allocate(bsemat_hh(sys%nv,sys%nv,sys%nv,sys%nv,sys%nk,sys%nk,sys%nQ))
         bsemat_hh = cmplx(0.0,0.0)
         call read_bsemat_hh(bsemat_hh)
         allocate(bsemat_d(sys%nv,sys%nv,sys%nc,sys%nc,sys%nk,sys%nk,sys%nQ))
         bsemat_d = cmplx(0.0,0.0)
         call read_bsemat_d(bsemat_d)
        !call read_bsemat_d(bsemat_d)
        !call read_bsemat_x(bsemat_x)
        print*,"hello5"
       
         c1 = 0 
         c2 = 0
        ! open(UNIT=10, FILE="lambda.dat", STATUS="REPLACE")
         do iQ_r = 1,sys%nQ
          do iQ_l = 1,sys%nQ
            do M = 1,sys%neb
              do N = 1,sys%neb
                 do I = 1,sys%neb
                    do J = 1,sys%neb
                       call compute_xi_ee1(I,J,M,N,iQ_r,iQ_l,bsemat_ee,xi_ee1_tttt,xi_ee1_ssss,xi_ee1_sstt,xi_ee1_ttss)
                       call compute_xi_ee2(I,J,M,N,iQ_r,iQ_l,bsemat_ee,xi_ee2_tttt,xi_ee2_ssss,xi_ee2_sstt,xi_ee2_ttss)
                       call compute_xi_hh1(I,J,M,N,iQ_r,iQ_l,bsemat_hh,xi_hh1_tttt,xi_hh1_ssss,xi_hh1_sstt,xi_hh1_ttss)
                       call compute_xi_hh2(I,J,M,N,iQ_r,iQ_l,bsemat_hh,xi_hh2_tttt,xi_hh2_ssss,xi_hh2_sstt,xi_hh2_ttss)
                       call compute_xi_in_ehd1(I,J,M,N,iQ_r,iQ_l,bsemat_d,xi_in_ehd1_tttt,xi_in_ehd1_ssss,xi_in_ehd1_sstt,xi_in_ehd1_ttss)
                       call compute_xi_in_ehd2(I,J,M,N,iQ_r,iQ_l,bsemat_d,xi_in_ehd2_tttt,xi_in_ehd2_ssss,xi_in_ehd2_sstt,xi_in_ehd2_ttss)
                       call compute_xi_out_ehd1(I,J,M,N,iQ_r,iQ_l,bsemat_d,xi_out_ehd1_tttt,xi_out_ehd1_ssss,xi_out_ehd1_sstt,xi_out_ehd1_ttss)
                       call compute_xi_out_ehd2(I,J,M,N,iQ_r,iQ_l,bsemat_d,xi_out_ehd2_tttt,xi_out_ehd2_ssss,xi_out_ehd2_sstt,xi_out_ehd2_ttss)
                    !call compute_xi_hh(I,J,M,N,bsemat_d,xi_hh1_tttt,xi_hh2_tttt,xi_hh1_ssss,xi_hh2_ssss,xi_hh1_sstt,xi_hh2_sstt,xi_hh1_ttss,xi_hh2_ttss)
                    !call compute_xi_in_ehd(I,J,M,N,bsemat_d,xi_in_ehd1_tttt,xi_in_ehd2_tttt,xi_in_ehd1_ssss,xi_in_ehd2_ssss,xi_in_ehd1_sstt,xi_in_ehd2_sstt,xi_in_ehd1_ttss,xi_in_ehd2_ttss)   
                    !call compute_xi_out_ehd(I,J,M,N,bsemat_d,xi_out_ehd1_tttt,xi_out_ehd2_tttt,xi_out_ehd1_ssss,xi_out_ehd2_ssss,xi_out_ehd1_sstt,xi_out_ehd2_sstt,xi_out_ehd1_ttss,xi_out_ehd2_ttss)
                    !call compute_xi_in_ehx(I,J,M,N,bsemat_x,xi_in_ehx1_tttt,xi_in_ehx2_tttt,xi_in_ehx1_ssss,xi_in_ehx2_ssss,xi_in_ehx1_sstt,xi_in_ehx2_sstt,xi_in_ehx1_ttss,xi_in_ehx2_ttss)
                    !call compute_xi_out_ehx(I,J,M,N,bsemat_x,xi_out_ehx1_tttt,xi_out_ehx2_tttt,xi_out_ehx1_ssss,xi_out_ehx2_ssss,xi_out_ehx1_sstt,xi_out_ehx2_sstt,xi_out_ehx1_ttss,xi_out_ehx2_ttss)
                    !print*,"hello7"
                    !ham%xi_d_ssss(I,J,M,N) = xi_ee1_ssss + xi_hh1_ssss  - xi_in_ehd1_ssss - xi_out_ehd1_ssss
                    !ham%xi_x_ssss(I,J,M,N) = xi_in_ehx1_ssss + xi_out_ehx1_ssss
                    !ham%xi_d_in_ssss(I,J,M,N) = xi_ee2_ssss + xi_hh2_ssss  - xi_in_ehd2_ssss - xi_out_ehd2_ssss
                    !ham%xi_x_in_ssss(I,J,M,N) = xi_in_ehx2_ssss + xi_out_ehx2_ssss
                     !print*,"hello8"
                    ham%xi_d_tttt(I,J,M,N,iQ_r,iQ_l) = xi_ee1_tttt + xi_hh1_tttt  - xi_in_ehd1_tttt - xi_out_ehd1_tttt
                    !ham%xi_x_tttt(I,J,M,N) = xi_in_ehx1_tttt + xi_out_ehx1_tttt
                    ham%xi_d_in_tttt(I,J,M,N,iQ_r,iQ_l) = xi_ee2_tttt + xi_hh2_tttt  - xi_in_ehd2_tttt - xi_out_ehd2_tttt
                    !print*,ham%xi_d_in_tttt(I,J,M,N,iQ_r,iQ_l),ham%xi_d_in_tttt(I,J,M,N,iQ_r,iQ_l)
                    !ham%xi_x_in_tttt(I,J,M,N) = xi_in_ehx2_tttt + xi_out_ehx2_tttt
                     !print*,"hello8"
                    !ham%xi_d_sstt(I,J,M,N) = xi_ee1_sstt + xi_hh1_sstt  - xi_in_ehd1_sstt - xi_out_ehd1_sstt
                    !ham%xi_x_sstt(I,J,M,N) = xi_in_ehx1_sstt + xi_out_ehx1_sstt
                    !ham%xi_d_in_sstt(I,J,M,N) = xi_ee2_sstt + xi_hh2_sstt  - xi_in_ehd2_sstt - xi_out_ehd2_sstt
                    !ham%xi_x_in_sstt(I,J,M,N) = xi_in_ehx2_sstt + xi_out_ehx2_sstt
                    !ham%xi_d_ttss(I,J,M,N) = xi_ee1_ttss + xi_hh1_ttss  - xi_in_ehd1_ttss - xi_out_ehd1_ttss
                    !ham%xi_x_ttss(I,J,M,N) = xi_in_ehx1_ttss + xi_out_ehx1_ttss
                    !ham%xi_d_in_ttss(I,J,M,N) = xi_ee2_ttss + xi_hh2_ttss  - xi_in_ehd2_ttss - xi_out_ehd2_ttss
                    !ham%xi_x_in_ttss(I,J,M,N) = xi_in_ehx2_ttss + xi_out_ehx2_ttss

                    !print*,I,J,M,N,ham%xi_x_ssss(I,J,M,N),ham%xi_x_tttt(I,J,M,N),ham%xi_x_sstt(I,J,M,N),ham%xi_x_ttss(I,J,M,N)
                    !print*,xi_d_in_tttt(I,J,M,N,iQ_r,iQ_l),xi_d_tttt(I,J,M,N,iQ_r,iQ_l)
                  end do
               end do
             end do
           end do
         end do
       end do
       print*,"hello6"
     
         !close(10)
     
     end subroutine compute_xi

      subroutine debug_xi()
          integer :: I,J,M,N,R,S
          integer :: iQ_r,iQ_l,iQ_m
          complex(kind=8) :: sum
         do iQ_r = 1,sys%nQ
           do iQ_l = 1,sys%nQ
            do M = 1,sys%neb
              do N = 1,sys%neb
                do I = 1,sys%neb
                  do J = 1,sys%neb
                     sum = cmplx(0.0,0.0)
                    do iQ_m = 1,sys%nQ
                      do R = 1,sys%neb
                        do S = 1,sys%neb
                          ! print*,"hey9"
                           sum = sum + ham%lambda_tttt(R,S,M,N,iQ_m,iQ_l)*ham%xi_d_tttt(I,J,R,S,iQ_r,iQ_m)
                           ! if(I==2 .and. J==1 .and. M==1 .and. N==2 .and. iQ_r==1 .and. iQ_l==1) then
                            !  print*,"RSQMN,lambda",R,S,iQ_m,M,n,ham%lambda_tttt(R,S,M,N,iQ_m,iQ_l)
                             ! print*,"IJRSQ,xi",I,J,R,S,iQ_m,ham%xi_d_tttt(I,J,R,S,iQ_r,iQ_m)
                        !print*,"debug xi_d_in_tttt",I,J,M,N,iQ_r,iQ_l,ham%xi_d_in_tttt(I,J,M,N,iQ_r,iQ_l),sum
                     ! print*,"debug xi_d_in_tttt",I,J,M,N,iQ_r,iQ_l,ham%xi_d_in_tttt(I,J,M,N,iQ_r,iQ_l),sum
                    !end if
                    !print*,"debug xi_d_in_tttt",ham%xi_d_tttt(I,J,M,N,iQ_r,iQ_l),ham%xi_d_in_tttt(I,J,M,N,iQ_r,iQ_l)
                    !print*,"hey",I,J,M,N,iQ_r,iQ_l,abs(ham%xi_d_in_tttt(I,J,M,N,iQ_r,iQ_l)-sum)
                     !end if
                        
                        end do
                        end do
                    end do
                    !print*,"debug xi_d_in_tttt",I,J,M,N,iQ_r,iQ_l,abs(ham%xi_d_in_tttt(I,J,M,N,iQ_r,iQ_l)-sum)
                    if(abs(ham%xi_d_in_tttt(I,J,M,N,iQ_r,iQ_l)-sum) > 1.0d-8) then
                      !print*,"yellow,yellow"
                     print*,"debug",I,J,M,N,iQ_r,iQ_l,ham%xi_d_tttt(I,J,M,N,iQ_r,iQ_l),ham%xi_d_in_tttt(I,J,M,N,iQ_r,iQ_l),sum
                    end if
                  end do
                end do  
               end do
            end do
           end do
         end do
                           
                       

      end subroutine debug_xi

!
   !  subroutine construct_A_TT()
   !  ! A_{MN,IJ} = (E_{M} + E_{N})\delta_{MN,IJ) +  xi_d_tttt(IJMN) + 3/2 * xi_x_tttt(IJMN) &
   !  ! + 1/2* (E_{I}+E_{J})*(lambda_h(IJMN)+lambda_e(IJMN)) &
   !  ! +1/2*(xi_d_in_tttt(IJMN) + xi_d_im_tttt(IJNM)
   !   
!
   !   implicit none
   !!   integer :: I,J,M,N,I_ex,J_ex,M_ex,N_ex
   !!   integer :: c1,c2
!!
   !!   allocate(ham%A_tttt(sys%nex*sys%nex,sys%nex*sys%nex))
   !!   allocate(ham%B_tttt(sys%nex*sys%nex,sys%nex*sys%nex))
   !!   ham%A_tttt = 0.0
   !!   ham%B_tttt = 0.0
   !!   c1 =0
   !!   c2 =0
   !!   do I = 1,sys%nex
   !!      do J = 1,sys%nex
   !!       c2 = c2 + 1
   !!       c1 = 0
   !!         do M = 1,sys%nex
   !!            do N = 1,sys%nex
   !!             c1 = c1 + 1
!!
   !!               if (I==M .and. J==N) then
   !!                  ham%A_tttt(c1,c2) = ham%A_tttt(c1,c2) + (exciton_sys%eigenvalues_t(I) + exciton_sys%eigenvalues_t(J))
   !!                  ham%B_tttt(c1,c2) = ham%B_tttt(c1,c2) + 1
   !!               end if
!!
   !!               if(I==N .and. J==M)then
   !!                  ham%A_tttt(c1,c2) = ham%A_tttt(c1,c2) + (exciton_sys%eigenvalues_t(I) + exciton_sys%eigenvalues_t(J))
   !!                  ham%B_tttt(c1,c2) = ham%B_tttt(c1,c2) + 1 
   !!               end if 
   !!                  ham%A_tttt(c1,c2) = ham%A_tttt(c1,c2) + (ham%xi_d_tttt(I,J,M,N) + ham%xi_d_tttt(I,J,N,M)) &
   !!                                  + 1.5 * (ham%xi_x_tttt(I,J,M,N)+ham%xi_x_tttt(I,J,N,M)) &
   !!                                  + 0.5 * (exciton_sys%eigenvalues_t(I) + exciton_sys%eigenvalues_t(J)) * (ham%lambda_tttt(I,J,M,N) + ham%lambda_tttt(I,J,N,M)) &
   !!                                  + 0.5 * (ham%xi_d_in_tttt(I,J,M,N) + ham%xi_d_in_tttt(I,J,N,M))
   !!               
   !!                  ham%B_tttt(c1,c2) = ham%B_tttt(c1,c2) + 0.5*(ham%lambda_tttt(I,J,M,N) +ham%lambda_tttt(I,J,N,M))
   !!               print*,I,J,M,N,ham%A_tttt(c1,c2)
   !!             end do
   !!         end do
   !!       end do
   !!     end do
   !!  end subroutine construct_A_TT
!!
   !!  subroutine construct_A_SS()
   !!   !A_{MN,IJ} =  (E_{I} + E_{J})*(delta_{I,M}* delta_{JN} &
   !!   !           + xi_d_ssss(I,J,M,N)
   !!   !           + 0.5*xi_x_ssss(I,J,M,N)
   !!   !           - 0.5*(E_{I} + E_{J})*(lambda(I,J,M,N) + lambda(I,J,N,M))
   !!   !           - 0.5*(xi_d_in_ssss(I,J,M,N) + xi_d_in_ssss(I,J,N,M))
   !!   !           - (xi_x_in_ssss(I,J,M,N) + xi_x_in_ssss(I,J,N,M))
!!
   !!   implicit none
   !!   integer :: I,J,M,N
   !!   integer :: c1,c2
!!
   !!   allocate(ham%A_ssss(sys%nex*sys%nex,sys%nex*sys%nex))
   !!   allocate(ham%B_ssss(sys%nex*sys%nex,sys%nex*sys%nex))
   !!   ham%A_ssss = 0.0
   !!   ham%B_ssss = 0.0
!!
   !!   c1 =0
   !!   c2 =0
   !!   do I = 1,sys%nex
   !!      do J = 1,sys%nex
   !!       c2 = c2 + 1
   !!       c1 = 0
   !!         do M = 1,sys%nex
   !!             do N = 1,sys%nex
   !!               c1 = c1 + 1
   !!               if (I==M .and. J==N) then
   !!                  ham%A_ssss(c1,c2) = ham%A_ssss(c1,c2) + (exciton_sys%eigenvalues_s(I) + exciton_sys%eigenvalues_s(J))
   !!                  ham%B_ssss(c1,c2) = ham%B_ssss(c1,c2) + 1
   !!               end if
   !!               if (I==N .and. J==M) then
   !!                  ham%A_ssss(c1,c2) = ham%A_ssss(c1,c2) + (exciton_sys%eigenvalues_s(I) + exciton_sys%eigenvalues_s(J))
   !!                  ham%B_ssss(c1,c2) = ham%B_ssss(c1,c2) + 1
   !!               end if
   !!               ham%A_ssss(c1,c2) = ham%A_ssss(c1,c2) + (ham%xi_d_ssss(I,J,M,N)+ham%xi_d_ssss(I,J,N,M)) &
   !!                                  + 0.5 * (ham%xi_x_ssss(I,J,M,N)+ham%xi_x_ssss(I,J,N,M)) &
   !!                                  - 0.5 * (exciton_sys%eigenvalues_s(I) + exciton_sys%eigenvalues_s(J)) * (ham%lambda_ssss(I,J,M,N) + ham%lambda_ssss(I,J,N,M)) &
   !!                                  - 0.5 * (ham%xi_d_in_ssss(I,J,M,N) + ham%xi_d_in_ssss(I,J,N,M)) &
   !!                                  - (ham%xi_x_in_ssss(I,J,M,N) + ham%xi_x_in_ssss(I,J,N,M))
   !!               ham%B_ssss(c1,c2) = ham%B_ssss(c1,c2) - 0.5*(ham%lambda_tttt(I,J,M,N) +ham%lambda_tttt(I,J,N,M))
   !!               print*,I,J,M,N,ham%A_ssss(c1,c2)
   !!             
   !!           end do
   !!         end do
   !!       end do
   !!     end do
!!
   !!   
!!
!!
   !! end subroutine construct_A_SS
!!
   !! subroutine construct_A_TTSS()
   !!      implicit none
   !!      integer :: I,J,M,N
   !!      integer :: c1,c2
!!
   !!      allocate(ham%A_ttss(sys%nex*sys%nex,sys%nex*sys%nex))
   !!      allocate(ham%B_ttss(sys%nex*sys%nex,sys%nex*sys%nex))
   !!      ham%A_ttss = 0.0
   !!      ham%B_ttss = 0.0
   !!      c1 =0
   !!      c2 =0
   !!      do I = 1,sys%nex
   !!         do J = 1,sys%nex
   !!          c2 = c2 + 1
   !!          c1 = 0
   !!            do M = 1,sys%nex
   !!                do N = 1,sys%nex
   !!                  c1 = c1 + 1
   !!                  ham%A_ttss(c1,c2) = ham%A_ttss(c1,c2) &
   !!                                     - sqrt(3.0)*(exciton_sys%eigenvalues_s(I) +exciton_sys%eigenvalues_s(J))*(ham%lambda_ttss(I,J,M,N) + ham%lambda_ttss(I,J,N,M))/2 &
   !!                                     - sqrt(3.0)*(ham%xi_d_in_ttss(I,J,M,N)+ham%xi_d_in_ttss(I,J,N,M))/2 &
   !!                                     + sqrt(3.0)*0.5*(ham%xi_x_ttss(I,J,M,N) + ham%xi_x_ttss(I,J,N,M))
   !!                  ham%B_ttss(c1,c2) = ham%B_ttss(c1,c2) -  sqrt(3.0)*0.5*(ham%lambda_ttss(I,J,M,N) +ham%lambda_ttss(I,J,N,M))                  
   !!                  print*,I,J,M,N,ham%A_ttss(c1,c2)
   !!                end do
   !!             
   !!            end do
   !!         end do
   !!       end do
!!
!!
!!
   !! end subroutine construct_A_TTSS
   ! 
  !!  subroutine construct_A_SSTT()
  !!       implicit none
  !!       integer :: I,J,M,N
  !!       integer :: c1,c2
!!
  !!       allocate(ham%A_sstt(sys%nex*sys%nex,sys%nex*sys%nex))
  !!       allocate(ham%B_sstt(sys%nex*sys%nex,sys%nex*sys%nex))
  !!       ham%A_sstt = 0.0
  !!       ham%B_sstt = 0.0
  !!       c1 =0
  !!       c2 =0
  !!       !do I = 1,sys%nex
  !!       !c1 =0
  !!       !c2 =0
  !!       do I = 1,sys%nex
  !!        do J = 1,sys%nex
  !!         c2 = c2 + 1
  !!         c1 = 0
  !!             do M = 1,sys%nex
  !!                 do N = 1,sys%nex
  !!                   c1 = c1 + 1
  !!                   ham%A_sstt(c1,c2) = ham%A_sstt(c1,c2) &
  !!                                      - sqrt(3.0)*(exciton_sys%eigenvalues_t(I) +exciton_sys%eigenvalues_t(J))*(ham%lambda_sstt(I,J,M,N) + ham%lambda_sstt(I,J,N,M))/2 &
  !!                                      - sqrt(3.0)*(ham%xi_d_in_sstt(I,J,M,N)+ham%xi_d_in_sstt(I,J,N,M))/2 &
  !!                                      + sqrt(3.0)*0.5*(ham%xi_x_sstt(I,J,M,N)+ham%xi_x_sstt(I,J,N,M)) &
  !!                                      - sqrt(3.0)*(ham%xi_x_in_sstt(I,J,M,N)+ham%xi_x_in_sstt(I,J,N,M))
  !!                  ham%B_sstt(c1,c2) = ham%B_sstt(c1,c2) -  sqrt(3.0)*0.5*(ham%lambda_sstt(I,J,M,N) +ham%lambda_sstt(I,J,N,M))
  !!                  print*,I,J,M,N,ham%A_sstt(c1,c2) 
  !!                 end do
  !!             end do
  !!         end do
  !!       end do
  !!  end subroutine construct_A_SSTT
!
   ! subroutine construct_AB_TT()
   !   ! A_{MN,IJ} = (E_{M} + E_{N})\delta_{MN,IJ) +  xi_d_tttt(IJMN) + 3/2 * xi_x_tttt(IJMN) &
   !   ! + 1/2* (E_{I}+E_{J})*(lambda_h(IJMN)+lambda_e(IJMN)) &
   !   ! +1/2*(xi_d_in_tttt(IJMN) + xi_d_im_tttt(IJNM)
   !    
 !
   !    implicit none
   !    integer :: I,J,M,N,I_ex,J_ex,M_ex,N_ex
   !    integer :: c1,c2
   !    print*,"hello7"
   !    allocate(ham%A_tttt(sys%nex,sys%nex,sys%nex,sys%nex))
   !    allocate(ham%B_tttt(sys%nex,sys%nex,sys%nex,sys%nex))
   !     print*,"hello8"
   !    ham%A_tttt = 0.0
   !    ham%B_tttt = 0.0
   !    c1 =0
   !    c2 =0
   !    do I = 1,sys%nex
   !       do J = 1,sys%nex
   !        c2 = c2 + 1
   !        c1 = 0
   !          do M = 1,sys%nex
   !             do N = 1,sys%nex
   !              c1 = c1 + 1
 !
   !                if (I==M .and. J==N) then
   !                   ham%A_tttt(I,J,M,N) = ham%A_tttt(I,J,M,N) + (exciton_sys%eigenvalues_t(I) + exciton_sys%eigenvalues_t(J))
   !                   ham%B_tttt(I,J,M,N) = ham%B_tttt(I,J,M,N) + 1
   !                end if
 !
   !                if(I==N .and. J==M)then
   !                   ham%A_tttt(I,J,M,N) = ham%A_tttt(I,J,M,N) + (exciton_sys%eigenvalues_t(I) + exciton_sys%eigenvalues_t(J))
   !                   ham%B_tttt(I,J,M,N) = ham%B_tttt(I,J,M,N) + 1 
   !                end if 
   !                   ham%A_tttt(I,J,M,N) = ham%A_tttt(I,J,M,N) + (ham%xi_d_tttt(I,J,M,N) + ham%xi_d_tttt(I,J,N,M)) &
   !                                   + 1.5 * (ham%xi_x_tttt(I,J,M,N)+ham%xi_x_tttt(I,J,N,M)) &
   !                                   + 0.5 * (exciton_sys%eigenvalues_t(I) + exciton_sys%eigenvalues_t(J)) * (ham%lambda_tttt(I,J,M,N) + ham%lambda_tttt(I,J,N,M)) &
   !                                   + 0.5 * (ham%xi_d_in_tttt(I,J,M,N) + ham%xi_d_in_tttt(I,J,N,M))
   !                
   !                   ham%B_tttt(I,J,M,N) = ham%B_tttt(I,J,M,N) + 0.5*(ham%lambda_tttt(I,J,M,N) +ham%lambda_tttt(I,J,N,M))
   !                print*,I,J,M,N,"TTTT",ham%A_tttt(I,J,M,N),ham%B_tttt(I,J,M,N)
   !              end do
   !          end do
   !        end do
   !      end do
   !   end subroutine construct_AB_TT
 !
   !   subroutine construct_AB_SS()
   !    !A_{MN,IJ} =  (E_{I} + E_{J})*(delta_{I,M}* delta_{JN} &
   !    !           + xi_d_ssss(I,J,M,N)
   !    !           + 0.5*xi_x_ssss(I,J,M,N)
   !    !           - 0.5*(E_{I} + E_{J})*(lambda(I,J,M,N) + lambda(I,J,N,M))
   !    !           - 0.5*(xi_d_in_ssss(I,J,M,N) + xi_d_in_ssss(I,J,N,M))
   !    !           - (xi_x_in_ssss(I,J,M,N) + xi_x_in_ssss(I,J,N,M))
 !
   !    implicit none
   !    integer :: I,J,M,N
   !    integer :: c1,c2
 !
   !    allocate(ham%A_ssss(sys%nex,sys%nex,sys%nex,sys%nex))
   !    allocate(ham%B_ssss(sys%nex,sys%nex,sys%nex,sys%nex))
   !    ham%A_ssss = 0.0
   !    ham%B_ssss = 0.0
 !
   !    c1 =0
   !    c2 =0
   !    do I = 1,sys%nex
   !       do J = 1,sys%nex
   !        c2 = c2 + 1
   !        c1 = 0
   !          do M = 1,sys%nex
   !              do N = 1,sys%nex
   !                c1 = c1 + 1
   !                if (I==M .and. J==N) then
   !                   ham%A_ssss(I,J,M,N) = ham%A_ssss(I,J,M,N) + (exciton_sys%eigenvalues_s(I) + exciton_sys%eigenvalues_s(J))
   !                   ham%B_ssss(I,J,M,N) = ham%B_ssss(I,J,M,N) + 1
   !                end if
   !                if (I==N .and. J==M) then
   !                   ham%A_ssss(I,J,M,N) = ham%A_ssss(I,J,M,N) + (exciton_sys%eigenvalues_s(I) + exciton_sys%eigenvalues_s(J))
   !                   ham%B_ssss(I,J,M,N) = ham%B_ssss(I,J,M,N) + 1
   !                end if
   !                ham%A_ssss(I,J,M,N) = ham%A_ssss(I,J,M,N) + (ham%xi_d_ssss(I,J,M,N)+ham%xi_d_ssss(I,J,N,M)) &
   !                                   + 0.5 * (ham%xi_x_ssss(I,J,M,N)+ham%xi_x_ssss(I,J,N,M)) &
   !                                   - 0.5 * (exciton_sys%eigenvalues_s(I) + exciton_sys%eigenvalues_s(J)) * (ham%lambda_ssss(I,J,M,N) + ham%lambda_ssss(I,J,N,M)) &
   !                                   - 0.5 * (ham%xi_d_in_ssss(I,J,M,N) + ham%xi_d_in_ssss(I,J,N,M)) &
   !                                   - (ham%xi_x_in_ssss(I,J,M,N) + ham%xi_x_in_ssss(I,J,N,M))
   !                ham%B_ssss(I,J,M,N) = ham%B_ssss(I,J,M,N) - 0.5*(ham%lambda_ssss(I,J,M,N) +ham%lambda_ssss(I,J,N,M))
   !                print*,I,J,M,N,"SSSS",ham%A_ssss(I,J,M,N),ham%B_ssss(I,J,M,N)
   !              
   !            end do
   !          end do
   !        end do
   !      end do
 !
   !    
 !
 !
   !  end subroutine construct_AB_SS
 !
   !  subroutine construct_AB_TTSS()
   !       implicit none
   !       integer :: I,J,M,N
   !       integer :: c1,c2
 !
   !       allocate(ham%A_ttss(sys%nex,sys%nex,sys%nex,sys%nex))
   !       allocate(ham%B_ttss(sys%nex,sys%nex,sys%nex,sys%nex))
   !       ham%A_ttss = 0.0
   !       ham%B_ttss = 0.0
   !       c1 =0
   !       c2 =0
   !       do I = 1,sys%nex
   !          do J = 1,sys%nex
   !           c2 = c2 + 1
   !           c1 = 0
   !             do M = 1,sys%nex
   !                 do N = 1,sys%nex
   !                   c1 = c1 + 1
   !                   ham%A_ttss(I,J,M,N) = ham%A_ttss(I,J,M,N) &
   !                                      - sqrt(3.0)*(exciton_sys%eigenvalues_s(I) +exciton_sys%eigenvalues_s(J))*(ham%lambda_ttss(I,J,M,N) + ham%lambda_ttss(I,J,N,M))/2 &
   !                                      - sqrt(3.0)*(ham%xi_d_in_ttss(I,J,M,N)+ham%xi_d_in_ttss(I,J,N,M))/2 &
   !                                      + sqrt(3.0)*0.5*(ham%xi_x_ttss(I,J,M,N) + ham%xi_x_ttss(I,J,N,M))
   !                   ham%B_ttss(I,J,M,N) = ham%B_ttss(I,J,M,N) -  sqrt(3.0)*0.5*(ham%lambda_ttss(I,J,M,N) +ham%lambda_ttss(I,J,N,M))                  
   !                   print*,I,J,M,N,"ttss",ham%A_ttss(I,J,M,N),ham%B_ttss(I,J,M,N)
   !                 end do
   !              
   !             end do
   !          end do
   !        end do
 !
 !
 !
   !  end subroutine construct_AB_TTSS
   !  
   !  subroutine construct_AB_SSTT()
   !       implicit none
   !       integer :: I,J,M,N
   !       integer :: c1,c2
 !
   !       allocate(ham%A_sstt(sys%nex,sys%nex,sys%nex,sys%nex))
   !       allocate(ham%B_sstt(sys%nex,sys%nex,sys%nex,sys%nex))
   !       ham%A_sstt = 0.0
   !       ham%B_sstt = 0.0
   !       c1 =0
   !       c2 =0
   !       !do I = 1,sys%nex
   !       !c1 =0
   !       !c2 =0
   !       do I = 1,sys%nex
   !        do J = 1,sys%nex
   !         c2 = c2 + 1
   !         c1 = 0
   !             do M = 1,sys%nex
   !                 do N = 1,sys%nex
   !                   c1 = c1 + 1
   !                   ham%A_sstt(I,J,M,N) = ham%A_sstt(I,J,M,N) &
   !                                      - sqrt(3.0)*(exciton_sys%eigenvalues_t(I) +exciton_sys%eigenvalues_t(J))*(ham%lambda_sstt(I,J,M,N) + ham%lambda_sstt(I,J,N,M))/2 &
   !                                      - sqrt(3.0)*(ham%xi_d_in_sstt(I,J,M,N)+ham%xi_d_in_sstt(I,J,N,M))/2 &
   !                                      + sqrt(3.0)*0.5*(ham%xi_x_sstt(I,J,M,N)+ham%xi_x_sstt(I,J,N,M)) &
   !                                      - sqrt(3.0)*(ham%xi_x_in_sstt(I,J,M,N)+ham%xi_x_in_sstt(I,J,N,M))
   !                  ham%B_sstt(I,J,M,N) = ham%B_sstt(I,J,M,N) -  sqrt(3.0)*0.5*(ham%lambda_sstt(I,J,M,N) +ham%lambda_sstt(I,J,N,M))
   !                  print*,"sstt",I,J,M,N,ham%A_sstt(I,J,M,N),ham%B_sstt(I,J,M,N) 
   !                 end do
   !             end do
   !         end do
   !       end do
   !  end subroutine construct_AB_SSTT
 !w!rite a subroutine to read lambda.dat file and use it to store lambda in ham%lambda and ham%lambda_e
  !!  subroutine read_lambda()
!
  !!      implicit none
  !!      integer :: I,J,M,N
  !!      integer :: I_ex,J_ex,M_ex,N_ex
  !!      integer :: c1,c2
  !!      double precision :: lambda1,lambda2
  !!      allocate(ham%lambda(sys%nex_max,sys%nex_max,sys%nex_max,sys%nex_max))
  !!      allocate(ham%lambda_e(sys%nex_max,sys%nex_max,sys%nex_max,sys%nex_max))
  !!      ham%lambda = 0.0
  !!      ham%lambda_e = 0.0
  !!      open(UNIT=10, FILE="lambda.dat", STATUS="OLD")
  !!      do M = 1,40
  !!          do N = 1,40
  !!              do I = 1,40
  !!                  do J = 1,40
  !!                      read(10,*) I_ex,J_ex,M_ex,N_ex,lambda1
  !!                      read(10,*) lambda2
  !!                      ham%lambda(I,J,M,N) = lambda1
  !!                      ham%lambda_e(I,J,M,N) = lambda2
  !!                      !print*,"lambda",I,J,M,N,ham%lambda(I,J,M,N),ham%lambda_e(I,J,M,N)
  !!                  end do
  !!              end do
  !!          end do
  !!      end do
  !!      close(10)
!!
!!
!!
  !!   end subroutine read_lambda
  !!   
   ! ! subroutine read_xi_xi_in()
!!
   ! !    implicit none
   ! !    integer :: I,J,M,N
   ! !    integer :: I_ex,J_ex,M_ex,N_ex
   ! !    integer :: c1,c2
   ! !    double precision :: xi1,xi2
   ! !    allocate(ham%xi(sys%nex_max,sys%nex_max,sys%nex_max,sys%nex_max))
   ! !    ham%xi = 0.0
   ! !    allocate(ham%xi_in(sys%nex_max,sys%nex_max,sys%nex_max,sys%nex_max))
   ! !    open(UNIT=10, FILE="xi.dat", STATUS="OLD")
   ! !    do M = 1,40
   ! !        do N = 1,40
   ! !            do I = 1,40
   ! !                do J = 1,40
!!
   ! !                    read(10,*) I_ex,J_ex,M_ex,N_ex,xi1
   ! !                    read(10,*) xi2
   ! !                    ham%xi(I_ex,J_ex,M_ex,N_ex) = xi1
   ! !                    ham%xi_in(I_ex,J_ex,M_ex,N_ex) = xi2
   ! !                   ! print*,"xi",I,J,M,N,ham%xi(I,J,M,N),ham%xi_in(I,J,M,N)
!!
   ! !                end do
   ! !            end do
   ! !        end do
   ! !    end do
   ! !    close(10)
   ! ! end subroutine read_xi_xi_in
   !
!
   ! 
   ! subroutine compute_A_and_B()
!
   !     !implicit none
   !     !integer :: I,J,M,N
   !     !integer :: c1,c2
!
   !     !allocate(ham%xi_in(sys%nex,sys%nex,sys%nex,sys%nex))
   !     !ham%xi_in = 0.0
!
   !     !call compute_xi_in()
   !     !allocate(ham%xi_in_min_xi_dir_x_lam_eph(sys%nex,sys%nex,sys%nex,sys%nex))
   !     !ham%xi_in_min_xi_dir_x_lam_eph = 0.0
   !     !print*,"hello"
   !     !call compute_xi_in_min_xi_dir_x_lam_eph()
   !     !allocate(ham%lambda_h_ERS_lambda_eph(sys%nex,sys%nex,sys%nex,sys%nex))
   !     !ham%lambda_h_ERS_lambda_eph = 0.0
   !     !call compute_lambda_h_ERS_lambda_eph()
   !     ! c1 = 0
   !     ! c2 = 0
   !     ! allocate(ham%A(ham%l,ham%l))
   !     ! allocate(ham%B(ham%l,ham%l))
   !     ! ham%A = 0.0
   !     ! ham%B = 0.0
   !     ! ham%xi = 0.0
   !     ! ham%xi_in=0.0
   !         
   !     ! do M = 1,sys%nex
   !     !     do N = M,sys%nex
   !     !         c2 = 0
   !     !         c1 = c1 + 1
   !     !         do I = 1,sys%nex
   !     !             do J = I,sys%nex
   !     !                 c2 = c2 + 1
   !     !                 if(I==M .and. J==N) then
   !     !                     if(I==J) then
   !     !                        ham%A(c1,c2) = 2*(exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
   !     !                        ham%B(c1,c2) = 2
   !     !                     end if
   !     !                     if(I.ne.J) then
   !     !                        ham%A(c1,c2) = (exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
   !     !                        ham%B(c1,c2) = 1
   !     !                     end if
   !     !                 end if
   !     !                 print*,I,J,M,N,ham%lambda(I,J,M,N),ham%lambda(I,J,N,M)
   !     !                 ham%A(c1,c2) = ham%A(c1,c2) + ham%xi(I,J,M,N) &
   !     !                                - ham%xi_in(I,J,M,N) - ham%xi_in(I,J,N,M) &
   !     !                                -((exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))*(ham%lambda(I,J,M,N)+ham%lambda(I,J,N,M)))                    
   !     !                 ham%B(c1,c2) = ham%B(c1,c2) - ham%lambda(I,J,M,N) -ham%lambda(I,J,N,M)
   !     !                 !print*,c1,c2,ham%A(c1,c2),ham%B(c1,c2)
   !                     
   !                     
   !     !             end do
   !     !         end do
   !     !     end do
   !     ! end do
   !     
   ! end subroutine compute_A_and_B
!
  !
  !!subroutine compute_A_and_B_double_basis1()
!!
  !!      implicit none
  !!      integer :: I,J,M,N
  !!      integer :: c1,c2
  !!      
!!
  !!      allocate(ham%xi_in(sys%nex,sys%nex,sys%nex,sys%nex))
  !!      !allocate(ham%xi_dir_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
  !!      !allocate(ham%xi_in_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
  !!      !allocate(ham%xi_out(sys%nex,sys%nex,sys%nex,sys%nex))
  !!      !allocate(ham%xi_out_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
  !!      !allocate(ham%E_MNRS_lambda_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
  !!      !allocate(ham%lambda_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
  !!      ham%xi_in = 0.0
  !!      !ham%xi_out = 0.0
  !!      !ham%xi_dir_lambda_e = 0.0
  !!      !ham%xi_in_lambda_e = 0.0
  !!      !ham%xi_out_lambda_e = 0.0
  !!      !ham%E_MNRS_lambda_lambda_e = 0.0
  !!      !ham%lambda_lambda_e = 0.0
  !!      call compute_xi_in()
  !!      c1 = 0
  !!      c2 = 0
  !!      allocate(ham%A(ham%l,ham%l))
  !!      allocate(ham%B(ham%l,ham%l))
  !!      ham%A = 0.0
  !!      ham%B = 0.0
  !!      !ham%xi = 0.0
  !!      !ham%xi_in=0.0
  !!          
  !!      do M = 1,sys%nex
  !!          do N = 1,sys%nex
  !!              c2 = 0
  !!              c1 = c1 + 1
  !!              do I = 1,sys%nex
  !!                  do J = 1,sys%nex
  !!                      c2 = c2 + 1
  !!                      if(I==M .and. J==N) then
  !!                          ham%A(c1,c2) = ham%A(c1,c2) + (exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
  !!                          ham%B(c1,c2) = ham%B(c1,c2) + 1
  !!                      end if
  !!                      !if(I==N .and. J==M) then
  !!                       !   ham%A(c1,c2) = ham%A(c1,c2) + (exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
  !!                        !  ham%B(c1,c2) = ham%B(c1,c2) + 1
  !!                      !end if
  !!                      !print*,I,J,M,N,ham%lambda(I,J,M,N),ham%lambda(I,J,N,M)
  !!                      ham%A(c1,c2) = ham%A(c1,c2) &
  !!                                     + ham%xi(I,J,M,N) &!+ ham%xi(I,J,N,M) &
  !!                                     - ham%xi_in(I,J,M,N) & !- ham%xi_in(I,J,N,M) &
  !!                                     -((exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))*(ham%lambda(I,J,M,N))) !+ham%lambda(I,J,N,M)))                    
  !!                      ham%B(c1,c2) = ham%B(c1,c2) - ham%lambda(I,J,M,N) !-ham%lambda(I,J,N,M)
  !!                      !print*,c1,c2,ham%A(c1,c2),ham%B(c1,c2)
  !!                      
  !!                      
  !!                  end do
  !!              end do
  !!          end do
  !!      end do
  !!      
  !!  end subroutine compute_A_and_B_double_basis1
!
   !! subroutine compute_A_and_B_double_basis()
!!
   !!     implicit none
   !!     integer :: I,J,M,N
   !!     integer :: c1,c2
   !!     
!!
   !!     allocate(ham%xi_in(sys%nex,sys%nex,sys%nex,sys%nex))
   !!     !allocate(ham%xi_dir_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
   !!     !allocate(ham%xi_in_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
   !!     !allocate(ham%xi_out(sys%nex,sys%nex,sys%nex,sys%nex))
   !!     !allocate(ham%xi_out_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
   !!     !allocate(ham%E_MNRS_lambda_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
   !!     !allocate(ham%lambda_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
   !!     ham%xi_in = 0.0
   !!     !ham%xi_out = 0.0
   !!     !ham%xi_dir_lambda_e = 0.0
   !!     !ham%xi_in_lambda_e = 0.0
   !!     !ham%xi_out_lambda_e = 0.0
   !!     !ham%E_MNRS_lambda_lambda_e = 0.0
   !!     !ham%lambda_lambda_e = 0.0
   !!     call compute_xi_in()
   !!     c1 = 0
   !!     c2 = 0
   !!     allocate(ham%A(ham%l,ham%l))
   !!     allocate(ham%B(ham%l,ham%l))
   !!     ham%A = 0.0
   !!     ham%B = 0.0
   !!     !ham%xi = 0.0
   !!     !ham%xi_in=0.0
   !!         
   !!     do M = 1,sys%nex
   !!         do N = 1,sys%nex
   !!             c2 = 0
   !!             c1 = c1 + 1
   !!             do I = 1,sys%nex
   !!                 do J = 1,sys%nex
   !!                     c2 = c2 + 1
   !!                     if(I==M .and. J==N) then
   !!                         ham%A(c1,c2) = ham%A(c1,c2) + (exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
   !!                         ham%B(c1,c2) = ham%B(c1,c2) + 1
   !!                     end if
   !!                     if(I==N .and. J==M) then
   !!                         ham%A(c1,c2) = ham%A(c1,c2) + (exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
   !!                         ham%B(c1,c2) = ham%B(c1,c2) + 1
   !!                     end if
   !!                     !print*,I,J,M,N,ham%lambda(I,J,M,N),ham%lambda(I,J,N,M)
   !!                     ham%A(c1,c2) = ham%A(c1,c2) &
   !!                                    + ham%xi(I,J,M,N) + ham%xi(I,J,N,M) &
   !!                                    - ham%xi_in(I,J,M,N) - ham%xi_in(I,J,N,M) &
   !!                                    -((exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))*(ham%lambda(I,J,M,N)+ham%lambda(I,J,N,M)))                    
   !!                     ham%B(c1,c2) = ham%B(c1,c2) - ham%lambda(I,J,M,N) -ham%lambda(I,J,N,M)
   !!                     print*,c1,c2,ham%A(c1,c2),ham%B(c1,c2)
   !!                     
   !!                     
   !!                 end do
   !!             end do
   !!         end do
   !!     end do
   !!     
   !! end subroutine compute_A_and_B_double_basis
!!
!
 ! !  subroutine compute_A_ortho()
!!
 ! !      implicit none
 ! !      integer :: I,J,M,N
 ! !      integer :: c1,c2
!!
 ! !      !allocate(ham%xi_in(sys%nex,sys%nex,sys%nex,sys%nex))
 ! !      !ham%xi_in = 0.0
 ! !      !call compute_xi_in()        
 ! !   
 ! !      c1 = 0
 ! !      c2 = 0
 ! !      allocate(ham%A(sys%nex_max,sys%nex_max,sys%nex_max,sys%nex_max))
 ! !      allocate(ham%B(sys%nex_max,sys%nex_max,sys%nex_max,sys%nex_max))
 ! !      ham%A = 0.0
 ! !      ham%B = 0.0
 ! !      do M = 1,sys%nex_max
 ! !          do N = 1,sys%nex_max
 ! !              do I = 1,sys%nex_max
 ! !                  do J = 1,sys%nex_max
 ! !                      c1 = c1 + 1
 ! !                      if(I==M .and. J==N) then
 ! !                          ham%A(I,J,M,N) = ham%A(I,J,M,N) + (exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
 ! !                          ham%B(I,J,M,N) = ham%B(I,J,M,N) + 1
 ! !                      end if
 ! !                      if(I==N .and. J==M) then
 ! !                          ham%A(I,J,M,N) = ham%A(I,J,M,N) + (exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
 ! !                          ham%B(I,J,M,N) = ham%B(I,J,M,N) + 1
 ! !                      end if
 ! !                      ham%A(I,J,M,N) = ham%A(I,J,M,N) &
 ! !                                  + ham%xi(I,J,M,N) + ham%xi(I,J,N,M) &
 ! !                                  - ham%xi_in(I,J,M,N) - ham%xi_in(I,J,N,M) &
 ! !                                  -((exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))*(ham%lambda(I,J,M,N)+ham%lambda(I,J,N,M))) 
 ! !                      ham%B(I,J,M,N) = ham%B(I,J,M,N) - ham%lambda(I,J,M,N) -ham%lambda(I,J,N,M)
 ! !                      !print*,c1,ham%A(c1)
 ! !                  end do
 ! !              end do
 ! !          end do
 ! !      end do
 ! !  end subroutine compute_A_ortho
!!
!!
   ! subroutine compute_H_ortho()
!
   !     implicit none
   !     integer :: v1,v2
   !     integer :: M_s,N_s,I_s,J_s
   !     integer :: c1,c2,s1,s2
   !     double precision,allocatable :: c_MNIJ(:)
!
   !     allocate(ham%H(sys%num_ortho,sys%num_ortho)) ! Hamiltonian projected on to the orthonormal basis
   !     ham%H = 0.0 !
   !    ! allocate(ham%H1(sys%northo,sys%northo))
   !   ! ham%H1 = 0.0
   !   ! allocate(c_MNIJ(ham%l*ham%l))
 ! !      
   !     do v1 = 1,sys%num_ortho
   !         do v2 = 1,sys%num_ortho       
   !            !call compute_c_MNIJ(v1,v2,c_MNIJ)
   !             !ham%H(v1,v2) =  ham%H(v1,v2) + dot_product(c_MNIJ,ham%A)
   !             !ham%H1(v1,v2) = ham%H1(v1,v2) + dot_product(c_MNIJ,ham%B)
   !             !print*,v1,v2,ham%H(v1,v2),ham%H1(v1,v2)
   !             do c1 = 1,sys%ncutoff
   !                 M_s = exciton_sys%I_arr(c1)  
   !                 N_s = exciton_sys%J_arr(c1)
   !                 s1 = exciton_sys%spin_ind(c1)
   !                 do c2 = 1,sys%ncutoff
   !                     I_s = exciton_sys%I_arr(c2)
   !                     J_s = exciton_sys%J_arr(c2)
   !                     s2 = exciton_sys%spin_ind(c2)
   !                     if (s1 == 1 .and. s2 ==1) then
   !                        ham%H(v1,v2) = ham%H(v1,v2) + ham%c(c1,v1)*ham%c(c2,v2)*ham%A_tttt(I_s,J_s,M_s,N_s)
   !                     end if
   !                     if (s1 == 2 .and. s2 ==2) then
   !                       ham%H(v1,v2) = ham%H(v1,v2) + ham%c(c1,v1)*ham%c(c2,v2)*ham%A_ssss(I_s,J_s,M_s,N_s)
   !                    end if
   !                    if (s1 ==1 .and. s2 ==2 )then
   !                       ham%H(v1,v2) = ham%H(v1,v2) + ham%c(c1,v1)*ham%c(c2,v2)*ham%A_ttss(I_s,J_s,M_s,N_s)
   !                    end if
   !                    if (s1 ==2 .and. s2 ==1 )then
   !                       ham%H(v1,v2) = ham%H(v1,v2) + ham%c(c1,v1)*ham%c(c2,v2)*ham%A_sstt(I_s,J_s,M_s,N_s)
   !                    end if 
   !                     !print*,"c",v2,I_s,J_s,ham%c(c2,v2)   
   !                 end do
   !             end do  
 ! !              !print*,"H",ham%H(v1,v2),v1,v2
   !         end do
   !     end do
!!
!!
!!
   !end subroutine compute_H_ortho
!
   !subroutine compute_H_ortho1()
!
   !   implicit none
   !   integer :: v1,v2
   !   integer :: M_s,N_s,I_s,J_s
   !   integer :: c1,c2,s1,s2
   !   double precision,allocatable :: c_MNIJ(:)
!
   !   allocate(ham%H(sys%num_ortho,sys%num_ortho)) ! Hamiltonian projected on to the orthonormal basis
   !   ham%H = 0.0 
!
   !  ! allocate(ham%H1(sys%northo,sys%northo))
   ! ! ham%H1 = 0.0
   ! ! allocate(c_MNIJ(ham%l*ham%l))
!  !     
   !   do v1 = 1,sys%num_ortho
   !       do v2 = 1,sys%num_ortho       
   !          !call compute_c_MNIJ(v1,v2,c_MNIJ)
   !           !ham%H(v1,v2) =  ham%H(v1,v2) + dot_product(c_MNIJ,ham%A)
   !           !ham%H1(v1,v2) = ham%H1(v1,v2) + dot_product(c_MNIJ,ham%B)
   !           !print*,v1,v2,ham%H(v1,v2),ham%H1(v1,v2)
   !         c1 =0
   !         c2 =0
   !   
   !          do M_s = 1,sys%nex
   !            do N_s = 1,4!sys%nex
   !               c1 = c1 +1
   !               c2 = 0
!
   !                do I_s = 1,4!sys%nex
   !                 do J_s = 1,4!sys%nex
   !                  c2 = c2 +1
   !                   !if (s1 == 1 .and. s2 ==1) then
   !                      ham%H(v1,v2) = ham%H(v1,v2) + ham%c(c1,v1)*ham%c(c2,v2)*ham%A_tttt(I_s,J_s,M_s,N_s)
   !                      !print*,"c",I_s,J_s,M_s,N_s,ham%c(c1,v1)*ham%c(c2,v2)*ham%A_tttt(I_s,J_s,M_s,N_s)
   !                   !end if
   !                 !  if (s1 == 2 .and. s2 ==2) then
   !                 !    ham%H(v1,v2) = ham%H(v1,v2) + ham%c(c1,v1)*ham%c(c2,v2)*ham%A_ssss(I_s,J_s,M_s,N_s)
   !                 ! end if
   !                 ! if (s1 ==1 .and. s2 ==2 )then
   !                  !   ham%H(v1,v2) = ham%H(v1,v2) + ham%c(c1,v1)*ham%c(c2,v2)*ham%A_ttss(I_s,J_s,M_s,N_s)
   !                  !end if
   !                  !if (s1 ==2 .and. s2 ==1 )then
   !                  !   ham%H(v1,v2) = ham%H(v1,v2) + ham%c(c1,v1)*ham%c(c2,v2)*ham%A_sstt(I_s,J_s,M_s,N_s)
   !                  !end if 
   !                   
   !                end do
   !               end do 
   !            end do
   !         end do
   !        
   !            print*,"H",ham%H(v1,v2),v1,v2
   !    end do
   !   end do
!!
!!
!!
 !end subroutine compute_H_ortho1
 ! !  
 ! !  subroutine compute_c_MNIJ(v1,v2,c_MNIJ)
!!
!!
 ! !      implicit none
!!
!!
!!
 ! !      integer,intent(in) :: v1,v2
 ! !      integer :: M,N,I,J,c1,c2,c3
 ! !      double precision,dimension(ham%l*ham%l),intent(out) :: c_MNIJ(:)
 ! !      c1 = 0
 ! !      c2 = 0
 ! !      c3 = 0
 ! !      do M =  1,sys%nex
 ! !          do N = 1,sys%nex
 ! !              c2 = c2 + 1
 ! !              c3 = 0
 !               do I = 1,sys%nex
 !                   do J = 1,sys%nex
 !                       c1 = c1 + 1
 !                       c3 = c3 + 1
 !                       c_MNIJ(c1) = ham%c(c2,v1)*ham%c(c3,v2)
 !                      ! print*,c1,c2,c3,c_MNIJ(c1)
 !                   end do
 !               end do
 !           end do
 !       end do
!
 !   end subroutine compute_c_MNIJ
!
 !   !subroutine compute_A_and_B_sym()
!!
 !   !    implicit none
 !   !    integer :: I,J,M,N
 !   !    integer :: c1,c2
    !    
!
    !    allocate(ham%xi_in(sys%nex,sys%nex,sys%nex,sys%nex))
    !    allocate(ham%xi_dir_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
    !    allocate(ham%xi_in_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
    !    allocate(ham%xi_out(sys%nex,sys%nex,sys%nex,sys%nex))
    !    allocate(ham%xi_out_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
    !    allocate(ham%E_MNRS_lambda_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
    !    allocate(ham%lambda_lambda_e(sys%nex,sys%nex,sys%nex,sys%nex))
    !    ham%xi_in = 0.0
    !    ham%xi_out = 0.0
    !    ham%xi_dir_lambda_e = 0.0
    !    ham%xi_in_lambda_e = 0.0
    !    ham%xi_out_lambda_e = 0.0
    !    ham%E_MNRS_lambda_lambda_e = 0.0
    !    ham%lambda_lambda_e = 0.0
    !    call compute_xi_in()
    !    call compute_xi_out()
    !    call compute_xi_dir_lambda_e()
    !    call compute_xi_in_lambda_e()
    !    call compute_xi_out_lambda_e()
    !    call compute_E_MNRS_lambda_lambda_e()
    !    call compute_lambda_lambda_e()
    !    c1 = 0
    !    c2 = 0
    !    allocate(ham%A(ham%l,ham%l))
    !    allocate(ham%B(ham%l,ham%l))
    !    ham%A = 0.0
    !    ham%B = 0.0
    !    !ham%xi = 0.0
    !    !ham%xi_in=0.0
    !        
    !    do M = 1,sys%nex
    !        do N = 1,sys%nex
    !            c2 = 0
    !            c1 = c1 + 1
    !            do I = 1,sys%nex
    !                do J = 1,sys%nex
    !                    c2 = c2 + 1
    !                    if(I==M .and. J==N) then
    !                        ham%A(c1,c2) = ham%A(c1,c2) + (exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))
    !                        ham%B(c1,c2) = ham%B(c1,c2) + 1
    !                    end if
    !                    
    !                    !print*,I,J,M,N,ham%lambda(I,J,M,N),ham%lambda(I,J,N,M)
    !                    ham%A(c1,c2) = ham%A(c1,c2) &
    !                                   -  ((exciton_sys%eigenvalues(M) + exciton_sys%eigenvalues(N))*ham%lambda_e(I,J,M,N)) &
    !                                   !-((exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J)+exciton_sys%eigenvalues(M) + exciton_sys%eigenvalues(N))*(ham%lambda_e(I,J,M,N)/2))&
    !                                   + ham%xi(I,J,M,N) - (1/2*ham%xi_in(I,J,M,N)) - (1/2*ham%xi_out(I,J,M,N)) &
    !                                   -((exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J)+exciton_sys%eigenvalues(M) + exciton_sys%eigenvalues(N))*(ham%lambda(I,J,M,N)/2))&
    !                                   - ham%xi_dir_lambda_e(I,J,M,N) + (1/2*ham%xi_in_lambda_e(I,J,M,N)) + (1/2*ham%xi_out_lambda_e(I,J,M,N)) &
    !                                   + (1/2*ham%E_MNRS_lambda_lambda_e(I,J,M,N))
    !                    ham%B(c1,c2) = ham%B(c1,c2) - ham%lambda(I,J,M,N) - ham%lambda_e(I,J,M,N) + ham%lambda_lambda_e(I,J,M,N) 
    !                    print*,c1,c2,ham%A(c1,c2),ham%B(c1,c2)
    !                                            
    !                end do
    !            end do
    !        end do
    !    end do
    !    
    !end subroutine compute_A_and_B_sym

 !   subroutine compute_A_and_B_with_parity()
!
 !       integer :: I,J,M,N
 !       integer :: c1,c2
!
!
 !       allocate(ham%xi_in(sys%nex,sys%nex,sys%nex,sys%nex))
 !       allocate(ham%xi_in_min_xi_dir_x_lam_eph(sys%nex,sys%nex,sys%nex,sys%nex))
 !       allocate(ham%lambda_h_ERS_lambda_eph(sys%nex,sys%nex,sys%nex,sys%nex))
!
 !       ham%xi_in = 0.0
 !       ham%xi_in_min_xi_dir_x_lam_eph = 0.0
 !       ham%lambda_h_ERS_lambda_eph = 0.0
!
 !       call compute_xi_in()
 !       call compute_xi_in_min_xi_dir_x_lam_eph()
 !       call compute_lambda_h_ERS_lambda_eph()
!
 !       allocate(ham%A(ham%l,ham%l))
 !       allocate(ham%B(ham%l,ham%l))
 !       ham%A = 0.0
 !       ham%B = 0.0
!
!
 !       c1 = 0
 !       c2 = 0
!
 !       open(unit=11,file='Ham.dat',status='replace')
 !       do M = 1,sys%nex
 !           do N = 1,sys%nex
 !               c2 = 0
 !               c1 = c1 + 1
 !               do I = 1,sys%nex
 !                   do J = 1,sys%nex
 !                       c2 = c2 + 1
 !                       if ( I == M .and. J == N ) then
 !                           ham%A(c1,c2) = ham%A(c1,c2) + 2*(exciton_sys%eigenvalues(M) + exciton_sys%eigenvalues(N))
 !                           ham%B(c1,c2) = ham%B(c1,c2) + 4    
 !                       end if
 !                       
 !                       ham%A(c1,c2) = ham%A(c1,c2) &
 !                                      - ((exciton_sys%eigenvalues(M) + exciton_sys%eigenvalues(N))*(ham%lambda(I,J,M,N)+ham%lambda_e(I,J,M,N))) &
 !                                      + 2*(ham%xi(I,J,M,N) - ham%xi_in(I,J,M,N)) &
 !                                      - 2*(exciton_sys%eigenvalues(I) + exciton_sys%eigenvalues(J))*(ham%lambda(I,J,M,N)) & 
 !                                      + ham%xi_in_min_xi_dir_x_lam_eph(I,J,M,N) &
 !                                      + ham%lambda_h_ERS_lambda_eph(I,J,M,N)
 !                       ham%B(c1,c2) = ham%B(c1,c2) - 3*ham%lambda(I,J,M,N) - ham%lambda_e(I,J,M,N)
 !                       write(11,*) c1,c2,ham%A(c1,c2),ham%B(c1,c2)
!
 !                   end do
 !               end do
 !           end do
 !       end do
!
 !       close(unit=11)
!
!
 !   end subroutine compute_A_and_B_with_parity
!
!
!
 !   
!
 !   subroutine compute_A_and_B_vectorized()
!
 !       implicit none
 !       integer :: I,J,M,N
 !       integer :: x,c,s,t,z,y
 !       double precision :: lambda,sum, vec_start,vec_end,loop_start,loop_end
 !       double precision,allocatable,dimension(:,:) :: left_eigp,right_eigp 
 !       double precision, allocatable,dimension(:,:) :: A_I,A_J
 !       double precision,allocatable,dimension(:):: B_I,B_J
 !       
 !       ! \lambda(M,N,I,J) = \sum_{c,c',v,v'} left_eigp(c,c',v,v')* right_eigp(c,c',v,v') 
 !       != inner product(left_eigp,right_eigp)
 !       !right_eigp = outer prduct(A^{I}_cv,A^{J}_c'v')
 !       
 !       allocate(right_eigp(sys%nc*sys%nv,sys%nc*sys%nv))
 !       !allocate(right_eigp(sys%nc,sys%nv,sys%nc,sys%nv))
 !      !allocate(A_I(sys%nc,sys%nv))
 !       !allocate(A_J(sys%nc,sys%nv))
 !       allocate(B_I(sys%nc*sys%nv))
 !       !allocate(B_J(sys%nc*sys%nv))
 !       !debug
 !       z =sys%nv*sys%nc
 !       sys%nex =1
 !       do M = 1,sys%nex
 !           do N = M,sys%nex
 !               do I = 1,sys%nex
 !                   do J = I,sys%nex
 !                       !A_I(:,:) = exciton_sys%eigenvectors(I,:,:)
 !                       !A_J(:,:) = exciton_sys%eigenvectors(J,:,:)
 !                       !print*,A_I,"A_I"
 !                       !B_I = reshape(exciton_sys%eigenvectors(I,:,:),[sys%nv*sys%nc])
 !                       !B_J = reshape(exciton_sys%eigenvectors(J,:,:),[sys%nv*sys%nc])
 !                       !print*,B_I,"B_I"
 !                       
 !                       !result = matmul(reshape(a, [n,1]), reshape(b, [1,n]))
 !                       !print*,size(matmul(reshape(reshape(exciton_sys%eigenvectors(I,:,:),[sys%nv*sys%nc]),[x,1]),reshape(B_J,[1,x])))
 !                       call cpu_time(vec_start)
 !                       right_eigp = matmul(reshape(reshape(exciton_sys%eigenvectors(I,:,:),[sys%nv*sys%nc]),[z,1]), &
 !                         reshape(reshape(exciton_sys%eigenvectors(J,:,:),[sys%nv*sys%nc]),[1,z]))
 !                      ! print*,size(right_eigp)
 !                       !print*,reshape(B_I,[x,1]),"right_eigp"
 !                       
 !                                           !reshape(exciton_sys%eigenvectors(J,:,:),(/1,sys%nv*sys%nc/)) )
 !                       !left_eigp = matmul(reshape(exciton_sys%eigenvectors(M,:,:),(/sys%nv*sys%nc,1/)), &
 !                                     !      reshape(exciton_sys%eigenvectors(N,:,:),(/1,sys%nv*sys%nc/)))
 !                       !right_eigp = reshape(right_eigp,(/sys%nc*sys%nv*sys%nc*sys%nv/))                    
 !                       !left_eigp = reshape(left_eigp,(/sys%nc*sys%nv*sys%nc*sys%nv/))
 !                       lambda = dot_product(reshape(reshape(right_eigp,[sys%nv,sys%nc,sys%nv,sys%nc],order=[1,2,3,4]),[sys%nv*sys%nc*sys%nv*sys%nc]), &
 !                                    reshape(reshape(right_eigp,[sys%nv,sys%nc,sys%nv,sys%nc],order=[3,2,1,4]),[sys%nv*sys%nc*sys%nv*sys%nc]))
 !                       call cpu_time(vec_end)
 !                       !print*,lambda,"lambda"
 !                       B_I = reshape(exciton_sys%eigenvectors(J,:,:),[sys%nv*sys%nc]) 
 !                       sum =0
 !                       call cpu_time(loop_start)
 !                       do t=1,sys%nc ! c'
 !                           do s = 1,sys%nv !v'
 !                                do y = 1,sys%nc !c 
 !                                  do x = 1, sys%nv !v
 !                                       sum = sum + (exciton_sys%eigenvectors(I,x,y)*exciton_sys%eigenvectors(J,s,t) &
 !                                                   *exciton_sys%eigenvectors(M,s,y)*exciton_sys%eigenvectors(N,x,t)) 
 !                                       !print*,exciton_sys%eigenvectors(I,x,t)*exciton_sys%eigenvectors(J,s,y),right_eigp(s,t,x,y)
 !                                       !print*,s,t,exciton_sys%eigenvectors(I,s,t),B_I(c)
 !                                   end do
 !                               end do
 !                           end do
 !                       end do
 !                       call cpu_time(loop_end)
 !                      ! print*,sum,"sum"
 !                   end do
 !               end do
 !           end do
 !       end do
 !       
 !       print*,loop_end - loop_start,"loop"
 !       print*,vec_end - vec_start,"vec"
 !   end subroutine compute_A_and_B_vectorized
 !   
end module compute_AB
