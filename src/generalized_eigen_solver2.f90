module solve_generalized_eigenvalue_problem
   use global_variables
   use compute_lambda_and_xi
   use compute_AB
   use kpoints

   !use lapacks generalized eigenvalue solver dggevx, with A = ham%A and B = ham%B and m=ham%l

   implicit none


   contains


   subroutine solve_eigenvalue_overlap()
    
   integer :: ilo,ihi
   integer :: I_s,J_s,M_s,N_s
   integer :: s1,s2   !ss or tt state
   integer :: v1,v2 
   integer :: i,j,alpha,bet,z,k,s,a,b,r

   integer :: I_ex,J_ex,M_ex,N_ex,n,m
   integer :: c1,c2,c3,c4,c
   double precision :: ellow,elup,abstol
   integer :: nfound,info,iup,ilow
   complex(kind=8) :: sum,sum2,sum3,sum4
   complex(kind=8) :: A_MN,A_MNe,A_MNh,A_NM,IJ_amp
   integer :: ik1,ik2,iq,ip,imQ_l,imQ_r,iQ_l,iQ_r
   double precision :: q(3),mQ_l(3),p(3),mQ_r(3)
   double precision,allocatable,dimension(:,:) :: cMN,cijab
   double precision,allocatable,dimension(:,:) :: old
   integer,allocatable :: iarray(:),jarray(:),aarray(:),barray(:)
   double precision , allocatable :: c_ijab(:)
   CHARACTER(100) :: line,res1
   CHARACTER, DIMENSION(100) :: res
   integer :: status
   integer :: p1,p2,ind,v
   double precision :: x,c_v_MN
   double precision  ::Q_l(3)
   integer,dimension(:),allocatable :: ifail
   complex(kind=8), allocatable :: work(:)
   double precision, allocatable :: rwork(:)
   integer, allocatable :: iwork(:)
   double precision,dimension(:),allocatable  :: evals_t
   complex(kind=8),dimension(:,:),allocatable :: evecs_r
   complex(kind=8), allocatable :: Hortho(:,:),Hortho1(:,:)
  
   ! Construct overlap matrices and diagonalize it.
  
   ! lwork = 8*ham%l
    !allocate(ifail(ham%l))
     !z = sys%nex*(sys%nex+1)/2
     !z = sys%ncutoff
     ham%l = sys%neb*sys%neb*sys%nQ 
     allocate(ham%O(ham%l,ham%l))
     allocate(ham%A(ham%l,ham%l))
     ham%A = cmplx(0.0,0.0)
 
    ham%O = cmplx(0.0,0.0)
    !ham%lambda_tttt = cmplx(0.0,0.0)
     ! ham%xi_d_tttt = cmplx(0.0,0.0)
     ! ham%xi_d_in_tttt = cmplx(0.0,0.0)
    c1 =0
    c2 =0
    c3 =0
    c4 = 0

    
    
   
    !do c1 = 1,sys%ncutoff
    do iQ_l =1,sys%nQ
      mQ_l = -sys%Qpts(:,iQ_l)
      !print*, "Q_l",iQ_l,Q_l(1),Q_l(2),Q_l(3)
      call Qpoint_to_index(mQ_l,imQ_l)
      do M_s = 1,sys%neb
         do N_s = 1,sys%neb
            c1 = c1 + 1
            c2 = 0
           do iQ_r =1,sys%nQ
              mQ_r = -sys%Qpts(:,iQ_r)
              call Qpoint_to_index(mQ_r,imQ_r)
               do I_s = 1,sys%neb
                 do J_s = 1,sys%neb
                    c2 = c2 + 1
                        !if (c1 .eq. c2)then
                           if(I_s == M_s .and. J_s == N_s .and. iQ_l==iQ_r .and. imQ_l==imQ_r)then 
                              ham%O(c1,c2) = ham%O(c1,c2) + cmplx(1,0)    
                              ham%A(c1,c2) = ham%A(c1,c2) + cmplx(exciton_sys%eigenvalues_t(I_s,iQ_r) + exciton_sys%eigenvalues_t(J_s,imQ_r),0.0)
                           end if

                           if(I_s == N_s .and. J_s == M_s .and. iQ_r==imQ_l .and. imQ_r==iQ_l)then
                              ham%O(c1,c2) = ham%O(c1,c2) + cmplx(1,0)
                              ham%A(c1,c2) = ham%A(c1,c2) + cmplx(exciton_sys%eigenvalues_t(I_s,iQ_r) + exciton_sys%eigenvalues_t(J_s,imQ_r), 0.0)
                           end if
                         
                        !end if
                        ham%O(c1,c2) = ham%O(c1,c2)  - ham%lambda_tttt(I_s,J_s,M_s,N_s,iQ_r,iQ_l) -  ham%lambda_tttt(I_s,J_s,N_s,M_s,iQ_r,imQ_l)
                        ham%A(c1,c2) = ham%A(c1,c2)+ ham%xi_d_tttt(I_s,J_s,M_s,N_s,iQ_r,iQ_l) + ham%xi_d_tttt(I_s,J_s,N_s,M_s,iQ_r,imQ_l) &
                                      -ham%xi_d_in_tttt(I_s,J_s,M_s,N_s,iQ_r,iQ_l) - ham%xi_d_in_tttt(I_s,J_s,N_s,M_s,iQ_r,imQ_l) &
                                       -((cmplx(exciton_sys%eigenvalues_t(I_s,iQ_r) + exciton_sys%eigenvalues_t(J_s,imQ_r),0.0))*(ham%lambda_tttt(I_s,J_s,M_s,N_s,iQ_r,iQ_l) + ham%lambda_tttt(I_s,J_s,N_s,M_s,iQ_r,imQ_l)))
                        print*,c1,c2,ham%O(c1,c2)
                        print*,c1,c2,ham%A(c1,c2),(exciton_sys%eigenvalues_t(I_s,iQ_r) + exciton_sys%eigenvalues_t(J_s,imQ_r))
                        
                  end do
               end do
            end do
         end do
      end do
   end do
   sum = cmplx(0.0,0.0)
   IJ_amp = 0.0
   A_MN = 0.0
   A_MNe = 0.0
   A_MNh = 0.0
   A_NM = 0.0
   !I = 4 
   !J =1
   
   ik1 = 1
   ik2 = 2
   !iQ_r = 2
  ! do ik1 = 1,sys%nk
  ! do ik2 = 1,sys%nk
  ! do I = 1,sys%neb 
  !    do J = 1,sys%neb
  !       do iQ_r = 1,sys%nQ
  !       sum = cmplx(0.0,0.0)
!
  ! 
  !    do M = 1,sys%neb
  !      do N = 1,sys%neb
  !       do iQ_l = 1,sys%nQ
  !          A_MN = 0.0
  !          A_MNe = 0.0
  !          A_MNh = 0.0
  !          A_NM = 0.0
!
  !          mQl = -sys%Qpts(:,iQ_l)
  !          call Qpoint_to_index(mQl,imQ_l)
  !          if(iQ_l == iQ_r)then
  !             A_MN = (exciton_sys%eigenvectors_t(1,1,1,ik1,M,iQ_l))*(exciton_sys%eigenvectors_t(1,1,1,ik2,N,imQ_l))
  !          end if
  !          p = sys%kpts(:,ik2) - sys%kpts(:,ik1) + sys%Qpts(:,iQ_r)
  !          call kpoint_to_index(p,ip)
  !          if(iQ_l == ip)then
  !             A_MNe = (exciton_sys%eigenvectors_t(1,1,1,ik2,M,iQ_l))*(exciton_sys%eigenvectors_t(1,1,1,ik1,N,imQ_l))
  !          end if
  !          q = sys%kpts(:,ik1) - sys%kpts(:,ik2) - sys%Qpts(:,iQ_r)
  !          call kpoint_to_index(q,iq)
  !          if(iQ_l==iq)then
  !             A_MNh = (exciton_sys%eigenvectors_t(1,1,1,ik1,M,iQ_l))*(exciton_sys%eigenvectors_t(1,1,1,ik2,N,imQ_l))
  !          end if
  !          if(imQ_l == iQ_r) then 
  !             A_NM =(exciton_sys%eigenvectors_t(1,1,1,ik2,M,iQ_l))*(exciton_sys%eigenvectors_t(1,1,1,ik1,N,imQ_l))
  !          end if
  !          if(I ==M .and.J==N .and. iQ_l == iQ_r)then
  !             IJ_amp = (A_MN - A_MNe - A_MNh + A_NM)
  !          end if
  !          sum = sum + ham%lambda_tttt(I,J,M,N,iQ_r,iQ_l)*(A_MN - A_MNe -A_MNh + A_NM)
  !          !- exciton_sys%eigenvectors(1,1,1,ik1,M,iQ_l)*exciton_sys%eigenvectors(1,1,1,ik2,M,imQ_l) &
  !          !- exciton_sys%eigenvectors(1,1,1,ik1,M,iQ_r)*exciton_sys%eigenvectors(1,1,1,ik2,N,iQ_l) &
  !          !+ exciton_sys%eigenvectors(1,1,1,ik1,N,iQ_r)*exciton_sys%eigenvectors(1,1,1,ik2,M,iQ_l))
  !        
  !       end do
  !    end do
  ! end do      
  !
  ! print*,"projection on IJ",sum,IJ_amp,abs(sum+IJ_amp)         
     
  !end do
  !end do
  !end do
  !end do
  !end do
      
    z = 64 !debugging stattemen
   ilow = 1
   iup = z
   abstol = 0.00000001
!
  !  print*,H(1,1),"this is me"
  ! print*,l,p,q,s
   allocate(evals_t(z))
   allocate(evecs_r(z,z)) !500
   allocate(work(20*z))
   allocate(iwork(10*z))
   allocate(ifail(z))
   allocate(rwork(14*z))
  !  print*,l,p,q,s,z
    
   call zheevx('V','I','U',z,ham%O,z,ellow,elup,ilow,iup, &
                          abstol,nfound,evals_t,evecs_r,z,work,8*z,rwork,iwork,ifail,info)

   c =0
   allocate(ham%c(ham%l,z))
   ham%c = cmplx(0.0,0.0)
   do r=1,z !500
      print*, evals_t(r)*13.6056980659
      if(r .ge. 55) then
         c= c+1
         ham%c(:,c) = evecs_r(:,r)/sqrt(evals_t(r))
      end if 
   end do 

   deallocate(evals_t)
   deallocate(evecs_r)
   deallocate(work)
   deallocate(iwork)
   deallocate(rwork)
   deallocate(ifail)

   z = 10
   iup = z
   allocate(Hortho1(z,z))
   Hortho1 = cmplx(0.0,0.0)
   do v1 = 1,10
      do v2 = 1,10
         do c1 = 1,64
            do c2 = 1,64
            Hortho1(v1,v2) = Hortho1(v1,v2) +conjg(ham%c(c1,v1))*(ham%c(c2,v2))*ham%A(c1,c2)
            end do
         end do
         print*,v1,v2,Hortho1(v1,v2)
      end do
   end do
   allocate(Hortho(z,z))
   Hortho = cmplx(0.0,0.0)
   Hortho = matmul(conjg(transpose(ham%c(:,1:z))),matmul(ham%A,ham%c(:,1:z)))
   do v1 = 1,z
      do v2 = 1,z
         print*,v1,v2,Hortho(v1,v2)!Hortho1(v1,v2)
      end do
   end do
   !z = 10
   iup = z
   allocate(evals_t(z))
   allocate(evecs_r(z,z)) !500
   allocate(work(20*z))
   allocate(iwork(10*z))
   allocate(ifail(z))
   allocate(rwork(14*z))
   call zheevx('V','I','U',z,Hortho,z,ellow,elup,ilow,iup, &
                          abstol,nfound,evals_t,evecs_r,z,work,8*z,rwork,iwork,ifail,info)
   

   do r=1,z
      print*,  evals_t(r)*13.6056980659
   end do
   
!
 !
    
            
   end subroutine solve_eigenvalue_overlap

   subroutine solve_eigenvalue_H
   ! diagonalize Hamiltonian ham%h using lapck routine dsyevx
      integer :: z
      complex(kind=8),allocatable,dimension(:) :: work
      integer,allocatable,dimension(:) :: iwork
      integer,allocatable,dimension(:) :: ifail
      integer :: lwork
      integer :: info,ilo,ihi
      integer :: nfound
      integer :: c1,i
      double precision :: abstol,ellow,elup,abnrm,bbnrm
      complex(kind=8),allocatable,dimension(:,:) :: vr
      double precision,allocatable,dimension(:) :: alphar
      double precision,allocatable,dimension(:) :: alphai
      double precision,allocatable,dimension(:) :: beta
      double precision,allocatable,dimension(:,:) :: vl
      double precision,allocatable,dimension(:) :: scale
      double precision,allocatable,dimension(:) :: rconde
      double precision,allocatable,dimension(:) :: rcondv
      double precision, allocatable,dimension(:) :: rwork
      logical,allocatable,dimension(:) :: bwork
      double precision,allocatable,dimension(:) :: lscale,rscale
      sys%num_ortho = ham%l
      abstol = -1.0
      allocate(rwork(7*sys%num_ortho))
      allocate(work(20*sys%num_ortho))
      allocate(iwork(5*sys%num_ortho))
      allocate(ifail(sys%num_ortho))
      allocate(vr(sys%num_ortho,sys%num_ortho))
      allocate(alphar(sys%num_ortho))
      allocate(alphai(sys%num_ortho))
      allocate(vl(sys%num_ortho,sys%num_ortho))
      allocate(scale(sys%num_ortho))
      allocate(rconde(sys%num_ortho))
      allocate(rcondv(sys%num_ortho))
      allocate(bwork(sys%num_ortho))
      allocate(beta(sys%num_ortho)) 
      allocate(lscale(sys%num_ortho))
      allocate(rscale(sys%num_ortho))
      lwork = 20*sys%num_ortho
      z = sys%num_ortho
      print*,"sys%num_ortho",sys%num_ortho
      !call zheevx('V','I','U',z,ham%O,z,ellow,elup,1,z, &
      !         abstol,nfound,alphar,vr,z,work,lwork,rwork,iwork,ifail,info)
      !call zgeevx('N', 'N', 'V', 'N', z, ham%H, z, &
    !alphar, alphai, vl, z, vr, z, &
    !ilo,ihi, RSCALE, bbnrm, &
    !rconde, rcondv, work, lwork, iwork, info)
      
      do i = 1,sys%num_ortho
                  !print*,alphar(i),alphai(i),beta(i)
        print*,i,"th eigenvector"
         print*,(alphar(i)*27.2114)
       
         !beta = matmul(ham%c,vr(:,i))
         !abnrm = norm2(beta) 
         !print*,beta/abnrm
         !print*,vr(:,i)
      end do
   end subroutine solve_eigenvalue_H

end module solve_generalized_eigenvalue_problem
