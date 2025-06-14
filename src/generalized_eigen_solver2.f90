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
    integer,allocatable,dimension(:) :: iwork
    logical,allocatable,dimension(:) :: bwork
    integer :: lwork,info,nfound
    integer,dimension(:),allocatable :: ifail
    integer :: i,j,alpha,bet,z,k,s,a,b
    double precision :: abnrm,bbnrm,abstol,norm
    integer :: I_ex,J_ex,M_ex,N_ex,n,m
    integer :: c1,c2,c3,c4,c
    double precision :: ellow,elup
    complex(kind=8) :: sum,sum2,sum3,sum4
    complex(kind=8) :: A_MN,A_MNe,A_MNh,A_NM,IJ_amp
    integer :: ik1,ik2,iq
    double precision :: q(3),mQl(3)
    double precision,allocatable,dimension(:,:) :: cMN,cijab
    double precision,allocatable,dimension(:,:) :: old
    integer,allocatable :: iarray(:),jarray(:),aarray(:),barray(:)
    double precision , allocatable :: c_ijab(:)
    CHARACTER(100) :: line,res1
    CHARACTER, DIMENSION(100) :: res
    integer :: status
    integer :: p1,p2,ind,v
    double precision :: x,c_v_MN
   integer :: iQ_l,iQ_r,imQ_l
   double precision  ::Q_l(3)
   ! Construct overlap matrices and diagonalize it.
  
   ! lwork = 8*ham%l
    !allocate(ifail(ham%l))
     !z = sys%nex*(sys%nex+1)/2
     !z = sys%ncutoff
     ham%l = sys%neb*sys%neb*sys%nQ 
     allocate(ham%O(ham%l,ham%l))
 
    ham%O = cmplx(0.0,0.0)
    c1 =0
    c2 =0
    c3 =0
    c4 = 0

    
    
   
    !do c1 = 1,sys%ncutoff
    do iQ_l =1,sys%nQ
      Q_l = -sys%Qpts(:,iQ_l)
      print*, "Q_l",iQ_l,Q_l(1),Q_l(2),Q_l(3)
      call Qpoint_to_index(Q_l,imQ_l)
      do M_s = 1,sys%neb
         do N_s = 1,sys%neb
            c1 = c1 + 1
            c2 = 0
           do iQ_r =1,sys%nQ
               do I_s = 1,sys%neb
                 do J_s = 1,sys%neb
                    c2 = c2 + 1
                        if (c1 .eq. c2)then
                           !print*,"c1,c2",c1,c2
                           ham%O(c1,c2) = ham%O(c1,c2) + 2.0
                        end if
                        ham%O(c1,c2) = ham%O(c1,c2)  - ham%lambda_tttt(I_s,J_s,M_s,N_s,iQ_r,iQ_l) -  ham%lambda_tttt(I_s,J_s,N_s,M_s,iQ_r,imQ_l)
                        print*,c1,c2,ham%O(c1,c2)
                        
                  end do
               end do
            end do
         end do
      end do
   end do
   sum = 0.0
   IJ_amp = 0.0
   A_MN = 0.0
   A_MNe = 0.0
   A_MNh = 0.0
   A_NM = 0.0
   I = 1 
   J =4
   
   ik1 = 1
   ik2 = 2
   iQ_r = 4

   
   do M = 1,sys%neb
      do N = 1,sys%neb
         do iQ_l = 1,sys%nQ

            mQl = -sys%Qpts(:,iQ_l)
            call Qpoint_to_index(mQl,imQ_l)
            A_MN = (exciton_sys%eigenvectors_t(1,1,1,ik1,M,iQ_l))*(exciton_sys%eigenvectors_t(1,1,1,ik2,N,imQ_l))
            if(ik1 == ik2)then
               A_MNe = (exciton_sys%eigenvectors_t(1,1,1,ik2,M,iQ_l))*(exciton_sys%eigenvectors_t(1,1,1,ik1,N,imQ_l))
            end if
            q = sys%kpts(:,ik1) - sys%kpts(:,ik2) - 2 *sys%Qpts(:,iQ_l)
            call kpoint_to_index(q,iq)
            if(iq==1)then
               A_MNh = (exciton_sys%eigenvectors_t(1,1,1,ik1,M,iQ_l))*(exciton_sys%eigenvectors_t(1,1,1,ik2,N,imQ_l))
            end if

            A_NM =(exciton_sys%eigenvectors_t(1,1,1,ik2,M,iQ_l))*(exciton_sys%eigenvectors_t(1,1,1,ik1,N,imQ_l))
            if(I ==M .and.J==N .and. iQ_l == iQ_r)then
               IJ_amp = (A_MN - A_MNe - A_MNh + A_NM)
            end if
            sum = sum + ham%lambda_tttt(I,J,M,N,iQ_r,iQ_l)*(A_MN - A_MNe -A_MNh + A_NM)
            !- exciton_sys%eigenvectors(1,1,1,ik1,M,iQ_l)*exciton_sys%eigenvectors(1,1,1,ik2,M,imQ_l) &
            !- exciton_sys%eigenvectors(1,1,1,ik1,M,iQ_r)*exciton_sys%eigenvectors(1,1,1,ik2,N,iQ_l) &
            !+ exciton_sys%eigenvectors(1,1,1,ik1,N,iQ_r)*exciton_sys%eigenvectors(1,1,1,ik2,M,iQ_l))
         end do
      end do
   end do      
  
   print*,"projection on IJ",sum,IJ_amp         
            
      
    !ham%A =- ham%B
    z = ham%l
   ! call dsyevx('V','I','U',z,ham%O,z,ellow,elup,1,z, &
    !                          -1,nfound,alphar,vr,z,work,8*z,iwork,ifail,info)
   !print*,ifail
   !print*,alphar
   ! do gram schmidt for IJ operators
    !allocate(cMN(ham%l,ham%l))
 !   sys%num_ortho = 3
 !   allocate(ham%c(sys%ncutoff,sys%num_ortho))
 !   allocate(old(sys%ncutoff,sys%num_ortho))
 !   ham%c = 0.0
 !   old = 0.0
 !   do i = 1,sys%ncutoff
 !      do j = 1,sys%num_ortho
 !         if(i == j)then
 !        
 !            old(i,j) = 1.0
 !         end if
 !      end do
 !  end do
!
 !  !gram-schmidt
!
 !    do j = 1,sys%num_ortho
 !      norm = dot_product(old(:,j),matmul(ham%O,old(:,j))) 
 !      print*,norm
 !      ham%c(:,j) = old(:,j)/sqrt(norm)
 !      do k = j+1,sys%num_ortho
 !          old(:,k) =  old(:,k) - (dot_product(ham%c(:,j),matmul(ham%O,old(:,k))) * ham%c(:,j))
 !      end do
 !   end do 
 !!  s = sys%nc*(sys%nc - 1)*sys%nv*(sys%nv - 1)/4
 !!  allocate(c_ijab(s))
 !!  allocate(iarray(s))
 !!  allocate(jarray(s))
 !!  allocate(aarray(s))
 !!  allocate(barray(s))
 !!  cijab = 0.0
 !!  iarray = 0
 !!  jarray = 0
 !!  aarray = 0
 !!  barray = 0
!!
 !  ! calculate ijab indices for each entry
!
!   c = 0
!   do j = 1,sys%nc-1
!      do i = j+1,sys%nc
!          do b = 1,sys%nv-1
!              do a = b+1,sys%nv
!                  c = c+1
!
!               iarray(c) = i
!               jarray(c) = j
!               aarray(c) = a
!               barray(c) = b
!              end do
!          end do
!      end do
!   end do
!   ! read the biexciton eigenvector in 4p basis
!    open (1, file='line.dat', status='old')
!     do v = 1,1
!        print*,"vth eigenvector",v
!        do ind = 1,s
!           READ(1, '(A)', IOSTAT=STATUS) line
!           res = TRANSFER(line,res)
!           res1 = TRANSFER(PACK(res,res/=' '),res1)
!           !print*,res1
!           !READ(trim(res1))char1,char2,int1,char3,int2,char4,char5,x,char6,y,char7
!           p1 = index(res1, '=') +1
!           p2 = index(res1, '+') -1 
!           !print*,p1,p2
!           read(res1(p1:p2), *, iostat=STATUS) x
!           if (STATUS /= 0) stop "error reading from internal list-directed read"
!           !if(iarray(ind) .le. 10 .and. jarray(ind) .le. 10 .and. aarray(ind) .le. 10 .and. barray(ind) .le. 10)then
!            ! print*,iarray(ind),jarray(ind),aarray(ind),barray(ind),x
!          ! end if
!            c_ijab(ind) = x
!     
!        end do
!     end do
!
!     close(1)
!



   ! do i = 1,sys%num_ortho
       !do j = 1,sys%num_ortho
           !print*,j,"jth eigen vector"
          ! print*,i,j,dot_product(ham%c(:,i),matmul(ham%O,ham%c(:,j)))
    !      if (i == 32) then
     !         print*,ham%c(:,i)
     !     end if
          ! print*,ham%c(:,j)
      !end do
   !end do

    !read the vector in 4p
    ! c1 = 0
    ! c3 = 0
    !  sum = 0.0
    !  sum4 = 0.0
    !  do i = 1,sys%num_ortho
    !     sum3 = 0.0
    !     do c1 = 1,sys%ncutoff
    !           M_s = exciton_sys%I_arr(c1)  
    !           N_s = exciton_sys%J_arr(c1)
    !           
    !        
    !           !if (M_s .ge. N_s)then
    !              sum2 = 0.0
    !              !c3 = c3 + 1
    !               !c_v_MN = sum2 
    !               do ind = 1,s
    !                  sum2 = sum2 +  c_ijab(ind)*(exciton_sys%eigenvectors(M_s,aarray(ind),iarray(ind))*exciton_sys%eigenvectors(N_s,barray(ind),jarray(ind)) &
    !                        - exciton_sys%eigenvectors(M_s,barray(ind),iarray(ind))*exciton_sys%eigenvectors(N_s,aarray(ind),jarray(ind))&
    !                        - exciton_sys%eigenvectors(M_s,aarray(ind),jarray(ind))*exciton_sys%eigenvectors(N_s,barray(ind),iarray(ind)) &
    !                        + exciton_sys%eigenvectors(M_s,barray(ind),jarray(ind))*exciton_sys%eigenvectors(N_s,aarray(ind),iarray(ind)))
    !              end do
    !              c_v_MN = sum2*ham%c(c1,i)                
    !              sum3  = sum3 + c_v_MN
    !           !end if
    !        
    !     end do  
    !     print*,"coefficient of",i,"th vector",sum3 
    !     sum4 = sum4 + (sum3)**2
    !  end do
    !   
    !  print*,"total sum of squares value",sum4
   



   ! call dggevx('N', 'N', 'V', 'N', ham%l, ham%A, ham%l, ham%B, ham%l, &
   ! !alphar, alphai, beta, vl, ham%l, vr, ham%l, &
   ! !ilo,ihi, LSCALE, RSCALE, abnrm, bbnrm, &
   ! !rconde, rcondv, work, lwork, iwork, bwork, info)
   ! 
   ! !allocate(cMN(ham%l,ham%l))
   ! !allocate(cijab(ham%l,ham%l))
!
   !! c_MN = 0.0
   ! !c_ijab = 0.0
   ! !c1 = 0
   ! !c2 = 0
   ! 
   !allocate(ham%c(ham%l,sys%northo))
   !ham%c = 0.0
   ! 
   ! !convert y_mn to c_mn
   !c1 =0
   ! do i = 1,ham%l
   !   !print*,alphar(i),alphai(i),beta(i)
   !   print*,i,"th eigenvector"
   !   print*,(alphar(i))
   !   if(i .ge. 781)then
   !      c1 = c1 +1
   !      ham%c(:,c1) = vr(:,i)
   !      !print*,dot_product(ham%c(:,c1),vr(:,36))
   !   end if 
   !  !print*,vr(:,i)
   !
   ! end do
   !c1 = 0
   ! !do i = 1,sys%northo
   ! !   c1 = 0
   ! !   do I_ex = 1,sys%nex
    !       do J_ex = 1,sys%nex
    !          c1 = c1 + 1
    !           if (I_ex .le. 20 .and. J_ex .le. 20)then
    !              if (ham%c(c1,i) .ge. 0.01)then
    !                  !if(I_ex .le. 3 .and. J_ex .le. 3)then
    !                    print *,"i,I,J,c",i,I_ex,J_ex,ham%c(c1,i)
    !                  !end if
    !              end if
    !           end if
    !        
    !       end do
    !   end do
    !end do
   
   !deallocate(ham%O,alphar,alphai,beta,vl,vr,lscale,rscale,rconde,rcondv,work,iwork,bwork,ifail)
   

 
    
            
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
      call zheevx('V','I','U',z,ham%O,z,ellow,elup,1,z, &
               abstol,nfound,alphar,vr,z,work,lwork,rwork,iwork,ifail,info)
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
