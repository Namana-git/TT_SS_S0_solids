module solve_generalized_eigenvalue_problem
   use global_variables
   use compute_lambda_and_xi
   use compute_AB

   !use lapacks generalized eigenvalue solver dggevx, with A = ham%A and B = ham%B and m=ham%l

   implicit none


   contains

      subroutine solve_eigenvalue_overlap()
    
    integer :: ilo,ihi,valid_count
    integer :: I_s,J_s,M_s,N_s
    !double precision :: ellow,elup
    double precision,allocatable,dimension(:) :: alphar,alphai,beta
    double precision,allocatable,dimension(:,:) :: vl,vr
    double precision,allocatable,dimension(:) :: lscale,rscale
    double precision,allocatable,dimension(:) :: rconde,rcondv
    double precision,allocatable,dimension(:) :: work
    integer,allocatable,dimension(:) :: iwork
    logical,allocatable,dimension(:) :: bwork
    integer :: lwork,info,nfound
    integer,dimension(:),allocatable :: ifail
    integer :: i,j,alpha,bet,z,k,s,a,b
    double precision :: abnrm,bbnrm,abstol,norm
    integer :: I_ex,J_ex,M_ex,N_ex,n,m
    integer :: c1,c2,c3,c4,c
    double precision :: sum,ellow,elup,sum1,sum2,sum3,sum4
    double precision,allocatable,dimension(:,:) :: cMN,cijab
    double precision,allocatable,dimension(:,:) :: old,new
    integer,allocatable :: iarray(:),jarray(:),aarray(:),barray(:)
    double precision , allocatable :: c_ijab(:)
    !CHARACTER(100) :: line,res1
    CHARACTER, DIMENSION(100) :: res
    integer :: status
    integer :: s1,s2
    !nteger :: p1,p2,ind,v
    !double precision :: x,c_v_MN
   ! Construct overlap matrices and diagonalize it.
  
    lwork = 8*ham%l
    print*,sys%ncutoff
    !sys%ncutoff =  sys%nex*(sys%nex+1)
    !sys%ncutoff = !sys%nex*(sys%nex+1)/2
     z = sys%ncutoff
     ham%l = sys%ncutoff
     print*,"hello7"
     allocate(ham%O(sys%ncutoff,sys%ncutoff))
     allocate(ham%H(sys%ncutoff,sys%ncutoff))
    print*,"Hello8"
   
    allocate(alphar(ham%l),alphai(ham%l),beta(ham%l))
    allocate(vl(ham%l,ham%l),vr(ham%l,ham%l))
    allocate(lscale(ham%l),rscale(ham%l))
    allocate(rconde(ham%l),rcondv(ham%l))
    allocate(work(lwork))
    allocate(iwork(ham%l+6))
    allocate(bwork(ham%l))
    allocate(ifail(ham%l))
    alphar = 0.0
    alphai = 0.0
    beta = 0.0
    vl = 0.0
    vr = 0.0
    BWORK = .false.
    lwork = 8*ham%l
   
    !ham%A = 0.0
    ham%O = 0.0
    ham%H =0.0
    c1 =0
    c2 =0
    c3 =0
    c4 = 0

    
    
   
    do c1 = 1,sys%ncutoff
      M_s = exciton_sys%I_arr(c1)  
      N_s = exciton_sys%J_arr(c1)
      s1 = exciton_sys%spin_ind(c1)
          print*,"c1,M_s,N_s,s1",c1,M_s,N_s,s1
      do c2 = 1,sys%ncutoff
         I_s = exciton_sys%I_arr(c2)  
         J_s = exciton_sys%J_arr(c2)
         s2 = exciton_sys%spin_ind(c2)
            

                 if (s1 == 1 .and. s2 ==1 )then
                    ham%O(c1,c2) = ham%O(c1,c2) + ham%B_tttt(I_s,J_s,M_s,N_s)!/sqrt(ham%B_tttt(I_s,J_s,I_s,J_s)*ham%B_tttt(M_s,N_s,M_s,N_s))
                    ham%H(c1,c2) = ham%H(c1,c2) + ham%A_tttt(I_s,J_s,M_s,N_s)!/sqrt(ham%B_tttt(I_s,J_s,I_s,J_s)*ham%B_tttt(M_s,N_s,M_s,N_s))
                     !print*,ham%H(c1,c2),ham%O(c1,c2),c1,c2,I_s,J_s,M_s,N_s
                     !print*,ham%O(c1,c2),c1,c2,s1,s2,ham%B_tttt(I_s,J_s,M_s,N_s),I_s,J_s,M_s,N_s
                   ! print*,ham%H(c1,c2),ham%O(c1,c2),c1,c2,I_s,J_s,M_s,N_s
                     !print*,ham%O(c1,c2),c1,c2,s1,s2,ham%B_tttt(I_s,J_s,M_s,N_s),I_s,J_s,M_s,N_s
                 end if
                 if (s1 == 2 .and. s2 ==2 )then
                    ham%O(c1,c2) = ham%O(c1,c2) + ham%B_ssss(I_s,J_s,M_s,N_s)!/sqrt(ham%B_ssss(I_s,J_s,I_s,J_s)*ham%B_ssss(M_s,N_s,M_s,N_s))
                    ham%H(c1,c2) = ham%H(c1,c2) + ham%A_ssss(I_s,J_s,M_s,N_s)!/sqrt(ham%B_ssss(I_s,J_s,I_s,J_s)*ham%B_ssss(M_s,N_s,M_s,N_s))
                     !print*,ham%H(c1,c2),ham%O(c1,c2),c1,c2,I_s,J_s,M_s,N_s
                 end if
                 if (s1 == 1 .and. s2==2 )then
                     ham%O(c1,c2) = ham%O(c1,c2) + ham%B_ttss(I_s,J_s,M_s,N_s)!/sqrt(ham%B_ssss(I_s,J_s,I_s,J_s)*ham%B_tttt(M_s,N_s,M_s,N_s))
                     ham%H(c1,c2) = ham%H(c1,c2) + ham%A_ttss(I_s,J_s,M_s,N_s)!/sqrt(ham%B_ssss(I_s,J_s,I_s,J_s)*ham%B_tttt(M_s,N_s,M_s,N_s))
                 end if
                 if (s1 == 2 .and. s2 ==1) then 
                     ham%O(c1,c2) = ham%O(c1,c2) + ham%B_sstt(I_s,J_s,M_s,N_s)!/sqrt(ham%B_tttt(I_s,J_s,I_s,J_s)*ham%B_ssss(M_s,N_s,M_s,N_s))
                     ham%H(c1,c2) = ham%H(c1,c2) + ham%A_sstt(I_s,J_s,M_s,N_s)!/sqrt(ham%B_tttt(I_s,J_s,I_s,J_s)*ham%B_ssss(M_s,N_s,M_s,N_s))
                     !print*,"heeeeee",ham%H(c1,c2),ham%O(c1,c2),c1,c2,I_s,J_s,M_s,N_s
                 end if
                 
                 !write(2,*)c1,c2,ham%O(c1,c2),ham%H(c1,c2)
               !end do
              !end do
            !end do
         !end do
       end do
     end do
    ! close(2)
   

   ! do gram schmidt for IJ operators
    !allocate(cMN(ham%l,ham%l))
   ! sys%num_ortho = 32
     !sys%ncutoff = sys%nex*(sys%nex+1)
      sys%num_ortho = 12
    allocate(ham%c(sys%ncutoff,sys%num_ortho))
    allocate(new(sys%ncutoff,sys%num_ortho))
    allocate(old(sys%ncutoff,sys%num_ortho))
    ham%c = 0.0
    old = 0.0
    do i = 1,sys%ncutoff
         do j = 1,sys%num_ortho
            if(i == j)then
           
               old(i,j) = 1.0
            end if
         end do
     end do
!!
!!
!
 !  !gram-schmidt
     valid_count = 0
     do j = 1,sys%num_ortho
       norm = dot_product(old(:,j),matmul(ham%O,old(:,j))) 
       print*,"norm",norm
       new(:,j) = old(:,j)/sqrt(norm)
     
       do k = j+1,sys%num_ortho
           old(:,k) =  old(:,k) - (dot_product(new(:,j),matmul(ham%O,old(:,k))) * new(:,j))
       end do
       if (norm > 0.0001) then
         valid_count = valid_count + 1
         ham%c(:, valid_count) = old(:,j) / sqrt(norm)
       else
            print *, "Vector ", j, " discarded (linearly dependent)"
       end if

     !  norm = dot_product(old(:,j), matmul(ham%O, old(:,j)))
      !print *, "Vector ", j, " norm after ortho = ", norm
!! 
      ! if (norm > 0.0001) then
      !   valid_count = valid_count + 1
      !   ham%c(:, valid_count) = old(:,j) / sqrt(norm)
      ! else
      !    print *, "Vector ", j, " discarded (linearly dependent)"
      ! end if
     end do 

     sys%lin_ind = valid_count
  ! valid_count=0
  !  do j = 1, sys%num_ortho
  !     !Orthogonalize against previous accepted vectors
  !    do k = 1, valid_count
  !       old(:,j) = old(:,j) - dot_product(ham%c(:,k), matmul(ham%O, old(:,j))) * ham%c(:,k)
  !   end do
!
  !   ! Check norm after orthogonalization
  !   norm = dot_product(old(:,j), matmul(ham%O, old(:,j)))
  !   print *, "Vector ", j, " norm after ortho = ", norm
!
  !   if (norm > 0.84) then
  !       valid_count = valid_count + 1
  !       ham%c(:, valid_count) = old(:,j) / sqrt(norm)
  !   else
  !       print *, "Vector ", j, " discarded (linearly dependent)"
  !   end if
  ! end do
  !  !print*,ham%c(1,:)
  !  sys%num_ortho = valid_count
  !  print*,valid_count
!
  !  
 ! ! s = sys%nc*(sys%nc - 1)*sys%nv*(sys%nv - 1)/4
 !  allocate(c_ijab(s))
 !  allocate(iarray(s))
 !  allocate(jarray(s))
 !  allocate(aarray(s))
 !  allocate(barray(s))
 !  cijab = 0.0
 !  iarray = 0
 !  jarray = 0
 !  aarray = 0
 !  barray = 0
!
   ! calculate ijab indices for each entry

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
!   end doG
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
   


   subroutine solve_eigenvalue_overlap1()
    
    integer :: ilo,ihi,valid_count
    integer :: I_s,J_s,M_s,N_s
    !double precision :: ellow,elup
    double precision,allocatable,dimension(:) :: alphar,alphai,beta
    double precision,allocatable,dimension(:,:) :: vl,vr
    double precision,allocatable,dimension(:) :: lscale,rscale
    double precision,allocatable,dimension(:) :: rconde,rcondv
    double precision,allocatable,dimension(:) :: work
    integer,allocatable,dimension(:) :: iwork
    logical,allocatable,dimension(:) :: bwork
    integer :: lwork,info,nfound
    integer,dimension(:),allocatable :: ifail
    integer :: i,j,alpha,bet,z,k,s,a,b
    double precision :: abnrm,bbnrm,abstol,norm
    integer :: I_ex,J_ex,M_ex,N_ex,n,m
    integer :: c1,c2,c3,c4,c
    double precision :: sum,ellow,elup,sum1,sum2,sum3,sum4
    double precision,allocatable,dimension(:,:) :: cMN,cijab
    double precision,allocatable,dimension(:,:) :: old
    integer,allocatable :: iarray(:),jarray(:),aarray(:),barray(:)
    double precision , allocatable :: c_ijab(:)
    !CHARACTER(100) :: line,res1
    CHARACTER, DIMENSION(100) :: res
    integer :: status
    integer :: s1,s2
    !nteger :: p1,p2,ind,v
    !double precision :: x,c_v_MN
   ! Construct overlap matrices and diagonalize it.
  
    lwork = 8*ham%l
    print*,sys%ncutoff
    !sys%ncutoff =  sys%nex*(sys%nex+1)
    !sys%ncutoff = !sys%nex*(sys%nex+1)/2
     z = sys%ncutoff
     ham%l = sys%ncutoff
     print*,"hello7"
     allocate(ham%O(sys%ncutoff,sys%ncutoff))
     allocate(ham%H(sys%ncutoff,sys%ncutoff))
    print*,"Hello8"
   
    allocate(alphar(ham%l),alphai(ham%l),beta(ham%l))
    allocate(vl(ham%l,ham%l),vr(ham%l,ham%l))
    allocate(lscale(ham%l),rscale(ham%l))
    allocate(rconde(ham%l),rcondv(ham%l))
    allocate(work(lwork))
    allocate(iwork(ham%l+6))
    allocate(bwork(ham%l))
    allocate(ifail(ham%l))
    alphar = 0.0
    alphai = 0.0
    beta = 0.0
    vl = 0.0
    vr = 0.0
    BWORK = .false.
    lwork = 8*ham%l
   
    !ham%A = 0.0
    ham%O = 0.0
    ham%H =0.0
    c1 =0
    c2 =0
    c3 =0
    c4 = 0

    
    
   
    !do c1 = 1,sys%ncutoff
      !M_s = exciton_sys%I_arr(c1)  
      !N_s = exciton_sys%J_arr(c1)
      !s1 = exciton_sys%spin_ind(c1)
    open(2,file='overlap.dat',status='replace')
   do s1 = 1,2
    do M_s = 1,sys%nex
      do N_s = M_s,sys%nex
         c1  = c1 + 1
         c2  = 0 
       do s2=1,2
         do I_s = 1,sys%nex
            do J_s = I_s,sys%nex
               c2 = c2 + 1
               !print*,"c1,c2",c1,c2,s1,s2

                 if (s1 == 1 .and. s2 ==1 )then
                    ham%O(c1,c2) = ham%O(c1,c2) + ham%B_tttt(I_s,J_s,M_s,N_s)!/sqrt(ham%B_tttt(I_s,J_s,I_s,J_s)*ham%B_tttt(M_s,N_s,M_s,N_s))
                    ham%H(c1,c2) = ham%H(c1,c2) + ham%A_tttt(I_s,J_s,M_s,N_s)!/sqrt(ham%B_tttt(I_s,J_s,I_s,J_s)*ham%B_tttt(M_s,N_s,M_s,N_s))
                     !print*,ham%H(c1,c2),ham%O(c1,c2),c1,c2,I_s,J_s,M_s,N_s
                     !print*,ham%O(c1,c2),c1,c2,s1,s2,ham%B_tttt(I_s,J_s,M_s,N_s),I_s,J_s,M_s,N_s
                   ! print*,ham%H(c1,c2),ham%O(c1,c2),c1,c2,I_s,J_s,M_s,N_s
                     !print*,ham%O(c1,c2),c1,c2,s1,s2,ham%B_tttt(I_s,J_s,M_s,N_s),I_s,J_s,M_s,N_s
                 end if
                 if (s1 == 2 .and. s2 ==2 )then
                    ham%O(c1,c2) = ham%O(c1,c2) + ham%B_ssss(I_s,J_s,M_s,N_s)!/sqrt(ham%B_ssss(I_s,J_s,I_s,J_s)*ham%B_ssss(M_s,N_s,M_s,N_s))
                    ham%H(c1,c2) = ham%H(c1,c2) + ham%A_ssss(I_s,J_s,M_s,N_s)!/sqrt(ham%B_ssss(I_s,J_s,I_s,J_s)*ham%B_ssss(M_s,N_s,M_s,N_s))
                     !print*,ham%H(c1,c2),ham%O(c1,c2),c1,c2,I_s,J_s,M_s,N_s
                 end if
                 if (s1 == 1 .and. s2==2 )then
                     ham%O(c1,c2) = ham%O(c1,c2) + ham%B_ttss(I_s,J_s,M_s,N_s)!/sqrt(ham%B_ssss(I_s,J_s,I_s,J_s)*ham%B_tttt(M_s,N_s,M_s,N_s))
                     ham%H(c1,c2) = ham%H(c1,c2) + ham%A_ttss(I_s,J_s,M_s,N_s)!/sqrt(ham%B_ssss(I_s,J_s,I_s,J_s)*ham%B_tttt(M_s,N_s,M_s,N_s))
                 end if
                 if (s1 == 2 .and. s2 ==1) then 
                     ham%O(c1,c2) = ham%O(c1,c2) + ham%B_sstt(I_s,J_s,M_s,N_s)!/sqrt(ham%B_tttt(I_s,J_s,I_s,J_s)*ham%B_ssss(M_s,N_s,M_s,N_s))
                     ham%H(c1,c2) = ham%H(c1,c2) + ham%A_sstt(I_s,J_s,M_s,N_s)!/sqrt(ham%B_tttt(I_s,J_s,I_s,J_s)*ham%B_ssss(M_s,N_s,M_s,N_s))
                     !print*,"heeeeee",ham%H(c1,c2),ham%O(c1,c2),c1,c2,I_s,J_s,M_s,N_s
                 end if
                 
                 write(2,*)c1,c2,ham%O(c1,c2),ham%H(c1,c2)
               end do
              end do
            end do
         end do
       end do
     end do
     close(2)
   end subroutine solve_eigenvalue_overlap1
   
   subroutine solve_eigenvalue_O
      ! diagonalize Hamiltonian ham%h using lapck routine dsyevx
         integer :: z
         double precision,allocatable,dimension(:) :: work
         integer,allocatable,dimension(:) :: iwork
         integer,allocatable,dimension(:) :: ifail
         integer :: lwork
         integer :: info,ilo,ihi
         integer :: nfound
         integer :: c1,i,c,j
         double precision :: abstol,ellow,elup,abnrm,bbnrm,norm,max_diff
         double precision,allocatable,dimension(:,:) :: vr
         double precision,allocatable,dimension(:) :: alphar
         double precision,allocatable,dimension(:) :: alphai
         double precision,allocatable,dimension(:) :: beta
         double precision,allocatable,dimension(:,:) :: vl
         double precision,allocatable,dimension(:) :: scale
         double precision,allocatable,dimension(:) :: rconde
         double precision,allocatable,dimension(:) :: rcondv
         logical,allocatable,dimension(:) :: bwork
         double precision,allocatable,dimension(:) :: lscale,rscale
         real,allocatable,dimension(:,:) :: temp_O
         integer, allocatable :: isuppz(:)
         abstol = -1
         sys%num_ortho = sys%ncutoff!sys%nex*(sys%nex+1)
         allocate(temp_O(sys%num_ortho,sys%num_ortho))
         temp_O = ham%O 
         allocate(work(30*sys%num_ortho))
         allocate(iwork(sys%num_ortho*12))
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
         allocate(isuppz(2*sys%num_ortho))
         lwork = 30*sys%num_ortho
      

         z = sys%ncutoff!sys%num_ortho
         !do i = 1, z
         !   ham%O(i,i) = ham%O(i,i) + 0.0000001
         !end do
         !print*,z
         !call dpotrf('U', z, ham%O, z, info)   ! LAPACK Cholesky
         !if (info /= 0) stop "Overlap matrix O is not positive-definite"

         max_diff = 0.0d0
         do i = 1, 20
          !do j = i+1, 20
            ham%O(i,i) = ham%O(i,i) + 0.00001
          !end do
         end do
         print *, "Max asymmetry in matrix: ", max_diff
         call dsyevr('V','A','U',z,ham%O,z,ellow,elup,1,z, &
                 abstol,nfound,alphar,vr,z,isuppz,work,lwork,iwork,sys%num_ortho*12,info)


         !call dgeevx('N', 'N', 'V', 'N', z, ham%O, z, &
       !alphar, alphai, vl, z, vr, z, &
       !ilo,ihi, RSCALE, bbnrm, &
       !rconde, rcondv, work, lwork, iwork, info)
         c =0
         do i = 1,sys%num_ortho
       !              !print*,alphar(i),alphai(i),beta(i)
           norm = dot_product(vr(:,i),matmul(temp_O,vr(:,i)))
           if(abs(norm) .ge. 0.0001)then 
              c = c+1
           end if
       !    print*,i,"th eigenvector"
       !    
            print*,(alphar(i))!*27.2114)
       !     print*,(norm)
           !print*,vr(:,i)
!
       !     
       !     !beta = matmul(ham%c,vr(:,i))
       !     !abnrm = norm2(beta) 
       !     !print*,beta/abnrm
       !     !print*,vr(:,i)
         end do
        ! z = sys%num_ortho
         sys%num_ortho = c
         z =c
         allocate(ham%c(sys%nex*(sys%nex+1),sys%num_ortho))
         c =0
         do i = 1,(sys%nex*(sys%nex+1))
            norm = dot_product(vr(:,i),matmul(temp_O,vr(:,i)))
            if(abs(alphar(i)) .ge. 0.001)then
              
              c = c+1
               
!
              ham%c(:,c) = vr(:,i)/(sqrt(alphar(i)))
              norm = dot_product(ham%c(:,c),matmul(temp_O,ham%c(:,c)))
              !print*,alphar(i)
              print*,norm
                 !print*,"c",c,ham%c(:,c)
            end if
         end do
                  
      end subroutine solve_eigenvalue_O

   subroutine solve_eigenvalue_H
   ! diagonalize Hamiltonian ham%h using lapck routine dsyevx
      integer :: z
      double precision,allocatable,dimension(:) :: work
      integer,allocatable,dimension(:) :: iwork
      integer,allocatable,dimension(:) :: ifail
      integer :: lwork
      integer :: info,ilo,ihi
      integer :: nfound
      integer :: c1,i
      double precision :: abstol,ellow,elup,abnrm,bbnrm
      double precision,allocatable,dimension(:,:) :: vr
      double precision,allocatable,dimension(:) :: alphar
      double precision,allocatable,dimension(:) :: alphai
      double precision,allocatable,dimension(:) :: beta
      double precision,allocatable,dimension(:,:) :: vl
      double precision,allocatable,dimension(:) :: scale
      double precision,allocatable,dimension(:) :: rconde
      double precision,allocatable,dimension(:) :: rcondv
      logical,allocatable,dimension(:) :: bwork
      double precision,allocatable,dimension(:) :: lscale,rscale
      double precision,allocatable,dimension(:,:) :: H1,O1
      abstol = -1
      sys%num_ortho = sys%lin_ind
      print*,"sys%num_ortho",sys%num_ortho
      !sys%num_ortho = sys%nex*(sys%nex+1)
      allocate(H1(sys%lin_ind,sys%lin_ind))
      H1 = 0.0
     ! allocate(O1(sys%num_ortho,sys%num_ortho))

         H1 = matmul(transpose(ham%c(:,1:sys%lin_ind)),matmul(ham%H,ham%c(:,1:sys%lin_ind)))
         print*,"Hihihi"
      !O1 = ham%O !+ transpose(ham%O)/2
     ! allocate(work(8*sys%num_ortho))
      allocate(iwork(sys%num_ortho+6))
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
      lwork = -1
      allocate(work(1))
      z = sys%num_ortho
      !print*,ham%O(1,1),ham%O(1,2),ham%O(2,1),ham%O(2,2)
      call dsyevx('V','I','U',z,H1,z,ellow,elup,1,z, &
              abstol,nfound,alphar,vr,z,work,lwork,iwork,ifail,info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork)) 

      call dsyevx('V','I','U',z,H1,z,ellow,elup,1,z, &
              abstol,nfound,alphar,vr,z,work,lwork,iwork,ifail,info)
      !call dsyevx('V','I','U',z,ham%O,z,ellow,elup,1,z, &
      !        abstol,nfound,alphar,vr,z,work,lwork,iwork,ifail,info)
      !print*,"ifail",ifail
   !   call dgeevx('N', 'N', 'V', 'N', z, ham%H, z, &
   ! alphar, alphai, vl, z, vr, z, &
   ! ilo,ihi, RSCALE, bbnrm, &
   ! rconde, rcondv, work, lwork, iwork, info)
     ! call dggev('N', 'V', z, ham%H, z, ham%O, z, alphar, alphai, beta, vl, z, vr, z, work, lwork, info)
     ! lwork = int(work(1))
     ! deallocate(work)
     ! allocate(work(lwork)) 
     ! call dggev('N', 'V', z, ham%H, z, ham%O, z, alphar, alphai, beta, vl, z, vr, z, work, lwork, info) 
     ! print*,info
      do i = 1,sys%num_ortho
                  !print*,alphar(i),alphai(i),beta(i)
        print*,i,"th eigenvector"
         print*,(alphar(i)*27.2114)
         !print*,(alphar(i)/beta(i)*27.2114)
         
         !beta = matmul(ham%c,vr(:,i))
         !abnrm = norm2(beta) 
         !print*,beta/abnrm
         !print*,vr(:,i)
      end do

   
   end subroutine solve_eigenvalue_H

end module solve_generalized_eigenvalue_problem
