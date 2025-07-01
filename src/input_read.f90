module input_read
  use global_variables
  use hdf5
  implicit none
    
  contains


  subroutine read_input()

    !reads the no of valence and conduction states from "bsemat.h5"

    
    character(len=9), parameter :: filename = "bsemat.h5"
    character(len=21), parameter :: dsetname_1 = "/bse_header/bands/nvb"
    character(len=21), parameter :: dsetname_2 = "/bse_header/bands/ncb"
    !character(len=25), parameter :: dsetname_3 = "/mf_header/crystal/celvol"
  
    integer(hid_t) :: file_id
    integer(hid_t) :: dset1_id,dset2_id,dset3_id,dset4_id
    integer     ::   error
    integer(hsize_t), DIMENSION(0) :: data1_dims,data2_dims,data3_dims,data4_dims
  
  
  
    call h5open_f(error)
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
  
    call h5dopen_f(file_id, dsetname_1, dset1_id, error)
    call h5dread_f(dset1_id, H5T_NATIVE_INTEGER, sys%nv, data1_dims, error)
    call h5dclose_f(dset1_id, error)
  
    call h5dopen_f(file_id, dsetname_2, dset2_id, error)
    call h5dread_f(dset2_id, H5T_NATIVE_INTEGER, sys%nc, data2_dims, error)
    call h5dclose_f(dset2_id, error)
 
    call h5fclose_f(file_id, error)
    call h5close_f(error)
  
    !debug
    print *,"conduction states,valence states",sys%nc,sys%nv
  end subroutine read_input
  subroutine read_triplet_exciton_eigenvalues()
        !reads eigenvalues from eigenvalues.dat file
        
        integer :: error
        character(len=256) :: line
        character(len=20) :: filename
        integer :: i, n, ios
         integer :: iQ, x, iunit
        double precision :: val1, val2, val3, val4  
        double precision, allocatable :: exciton_eigenvalues(:),exciton_eigenvalues_Q(:,:)

        !allocate(exciton_eigenvalues(maxlines))
        !allocate(exciton_eigenvalues_Q(sys%neb,sys%nQ))
        allocate(exciton_sys%eigenvalues_t(sys%neb,sys%nQ))

        do iQ = 1, sys%nQ
          write(filename, '("eigenvalues_", I0, ".dat")') iQ
          iunit = 10 + iQ  ! Avoid reusing unit numbers
          open(unit=iunit, file=filename, status='old', action='read', iostat=error)
          if(error /= 0) then
              print *, "Error opening file: ", filename
              exit
          end if
              READ(iunit,'(A)', IOSTAT=ios) line
              READ(iunit,'(A)', IOSTAT=ios) line
              READ(iunit,'(A)', IOSTAT=ios) line
               READ(iunit,'(A)', IOSTAT=ios) line
               if (ios /= 0) then
                     print *, "Error reading header from file hey hey: ", filename
                     exit
               end if
              do x = 1,sys%neb
                  READ(iunit, *, IOSTAT=ios) val1, val2, val3, val4
                  if (ios /= 0) then
                      print *, "Error reading eigenvalue from file: ", filename
                      exit
                  end if
                  exciton_sys%eigenvalues_t(x,iQ) = val1 / 27.2114079527 ! Convert to eV
                  print *, "Eigenvalue for Q index ", iQ, " band ", x, ": ", exciton_sys%eigenvalues_t(x,iQ)
              end do
          close(iunit)
         end do
   end subroutine read_triplet_exciton_eigenvalues

    subroutine read_triplet_exciton_eigenvectors()
        !reads eigenvectors from eigenvec.h5 file      
        character(len=20) :: filename 
        character(len=25), parameter :: dsetname= "exciton_data/eigenvectors"
        integer(HID_T) :: file_id, dset_id, space_id
        integer :: error,rank
        integer(HSIZE_T), allocatable :: dims(:), maxdims(:)
        double precision,allocatable :: tmp(:,:,:,:,:,:,:)
        integer :: i
        integer :: i1, i2, i3, i4, i5, i6, i7
        allocate(exciton_sys%eigenvectors_t(1,sys%nv,sys%nc,sys%nk,sys%neb,sys%nQ))
        exciton_sys%eigenvectors_t = cmplx(0.0, 0.0) ! Initialize to zero
        call h5open_f(error)
        allocate(dims(7))
        dims(1) = 2 ! real and imaginary parts
        dims(2) = 1
        dims(3) = sys%nv ! valence band
        dims(4) = sys%nc ! conduction band
        dims(5) = sys%nk ! kpoints
        dims(6) = sys%neb ! number of exciton bands
        dims(7) = 1 ! number of Q points
        allocate(tmp(2,1,sys%nv,sys%nc,sys%nk,sys%neb,1))
        do i = 1, sys%nQ
        !open hdf5 file name eigenvectors_$i.h5
          write(filename, '("eigenvectors_", I0, ".h5")') i
         
          call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)   
          if(error /= 0) then
            print *, "Error opening file: ", filename
            exit
          end if
         
          !call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
          call h5dopen_f(file_id, dsetname, dset_id, error)
            if(error /= 0) then
                print *, "Error opening dataset: ", dsetname
                exit
            end if
         ! call h5dget_space_f(dset_id, space_id, error)
         !   if(error /= 0) then
         !       print *, "Error getting space for dataset: ", dsetname
         !       exit
         !   end if
         ! call h5sget_simple_extent_ndims_f(space_id, rank, error)
         !   if(error /= 0) then
         !       print *, "Error getting rank of dataset: ", dsetname
         !       exit
         !   end if
         ! print *, "Rank of the dataset: ", rank
         ! allocate(dims(rank))
          !allocate(maxdims(rank))
          !call h5sget_simple_extent_dims_f(space_id, dims,maxdims, error)
        !print *, "Dimensions of the dataset: ", dims(7), dims(6), dims(5), dims(4), dims(3), dims(2), dims(1)
          
           
           !allocate(exciton_sys%eigenvectors_t(dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))
           print *, "Dimensions of the eigenvectors: ", dims(1),dims(2), dims(3), dims(4), dims(5), dims(6), dims(7)
        
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmp, dims, error)
            if(error /= 0) then
                print *, "Error reading eigenvectors from dataset: ", dsetname
                exit
            end if
         !call h5sclose_f(space_id, error)
         call h5dclose_f(dset_id, error)
         call h5fclose_f(file_id, error)
        
         
        print *, "tmp :", tmp(1,1,2,1,1,1,1), tmp(2,1,2,1,1,1,1)

         do i2= 1, dims(2) ! spin
            do i3= 1, dims(3) ! valence band
                do i4= 1, dims(4) ! conduction band
                    do i5= 1, dims(5) ! kpoints
                        do i6= 1, dims(6) !neb
                           !do i7 = 1, dims(7)
                                 exciton_sys%eigenvectors_t(i2,i3,i4,i5,i6,i) = cmplx(tmp(1,i2,i3,i4,i5,i6,1), &
                                      tmp(2,i2,i3,i4,i5,i6,1), kind=16)
                                                     
                                !exciton_eigenvector(i2,i3,i4,i5,i6) = cmplx(tmp(1,i2,i3,i4,i5,i6,1), &
                                 !   tmp(2,i2,i3,i4,i5,i6,1), kind=8)
                                print *, "Eigenvector for Q index ", i, " spin ", i2, " valence ", i3, " conduction ", i4, " kpoints ", i5,"band",i6, ": ", exciton_sys%eigenvectors_t(i2,i3,i4,i5,i6,i)
                           !end do
                       end do
                   end do
               end do
            end do
        end do
        !deallocate(tmp)
        !deallocate(dims)
        !deallocate(maxdims)
       
        
      end do
      call h5close_f(error)  
   

    end subroutine read_triplet_exciton_eigenvectors



 ! subroutine read_singlet_exciton_eigenvalues()
 !   !reads eigenvalues from eig.dat file
 !   !and reads the no of excitons
 !   integer :: i,temp_ex
 !   integer :: error
 !   double precision :: x
 !   OPEN(UNIT=10, FILE="eig_s.dat", STATUS="OLD", ACTION="READ", IOSTAT=error)
 !   IF(error/=0) STOP "Error opening file eig.dat"
 !   READ(10, *, IOSTAT=error) temp_ex
 !   IF(error/=0) STOP "Error reading number of excitons"
 !   
 !    
 !   !debug
 !    PRINT *, "Number of excitons: ", temp_ex
!
 !   allocate(exciton_sys%eigenvalues_s(sys%nex))
 !   
 !   do i = 1, sys%nex
 !     !READ(10, *, IOSTAT=error)
 !     !if (error /= 0) stop "Error reading eigenvalue"
 !     READ(10, *, IOSTAT=error),x
 !     exciton_sys%eigenvalues_s(i) = (x/27.2114079527) 
 !     if (error /= 0) stop "Error reading eigenvalue"
 !   end do
 !   CLOSE(10)
!
!end subroutine read_singlet_exciton_eigenvalues
!
 ! subroutine read_singlet_exciton_eigenvectors()
 !   !reads eigenvectors from eigenvec.dat file
 !   integer :: i,j,k
 !   integer :: error
!
 !   open(UNIT=10, FILE="eigvec_s.dat", STATUS="OLD", ACTION="READ", IOSTAT=error)
 !   if(error/=0) stop "Error opening file eigvec.dat"
!
 !   allocate(exciton_sys%eigenvectors_s(sys%nex,sys%nv,sys%nc))
 !  print*,"hello3"
 !   do i = 1,sys%nex
 !     do j = 1,(sys%nv)
 !        do k = 1,(sys%nc)
 !           read(10, *, IOSTAT=error) exciton_sys%eigenvectors_s(i,j,k)
 !           !print*,i,j,k,exciton_sys%eigenvectors(i,j,k)
 !           if (error /= 0) stop "Error reading eigenvector"
 !        end do
 !     end do
 !   end do
!
 !   !print*,exciton_sys%eigenvectors(5,20)
 !   close(10)
!
!
!
 ! end subroutine read_singlet_exciton_eigenvectors
 ! subroutine read_triplet_exciton_eigenvalues()
 !  !reads eigenvalues from eig.dat file
 !  !and reads the no of excitons
 !  integer :: i,temp_ex
 !  integer :: error
 !  double precision :: x
 !  OPEN(UNIT=10, FILE="eig_t.dat", STATUS="OLD", ACTION="READ", IOSTAT=error)
 !  IF(error/=0) STOP "Error opening file eig.dat"
 !  READ(10, *, IOSTAT=error) temp_ex
 !  IF(error/=0) STOP "Error reading number of excitons"
 !  
 !   
 !  !debug
 !   PRINT *, "Number of excitons: ", temp_ex
!
 !  allocate(exciton_sys%eigenvalues_t(sys%nex))
 !  
 !  do i = 1, sys%nex
 !    !READ(10, *, IOSTAT=error)
 !    !if (error /= 0) stop "Error reading eigenvalue"
 !    READ(10, *, IOSTAT=error),x
 !    exciton_sys%eigenvalues_t(i) = (x/27.2114079527) 
 !    if (error /= 0) stop "Error reading eigenvalue"
 !  end do
 !  CLOSE(10)
!
!end subroutine read_triplet_exciton_eigenvalues
!
 !subroutine read_triplet_exciton_eigenvectors()
 !  !reads eigenvectors from eigenvec.dat file
 !  integer :: i,j,k
 !  integer :: error
!
 !  open(UNIT=10, FILE="eigvec_t.dat", STATUS="OLD", ACTION="READ", IOSTAT=error)
 !  if(error/=0) stop "Error opening file eigvec.dat"
!
 !  allocate(exciton_sys%eigenvectors_t(sys%nex,sys%nv,sys%nc))
 ! print*,"hello3"
 !  do i = 1,sys%nex
 !    do j = 1,(sys%nv)
 !       do k = 1,(sys%nc)
 !          read(10, *, IOSTAT=error) exciton_sys%eigenvectors_t(i,j,k)
 !          !print*,i,j,k,exciton_sys%eigenvectors(i,j,k)
 !          if (error /= 0) stop "Error reading eigenvector"
 !       end do
 !    end do
 !  end do
!
 !  !print*,exciton_sys%eigenvectors(5,20)
 !  close(10)
!
!
!
 !end subroutine read_triplet_exciton_eigenvectors
 !
  subroutine read_bsemat_d(bsemat)
    complex(kind=8), dimension(sys%nv,sys%nv,sys%nc,sys%nc,sys%nk,sys%nk,sys%nQ),intent(inout) :: bsemat
    character(len=20) :: filename 
    character(len=10), parameter :: name_1 = "/mats/head"
    character(len=10), parameter :: name_2 = "/mats/wing"
    character(len=10), parameter :: name_3 = "/mats/body"
    character(len=14), parameter :: name_4 = "/mats/exchange"
    double precision, DIMENSION(:,:,:,:,:,:,:), allocatable :: data1_out
    integer :: i,iQ,ik,ikp
    integer(hid_t) :: file_id
    integer(hid_t) :: dset1_id,dset2_id
    integer     ::   error,iv,ivp,ic,icp
    integer :: imatrix
    integer(hsize_t), DIMENSION(:), allocatable :: data1_dims
  
    allocate(data1_dims(7))
    allocate(data1_out(2,sys%nv,sys%nv,sys%nc,sys%nc,sys%nk,sys%nk))
    data1_dims(1) = 2
    data1_dims(2) = sys%nv
    data1_dims(3) = sys%nv
    data1_dims(4) = sys%nc
    data1_dims(5) = sys%nc
    data1_dims(6) = sys%nk
    data1_dims(7) = sys%nk
    
    call h5open_f(error)
    do iQ = 1,sys%nQ
      write(filename, '("bsemat_", I0, ".h5")') iQ
      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
      do imatrix = 3,3
        if (imatrix== 1) then
         call h5dopen_f(file_id, name_1, dset1_id, error)
        elseif (imatrix== 2) then
          call h5dopen_f(file_id, name_2, dset1_id, error)
        elseif (imatrix== 3) then
         call h5dopen_f(file_id, name_3, dset1_id, error)
        elseif (imatrix== 4) then
         call h5dopen_f(file_id, name_4, dset1_id, error)
        end if
      
        call h5dread_f(dset1_id, H5T_NATIVE_DOUBLE, data1_out, data1_dims, error)
        call h5dclose_f(dset1_id, error)
        do icp=1,sys%nc
         do ic=1,sys%nc
           do ivp=1,sys%nv
              do iv=1,sys%nv
                 do ik = 1,sys%nk
                    do ikp = 1,sys%nk
                      bsemat(iv,ivp,ic,icp,ik,ikp,iQ) = bsemat(iv,ivp,ic,icp,ik,ikp,iQ) +cmplx(data1_out(1,iv,ivp,ic,icp,ik,ikp), &
                                                              data1_out(2,iv,ivp,ic,icp,ik,ikp), kind=8)
                      !print*,bsemat(iv,ivp,ic,icp,ik,ikp,iQ),iv,ivp,ic,icp,ik,ikp,iQ
                    end do
                  end do
            
              enddo
          enddo
       enddo
      enddo
    !print*,bse_mat(1,1,2,1),imatrix
    end do !imatrix loop
    call h5fclose_f(file_id, error)
    !)
    end do !iQ loop

    call h5close_f(error)

    deallocate(data1_out)
  
  end subroutine read_bsemat_d

  subroutine read_bsemat_ee(bsemat_ee)
    complex(kind=8), dimension(sys%nc,sys%nc,sys%nc,sys%nc,sys%nk,sys%nk,sys%nQ),intent(inout) :: bsemat_ee
    character(len=20) :: filename 
    character(len=13), parameter :: name_1 = "/mats/head_ee"
    character(len=13), parameter :: name_2 = "/mats/wing_ee"
    character(len=13), parameter :: name_3 = "/mats/body_ee"
    character(len=14), parameter :: name_4 = "/mats/exchange"
    double precision, DIMENSION(:,:,:,:,:,:,:), allocatable :: data1_out
    integer :: i,iQ,ik,ikp
    integer :: imatrix
    integer(hid_t) :: file_id
    integer(hid_t) :: dset1_id,dset2_id
    integer     ::   error,iv,ivp,ic,icp
    integer(hsize_t), DIMENSION(:), allocatable :: data1_dims
  
    allocate(data1_dims(7))
    data1_dims(1) = 2
    data1_dims(2) = sys%nc
    data1_dims(3) = sys%nc
    data1_dims(4) = sys%nc
    data1_dims(5) = sys%nc
    data1_dims(6) = sys%nk
    data1_dims(7) = sys%nk
     allocate(data1_out(2,sys%nc,sys%nc,sys%nc,sys%nc,sys%nk,sys%nk))
    call h5open_f(error)
    do iQ = 1,sys%nQ
      
        !open hdf5 file name eigenvectors_$i.h5
    write(filename, '("bsemat_", I0, ".h5")') iQ
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
    do imatrix = 1,3
      if (imatrix== 1) then
       call h5dopen_f(file_id, name_1, dset1_id, error)
      elseif (imatrix== 2) then
       call h5dopen_f(file_id, name_2, dset1_id, error)
      elseif (imatrix== 3) then
       call h5dopen_f(file_id, name_3, dset1_id, error)
      elseif (imatrix== 4) then
       call h5dopen_f(file_id, name_4, dset1_id, error)
      end if
     
  
      call h5dread_f(dset1_id, H5T_NATIVE_DOUBLE, data1_out, data1_dims, error)
      call h5dclose_f(dset1_id, error)
      do icp=1,sys%nc
       do ic=1,sys%nc
          do ivp=1,sys%nc
              do iv=1,sys%nc
                 do ik = 1,sys%nk
                    do ikp = 1,sys%nk
                      bsemat_ee(iv,ivp,ic,icp,ik,ikp,iQ) = bsemat_ee(iv,ivp,ic,icp,ik,ikp,iQ) +cmplx(data1_out(1,iv,ivp,ic,icp,ik,ikp), &
                                                              data1_out(2,iv,ivp,ic,icp,ik,ikp), kind=8)
                     ! print*,bsemat_ee(iv,ivp,ic,icp,ik,ikp,iQ),iv,ivp,ic,icp,ik,ikp,iQ
                    end do
                  end do
            
              enddo
          enddo
       enddo
    enddo
    !print*,bse_mat(1,1,2,1),imatrix
    end do !imatrix loop
    call h5fclose_f(file_id, error)
   
    end do !iQ loop
    call h5close_f(error)
    deallocate(data1_out)
  
  end subroutine read_bsemat_ee

    subroutine read_bsemat_hh(bsemat_hh)
    complex(kind=8), dimension(sys%nv,sys%nv,sys%nv,sys%nv,sys%nk,sys%nk,sys%nQ),intent(inout) :: bsemat_hh
    character(len=20) :: filename 
    character(len=13), parameter :: name_1 = "/mats/head_hh"
    character(len=13), parameter :: name_2 = "/mats/wing_hh"
    character(len=13), parameter :: name_3 = "/mats/body_hh"
    character(len=13), parameter :: name_4 = "/mats/exchange"
    double precision, DIMENSION(:,:,:,:,:,:,:), allocatable :: data1_out
    integer :: i,iQ,ik,ikp
    integer :: imatrix
    integer(hid_t) :: file_id
    integer(hid_t) :: dset1_id,dset2_id
    integer     ::   error,iv,ivp,ic,icp
    integer(hsize_t), DIMENSION(:), allocatable :: data1_dims
    allocate(data1_out(2,sys%nv,sys%nv,sys%nv,sys%nv,sys%nk,sys%nk))
    allocate(data1_dims(7))
    data1_dims(1) = 2
    data1_dims(2) = sys%nv
    data1_dims(3) = sys%nv
    data1_dims(4) = sys%nv
    data1_dims(5) = sys%nv
    data1_dims(6) = sys%nk
    data1_dims(7) = sys%nk
    
    call h5open_f(error)
    do iQ = 1,sys%nQ
      
        !open hdf5 file name eigenvectors_$i.h5
    write(filename, '("bsemat_", I0, ".h5")') iQ
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
    do imatrix = 1,3
      if (imatrix== 1) then
       call h5dopen_f(file_id, name_1, dset1_id, error)
      elseif (imatrix== 2) then
       call h5dopen_f(file_id, name_2, dset1_id, error)
      elseif (imatrix== 3) then
       call h5dopen_f(file_id, name_3, dset1_id, error)
      elseif (imatrix== 4) then
       call h5dopen_f(file_id, name_4, dset1_id, error)
      end if
      
  
      call h5dread_f(dset1_id, H5T_NATIVE_DOUBLE, data1_out, data1_dims, error)
      call h5dclose_f(dset1_id, error)
      do icp=1,sys%nv
       do ic=1,sys%nv
          do ivp=1,sys%nv
              do iv=1,sys%nv
                 do ik = 1,sys%nk
                    do ikp = 1,sys%nk
                      bsemat_hh(iv,ivp,ic,icp,ik,ikp,iQ) = bsemat_hh(iv,ivp,ic,icp,ik,ikp,iQ) +cmplx(data1_out(1,iv,ivp,ic,icp,ik,ikp), &
                                                              data1_out(2,iv,ivp,ic,icp,ik,ikp), kind=8)
                      !print*,bsemat(iv,ivp,ic,icp),iv,ivp,ic,icp
                    end do
                  end do
            
              enddo
          enddo
       enddo
    enddo
    !print*,bse_mat(1,1,2,1),imatrix
    end do !imatrix loop
    call h5fclose_f(file_id, error)
    !call h5close_f(error)
    end do !iQ loop
    call h5close_f(error)
    deallocate(data1_out)
  
  end subroutine read_bsemat_hh

  
  

!
 ! subroutine read_bsemat_x(bsemat)
 !  double precision, dimension(sys%nb,sys%nb,sys%nb,sys%nb),intent(inout) :: bsemat
 !  character(len=9), parameter :: filename = "bsemat.h5"
 !  character(len=10), parameter :: name_1 = "/mats/head"
 !  character(len=10), parameter :: name_2 = "/mats/wing"
 !  character(len=10), parameter :: name_3 = "/mats/body"
 !  character(len=14), parameter :: name_4 = "/mats/exchange"
 !  double precision, DIMENSION(:,:,:,:,:,:,:), allocatable :: data1_out
 !
 !  integer(hid_t) :: file_id
 !  integer(hid_t) :: dset1_id,dset2_id
 !  integer     ::   error,iv,ivp,ic,icp,imatrix
 !  integer(hsize_t), DIMENSION(:), allocatable :: data1_dims
 !
 !  allocate(data1_dims(7))
 !  data1_dims(1) = 2
 !  data1_dims(2) = sys%nb
 !  data1_dims(3) = sys%nb
 !  data1_dims(4) = sys%nb
 !  data1_dims(5) = sys%nb
 !  data1_dims(6) = 1
 !  data1_dims(7) = 1
 !  imatrix = 4
 !  call h5open_f(error)
 !  call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
 !  if (imatrix== 1) then
 !     call h5dopen_f(file_id, name_1, dset1_id, error)
 !  elseif (imatrix== 2) then
 !     call h5dopen_f(file_id, name_2, dset1_id, error)
 !  elseif (imatrix== 3) then
 !     call h5dopen_f(file_id, name_3, dset1_id, error)
 !  elseif (imatrix== 4) then
 !     call h5dopen_f(file_id, name_4, dset1_id, error)
 !  end if
 !  allocate(data1_out(2,sys%nb,sys%nb,sys%nb,sys%nb,1,1))
 !  print*,sys%nb
 !  call h5dread_f(dset1_id, H5T_NATIVE_DOUBLE, data1_out, data1_dims, error)
 !  call h5dclose_f(dset1_id, error)
 !  do icp=1,sys%nb
 !     do ic=1,sys%nb
 !        do ivp=1,sys%nb
 !            do iv=1,sys%nb
 !                bsemat(iv,ivp,ic,icp) = data1_out(1,iv,ivp,ic,icp,1,1)
 !                !print*,bsemat(iv,ivp,ic,icp),iv,ivp,ic,icp
 !            enddo
 !        enddo
 !     enddo
 !  enddo
 !  !print*,bse_mat(1,1,2,1),imatrix
 !  call h5fclose_f(file_id, error)
 !  call h5close_f(error)
 !  deallocate(data1_out)
 !
 !end subroutine read_bsemat_x

  
    
end module input_read



  
    
