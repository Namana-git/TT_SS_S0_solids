module kpoints
    use global_variables
    use hdf5
    implicit none

    contains

    subroutine read_Qpoints()
    ! read kpoint sfrom eigenvectors.h5 file
  
    integer(hid_t) :: file_id
    integer(hid_t) :: dset1_id, dset2_id
    integer :: i
    integer :: error
    integer(hsize_t), DIMENSION(0) :: data1_dims
    integer(hsize_t), DIMENSION(2) :: data2_dims
    double precision, allocatable :: Qpts(:,:)
    character(len=40) :: dsetname2="/exciton_header/kpoints/exciton_Q_shifts"
    character(len=20) :: filename
    allocate(sys%Qpts(3, sys%nQ))
    allocate(Qpts(3, 1))
    call h5open_f(error)
    do i = 1, sys%nQ
        !open hdf5 file name eigenvectors_$i.h5
        write(filename, '("eigenvectors_", I0, ".h5")') i
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)   
        if(error /= 0) then
            print *, "Error opening file: ", filename
            exit
        end if
        !allocate Qpts
    
        !Qpts = 0.0d0
        data2_dims(1) = 3
        data2_dims(2) = 1
        !open,read,the dateset that contains tohe Q points
        call h5dopen_f(file_id, dsetname2, dset2_id, error) 
        call h5dread_f(dset2_id, H5T_NATIVE_DOUBLE, Qpts, data2_dims, error)
        call h5dclose_f(dset2_id, error)
        sys%Qpts(:, i) = Qpts(:, 1)
        print*, "Q-point ", i, ": ", sys%Qpts(1, i), sys%Qpts(2, i), sys%Qpts(3, i)
        !close the file
        call h5fclose_f(file_id, error)
    end do

    end subroutine read_Qpoints

        subroutine read_kpoints()
    ! read kpoint sfrom eigenvectors.h5 file
    character(len=15), parameter :: filename = "eigenvectors.h5"
    character(len=26), parameter :: dsetname1 = "/exciton_header/kpoints/nk"
    character(len=28), parameter :: dsetname2 = "/exciton_header/kpoints/kpts"
    integer(hid_t) :: file_id
    integer(hid_t) :: dset1_id, dset2_id
    integer :: i
    integer :: error
    integer(hsize_t), DIMENSION(0) :: data1_dims
    integer(hsize_t), DIMENSION(2) :: data2_dims

    ! Open the HDF5 file
    call h5open_f(error)
    call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

    !read the number of kpoints
    call h5dopen_f(file_id, dsetname1, dset1_id, error)
    call h5dread_f(dset1_id, H5T_NATIVE_INTEGER, sys%nk, data1_dims, error)
    call h5dclose_f(dset1_id, error)

    !read the kpoints

    data2_dims(1) = 3
    data2_dims(2) = sys%nk
    allocate(sys%kpts(data2_dims(1), data2_dims(2)))
    call h5dopen_f(file_id, dsetname2, dset2_id, error) 
    call h5dread_f(dset2_id, H5T_NATIVE_DOUBLE, sys%kpts, data2_dims, error)
    call h5dclose_f(dset2_id, error)


    ! Close the HDF5 file   

    do i = 1, sys%nk
        ! fold the k point into first brillouin zone
        print*, "K-point ", i, ": ", sys%kpts(1, i), sys%kpts(2, i), sys%kpts(3, i)
    end do

    end subroutine read_kpoints
   
    ! soubroutine to convert a k point value to index. also if the 
 
    subroutine kpoint_to_index(kpt, index)
        double precision, intent(inout) :: kpt(3)
        integer, intent(out) :: index
        double precision :: dist
        integer :: i

        ! Initialize the index to -1 (not found)
        index = -1
        ! fold the point into firtst brillion zone from 0 to 1
        !print*,anint(kpt)
        kpt = kpt - floor(kpt+0.5)
        !kpt = kpt - anint(kpt)
        !kpt = mod(kpt + 0.5, 1.0) - 0.5
        !print*, "kpt", kpt(1),kpt(2),kpt(3)
        ! Loop through all kpoints to find the closest match
        do i = 1, sys%nk
            dist = sqrt(sum((sys%kpts(:, i) - kpt)**2))
            if (dist < 1.0e-6) then  ! tolerance for matching
                index = i
                return
            end if
        end do
        if(index==-1)then
            print*,"kpoint not found"
        end if

    end subroutine kpoint_to_index

    subroutine Qpoint_to_index(Qpt, index)
        double precision, intent(inout) :: Qpt(3)
        integer, intent(out) :: index
        double precision :: dist
        integer :: i

        ! Initialize the index to -1 (not found)
        index = -1
        ! fold the point into firtst brillion zone from 0 to 1
        Qpt = Qpt - floor(Qpt+0.5)
        ! Loop through all kpoints to find the closest match
        do i = 1, sys%nQ
            dist = sqrt(sum((sys%Qpts(:, i) - Qpt)**2))
            if (dist < 1.0e-6) then  ! tolerance for matching
                index = i
                return
            end if
        end do

    end subroutine Qpoint_to_index







end module kpoints