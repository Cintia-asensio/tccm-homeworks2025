program dynamics
    implicit none

    integer :: input_file
    integer :: Natoms
    integer :: ios
    integer :: i, j

    double precision, allocatable :: coord(:, :)
    double precision, allocatable :: mass(:)
    double precision, allocatable :: distance(:, :)

    ! Assign file unit
    input_file = 10

    ! Open the file
    open(unit=input_file, file="inp.txt", status="old", action="read", iostat=ios)

    if (ios /= 0) then
        print *, "Error: The file inp.txt cannot be opened!"
        stop
    end if

    ! Read number of atoms
    Natoms = read_Natoms(input_file)
    print *, "Number of atoms:", Natoms
    print *

    ! Allocate arrays
    allocate(coord(Natoms, 3))
    allocate(mass(Natoms))
    allocate(distance(Natoms, Natoms))

    ! Read coordinates and masses
    call read_molecule(input_file, Natoms, coord, mass)

    close(input_file)

    ! ---- PRINT COORD AND MASS ----
    print *, "Atomic coordinates and masses:"
    print *, "Atom        x           y           z         mass"

    do i = 1, Natoms
        write(*,'(i4,4f12.6)') i, coord(i,1), coord(i,2), coord(i,3), mass(i)
    end do
    print *

    ! ---- COMPUTE DISTANCES ----
    call compute_distances(Natoms, coord, distance)

    ! ---- PRINT DISTANCE MATRIX ----
    print *, "Internuclear distance matrix:"
    do i = 1, Natoms
        write(*,'(100f10.4)') (distance(i,j), j=1,Natoms)
    end do

contains

    integer function read_Natoms(input_file) result(Natoms)
        implicit none
        integer, intent(in) :: input_file

        read(input_file, *) Natoms
    end function read_Natoms


    subroutine read_molecule(input_file, Natoms, coord, mass)
        implicit none
        integer, intent(in) :: input_file
        integer, intent(in) :: Natoms
        double precision, intent(out) :: coord(Natoms,3)
        double precision, intent(out) :: mass(Natoms)

        integer :: i

        do i = 1, Natoms
            read(input_file, *) coord(i,1), coord(i,2), coord(i,3), mass(i)
        end do

    end subroutine read_molecule


    subroutine compute_distances(Natoms, coord, distance)
        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: coord(Natoms,3)
        double precision, intent(out) :: distance(Natoms,Natoms)

        integer :: i, j
        double precision :: dx, dy, dz

        do i = 1, Natoms
            do j = 1, Natoms
                dx = coord(i,1) - coord(j,1)
                dy = coord(i,2) - coord(j,2)
                dz = coord(i,3) - coord(j,3)

                distance(i,j) = sqrt(dx*dx + dy*dy + dz*dz)
            end do
        end do

    end subroutine compute_distances

end program dynamics

