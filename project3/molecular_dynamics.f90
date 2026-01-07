program dynamics
    implicit none

    integer :: input_file
    integer :: Natoms
    integer :: ios

    double precision, allocatable :: coord(:, :)
    double precision, allocatable :: mass(:)

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

    ! Allocate arrays
    allocate(coord(Natoms, 3))
    allocate(mass(Natoms))

    ! Read coordinates and masses
    call read_molecule(input_file, Natoms, coord, mass)

    close(input_file)

    ! (Optional) print to check
    print *, "Atom   x       y       z       mass"
    do ios = 1, Natoms
        print *, ios, coord(ios,1), coord(ios,2), coord(ios,3), mass(ios)
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

        ! Read coordinates and masses
        do i = 1, Natoms
            read(input_file, *) coord(i,1), coord(i,2), coord(i,3), mass(i)
        end do

    end subroutine read_molecule

end program dynamics

