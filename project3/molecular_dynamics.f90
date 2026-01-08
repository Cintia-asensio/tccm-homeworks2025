program dynamics
    implicit none

    integer :: input_file
    integer :: Natoms
    integer :: ios
    integer :: i, j

    double precision, parameter :: angstrom_to_nm = 0.1d0

    double precision :: epsilon, sigma
    double precision :: Vtot, Ttot, Etot

    double precision, allocatable :: coord(:, :)
    double precision, allocatable :: mass(:)
    double precision, allocatable :: distance(:, :)
    double precision, allocatable :: velocity(:, :)

    ! -----------------------
    ! Lennard-Jones parameters
    ! -----------------------
    epsilon = 0.997d0
    sigma   = 3.405d0 * angstrom_to_nm

    ! -----------------------
    ! Open input file
    ! -----------------------
    input_file = 10
    open(unit=input_file, file="inp.txt", status="old", action="read", iostat=ios)

    if (ios /= 0) then
        print *, "Error: cannot open inp.txt"
        stop
    end if

    ! -----------------------
    ! Read number of atoms
    ! -----------------------
    Natoms = read_Natoms(input_file)
    print *, "Number of atoms:", Natoms
    print *

    ! -----------------------
    ! Allocate arrays
    ! -----------------------
    allocate(coord(Natoms,3))
    allocate(mass(Natoms))
    allocate(distance(Natoms,Natoms))
    allocate(velocity(Natoms,3))

    ! -----------------------
    ! Read molecule
    ! -----------------------
    call read_molecule(input_file, Natoms, coord, mass)
    close(input_file)

    ! -----------------------
    ! Initialize velocities to zero
    ! -----------------------
    velocity = 0.0d0

    ! -----------------------
    ! Print coordinates and masses
    ! -----------------------
    print *, "Atomic coordinates (nm) and masses (g/mol):"
    print *, "Atom        x(nm)       y(nm)       z(nm)      mass"
    do i = 1, Natoms
        write(*,'(i4,4f12.6)') i, coord(i,1), coord(i,2), coord(i,3), mass(i)
    end do
    print *

    ! -----------------------
    ! Compute distances
    ! -----------------------
    call compute_distances(Natoms, coord, distance)

    ! -----------------------
    ! Compute energies
    ! -----------------------
    Vtot = V(epsilon, sigma, Natoms, distance)
    Ttot = T(Natoms, velocity, mass)
    Etot = Ttot + Vtot

    ! -----------------------
    ! Print energies
    ! -----------------------
    print *, "Total kinetic energy T (kJ/mol):", Ttot
    print *, "Total potential energy V (kJ/mol):", Vtot
    print *, "Total energy E = T + V (kJ/mol):", Etot

contains

    ! =====================================================
    integer function read_Natoms(input_file) result(Natoms)
        implicit none
        integer, intent(in) :: input_file
        read(input_file, *) Natoms
    end function read_Natoms
    ! =====================================================


    ! =====================================================
    subroutine read_molecule(input_file, Natoms, coord, mass)
        implicit none
        integer, intent(in) :: input_file
        integer, intent(in) :: Natoms
        double precision, intent(out) :: coord(Natoms,3)
        double precision, intent(out) :: mass(Natoms)

        integer :: i
        double precision, parameter :: angstrom_to_nm = 0.1d0

        do i = 1, Natoms
            read(input_file, *) coord(i,1), coord(i,2), coord(i,3), mass(i)
            coord(i,:) = coord(i,:) * angstrom_to_nm
        end do
    end subroutine read_molecule
    ! =====================================================


    ! =====================================================
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
    ! =====================================================


    ! =====================================================
    double precision function V(epsilon, sigma, Natoms, distance)
        implicit none
        double precision, intent(in) :: epsilon, sigma
        integer, intent(in) :: Natoms
        double precision, intent(in) :: distance(Natoms,Natoms)

        integer :: i, j
        double precision :: r, sr6, sr12

        V = 0.0d0
        do i = 1, Natoms
            do j = i+1, Natoms
                r = distance(i,j)
                sr6  = (sigma / r)**6
                sr12 = sr6 * sr6
                V = V + 4.0d0 * epsilon * (sr12 - sr6)
            end do
        end do
    end function V
    ! =====================================================


    ! =====================================================
    double precision function T(Natoms, velocity, mass)
        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: velocity(Natoms,3)
        double precision, intent(in) :: mass(Natoms)

        integer :: i
        double precision :: v2

        T = 0.0d0
        do i = 1, Natoms
            v2 = velocity(i,1)**2 + velocity(i,2)**2 + velocity(i,3)**2
            T = T + 0.5d0 * mass(i) * v2
        end do
    end function T
    ! =====================================================

end program dynamics

