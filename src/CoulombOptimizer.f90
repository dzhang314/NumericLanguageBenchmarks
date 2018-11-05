program CoulombOptimizer

    use iso_fortran_env, only : real64
    implicit none

    integer, parameter :: dp = real64
    integer :: arglen, num_points, dim, num_steps, step, i, u
    character(len=:), allocatable :: filename, num_steps_str
    real(dp) :: energy, step_size
    real(dp), allocatable :: points(:,:), forces(:,:), workspace(:,:)

    if (command_argument_count() < 2) then
        error stop "Needs 2 arguments: <input-file-name> <num-steps>"
    end if

    call get_command_argument(1, length=arglen)
    allocate (character(len=arglen) :: filename)
    call get_command_argument(1, value=filename)

    open(newunit=u, file=filename, action='READ')
    read(u,*) num_points, dim
    allocate(points(dim, num_points), forces(dim, num_points), &
           & workspace(dim, num_points))
    read(u,*) points
    close(u)

    call get_command_argument(2, length=arglen)
    allocate (character(len=arglen) :: num_steps_str)
    call get_command_argument(2, value=num_steps_str)

    read(num_steps_str,*) num_steps
    step_size = 1.0e-3_dp;
    call coulomb_forces(points, energy, forces)
    call normalize_forces(points, forces)

    do step = 1, num_steps
        call quadratic_line_search(points, energy, forces, step_size, workspace)
        points = points + step_size * forces
        call normalize_points(points)
        call coulomb_forces(points, energy, forces)
        call normalize_forces(points, forces)
    end do

    write(*,*) energy
    doi = 1, num_points
        write(*,*) points(:,i)
    end do

contains

    pure subroutine normalize_points(points)
        real(dp), intent(inout) :: points(:,:)
        integer :: i
        do i = 1, size(points, 2)
            points(:,i) = points(:,i) / norm2(points(:,i))
        end do
    end subroutine normalize_points

    pure subroutine normalize_forces(points, forces)
        real(dp), intent(in) :: points(:,:)
        real(dp), intent(inout) :: forces(:,:)
        integer :: i
        do i = 1, size(forces, 2)
            forces(:,i) = forces(:,i) - &
                & dot_product(forces(:,i), points(:,i)) * points(:,i)
        end do
    end subroutine normalize_forces

    pure real(dp) function coulomb_energy(points) result (energy)
        real(dp), intent(in) :: points(:,:)
        integer :: i, j
        energy = 0.0
        do i = 2, size(points, 2)
            do j = 1, i - 1
                energy = energy + 1.0_dp / norm2(points(:,i) - points(:,j))
            end do
        end do
    end function coulomb_energy

    pure subroutine coulomb_forces(points, energy, forces)
        real(dp), intent(in) :: points(:,:)
        real(dp), intent(out) :: energy, forces(:,:)
        real(dp) :: displ(size(points, 1)), inv_dist
        integer :: i, j
        energy = 0.0_dp
        forces = 0.0_dp
        do i = 2, size(points, 2)
            do j = 1, i - 1
                displ = points(:,i) - points(:,j)
                inv_dist = 1.0_dp / norm2(displ)
                energy = energy + inv_dist
                displ = displ * inv_dist**3
                forces(:,i) = forces(:,i) + displ
                forces(:,j) = forces(:,j) - displ
            end do
        end do
    end subroutine coulomb_forces

    pure subroutine quadratic_line_search(points, energy, &
                                        & step_direction, step_size, workspace)
        real(dp), intent(in) :: points(:,:), energy, step_direction(:,:)
        real(dp), intent(inout) :: step_size
        real(dp), intent(out) :: workspace(:,:)
        real(dp) :: new_energy, newer_energy
        workspace = points + step_size * step_direction
        call normalize_points(workspace)
        new_energy = coulomb_energy(workspace)
        if (new_energy < energy) then
            do
                step_size = 2.0_dp * step_size
                workspace = points + step_size * step_direction
                call normalize_points(workspace)
                newer_energy = coulomb_energy(workspace)
                if (newer_energy >= new_energy) then
                    step_size = 0.25_dp * step_size * &
                        & (4.0_dp * new_energy - newer_energy - 3.0_dp * energy) / &
                        & (2.0_dp * new_energy - newer_energy - energy)
                    return
                end if
                new_energy = newer_energy
            end do
        else
            do
                step_size = 0.5_dp * step_size
                workspace = points + step_size * step_direction
                call normalize_points(workspace)
                newer_energy = coulomb_energy(workspace)
                if (newer_energy < new_energy) then
                    step_size = 0.5_dp * step_size * &
                        & (new_energy - 4.0_dp * newer_energy + 3.0_dp * energy) / &
                        & (new_energy - 2.0_dp * newer_energy + energy)
                    return
                end if
                new_energy = newer_energy
            end do
        end if
    end subroutine quadratic_line_search

end program CoulombOptimizer
