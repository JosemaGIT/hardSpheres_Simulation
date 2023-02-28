program initial_conditions
    implicit none

    real, dimension(:,:), allocatable :: r,v
    !r = Matrix for position initial conditions
    !v = Matrix for velocities initial conditions

    real :: L, sigma, a, b

    integer :: n, m, i, d

!------------------------------------------------------------------------------

    write(*, fmt="(1x,a25)",advance="no") 'Diameter of the spheres: '
    read(*,*) sigma
    
    write(*,*) ''
    
    write(*, fmt="(1x,a32)",advance="no") 'Length of the side of the well: '
    read(*,*) L

    write(*,*) ''

    write(*, fmt="(1x,a31)",advance="no") 'Number of particles (squared): '
    read(*,*) n

    write(*,*) ''

    write(*, fmt="(1x,a22)",advance="no") 'Number of collisions: '
    read(*,*) m

    allocate(r(1:n**2,1:2),v(1:n**2,1:2))

!! We write it in the parameters archive

    open(10,file='parameters.dat')
    
    write(10,*) sigma
    write(10,*) L
    write(10,*) n**2
    write(10,*) m
    
    close(10)

    d=1

!! We generate the initial condition matrix
    do i=1,n**2
        if (mod(i,n).eq.1) then
            r(i,:)=[L/(n+1),d*L/(n+1)]
            d=d+1
        else
            r(i,:)=r(i-1,:)+[L/(n+1),0.0]
        end if
    end do

    do i=1,n**2
        call random_number(a)
        call random_number(b)

        v(i,1)=a-0.5
        v(i,2)=b-0.5
    end do

    open(20,file='initial_conditions.dat')
    do i=1,n**2
        write(20,*) r(i,1), r(i,2), v(i,1), v(i,2)
    end do
    close(20)

    deallocate(r,v)
end program initial_conditions