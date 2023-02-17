program HSS
!------------------------------------------------------------------------------------------------------
!
!       Hard Spheres Symulation
!       Author: José Manuel Begines Sánchez
!       This program will simulate a system of hard spheres enclosed in an squared infinite potential 
!       well of legth L.
!       Taking the Initial conditions from an external archive it will compute the collisions that  
!       take place inside the well, and write down the position and velocity (state of the sistem)
!       of each particle in each iteration (in another external archive).
!
!       We'll consider equal masses and that the impulse in each collision is normal to the surface of
!       contact between spheres (or spheres and the well boundaries)
!
!------------------------------------------------------------------------------------------------------
    implicit none
    real, dimension(:,:), allocatable :: r,v
    !These will be the position (r) and velocity (v) matrixes. The first index indicates which particle 
    !are we considering while the second indicates the x and y component of the magnitude.

    real, dimension(:), allocatable :: t_col
    !The minimum time value for each particle for the next collision, 

    real, dimension(2) :: Rij, Vij

    !Difference vectors

    integer, dimension(3) :: boundary

    !boundary is an identification of wich particle reaches a boundary and what boundary is 
    !reached 

    real :: L, t_sym, t, sigma,b, Mod_Rij, Mod_Vij, t_cont, t_step, treal_col
    !L is the length of our squared well, t_sym will be the time at wich our simulation 
    !is currently running, t is time dummy variable, sigma the diameters of the spheres,
    !b is a parameter that tells us if the collision is going to take place, ModRij and 
    !ModVij are the modulus of the difference vectors, time of the arrival to the boundary,
    !t_step is the time you have to add for the next step, treal_col is the
    !minimal value of the collision times, 

    integer, dimension(:), allocatable :: partner
    !partner is the partner of collision of each particle

    integer :: n,m,i,j,k,part_col, real_partner, updown
    !n=number of particles, m=number of collisions, i,j,k=dummy variables, 
    !part_col is the particle we will take into account, updown is a dummy 
    !variable to identify if you have reached the top or the bottom 
    !of the well


    t_sym=0 !Origin of time at 0

    !! We'll first take the data from the parameters archive
    open(10,file='parameters.dat')
    read(10,*) sigma
    read(10,*) L
    read(10,*) n
    read(10,*) m
    close(10)

    ! And we'll use the parameter n to describe the legnth of our coordinate matrixes
    allocate(r(1:n,1:2),v(1:n,1:2),t_col(1:n),partner(1:n))
    
    !! We'll set our initial conditions
    
    open(20,file='initial_conditions.dat')

    do i=1,n
        read(20,*) r(i,1), r(i,2), v(i,1), v(i,2)
    end do

    close(20)

    !!We now start the collisions

    do k=1,m !Repeat for m collisions
        
        !Renew the variable t_col and t_cont for the next collision
        t_col=0
        partner=0
        t_cont=0
        boundary=[1,1,1]
        !First let's see at what time does the particles
        do i=1,n !Repeat for n particles
            !4 Boundary Conditions
            do j=1,2
                if(v(i,j)>0) then
                    t=(L-sigma-r(i,j))/v(i,j)
                    updown=1
                else if(v(i,j)<0) then
                    t=(sigma-r(i,j))/v(i,j)
                    updown=-1
                end if

                if(t_cont==0) then
                    t_cont=t
                    boundary(2)=j
                    boundary(3)=updown
                else if(t<t_cont) then
                    t_cont=t
                    boundary(1)=i
                    boundary(2)=j
                    boundary(3)=updown
                end if
            end do
        end do

        !Now we'll se the collision beetween each particle

        do i=1,n !Repeat for n particles
            do j=1,n !Collision with he n-1 particles
                !Exclude itself from the calculus
                if (i==j) then
                    exit
                end if
                
                !Calculate the data for the time

                Rij=r(j,:)-r(i,:)
                Vij=v(j,:)-v(i,:)
                Mod_Rij=norm2(Rij)
                Mod_Vij=norm2(Vij)
                b=dot_product(Rij,Vij)
                
                !Exclude particles going away from each other

                if(b>0) then
                    exit
                end if

                !We calculate the time of the current collision
                t=(-b-sqrt(b**2-Mod_Vij**2*(Mod_Rij**2-sigma**2)))/(Mod_Vij**2)
                
                !If it is the first value it sets t_col(i) to t
                !If t is less than the last t_col(i) it sets it to t
                !Else it leaves t_col(i) as it was
                !It also indicates the partner of collision

                if(t_col(i)==0) then
                    t_col(i)=t
                    partner(i)=j
                else if(t<t_col(i)) then
                    t_col(i)=t
                    partner(i)=j
                end if
            end do
        end do

        treal_col=minval(t_col)
        part_col=minval(minloc(t_col))
        real_partner=partner(part_col)

        if (treal_col<t_cont) then
            t_step=treal_col
            r=r+v*t_step

            v(part_col,:)=v(part_col,:)-dot_product(v(part_col,:)-v(part_col,:),r(part_col,:)-r(part_col,:))
            v(real_partner,:)=v(real_partner,:)-dot_product(v(real_partner,:)-v(real_partner,:),r(real_partner,:)-r(real_partner,:))
            
        else if (treal_col>t_cont) then
            t_step=t_cont
            r=r+v*t_step

            v(boundary(1),boundary(2))=-v(boundary(1),boundary(2))
        end if


        t_sym=t_sym+t_step
    end do

end program HSS