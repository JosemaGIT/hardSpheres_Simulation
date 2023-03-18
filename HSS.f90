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
    

    !--------------------------------------------------------------------------------------------------
    real, dimension(:,:), allocatable :: r,v

    !r = Position Matrix.
    !v = Velocity Matrix.
    !
    !First Index indicates which particle we are considering.
    !Second Index indicates the (x or y) componet we are considering.
    
    real, dimension(:), allocatable :: t_col
    
    !t_col = The minimum time value for each particle for the next collision.
    !
    !
    !The index indicates which particle we are considering

    real, dimension(2) :: Rij, Vij

    !Rij = r(i,:)-r(j,:)
    !Vij = v(i,:)-v(j,:)
    !
    !
    !Difference between postion and velocity vectors

    integer, dimension(3) :: boundary

    !boundary(1) = Index of the particle reaching the boundary
    !boundary(2) = Type of boundary (x .eq. 1 or y.eq.2 direction)
    !boundary(3) = Type of boundary (Plus.eq.1 or Minus.eq.-1)
    !
    !
    !Identification of the boundary reached.

    real :: L, sigma, t_sym, t, b, Mod_Rij, Mod_Vij, t_bound, t_real_col, t_step

    !L          = Length of the square well side
    !sigma      = Diameter of the spheres
    !t_sym      = Time at which our simulation is running (Set to 0 at the beginning of the simulation)
    !t          = Dummy variable, to calculate a particular time of collision between to particles
    !b          = Parameter that tells us if the collision is going to take place
    !ModRij     = Modulus of a particular Rij
    !ModVij     = Modulus of a particular Vij
    !t_bound    = Minimum time when a particular particle reaches a boundary
    !t_real_col = Minimum time of collision between particles
    !t_step     = Time that will pass to get to the next step
    !
    !
    !Real values that will be used in our simulation

    integer, dimension(:), allocatable :: partner
    !partner(i) = Index of the particle colliding with the i-particle.
    !
    !
    !Used to be able to identificate whiche particles have collided.

    integer :: n,m, i,j,k, updown, part_col, real_partner
    !n=number of particles, m=number of collisions, i,j,k=dummy variables, 
    !part_col is the particle we will take into account, updown is a dummy 
    !variable to identify if you have reached the top or the bottom 
    !of the well

    !n            = Number of particles
    !m            = Number of collisions
    !i,j,k        = Dummy loop variables
    !updown       = Dummy variable to identify some parameters in the boundary collision.
    !part_col     = Index of the particle colliding
    !real_partner = Index of the partner of the particle colliding
    !
    !
    !Integer values that will be used in our simulation.

    !--------------------------------------------------------------------------------------------------

    t_sym=0 !Origin of time at 0

    !! Parameters taken from an external Archive.
    
    open(10,file='parameters.dat')

    read(10,*) sigma
    read(10,*) L
    read(10,*) n
    read(10,*) m

    close(10)

    
    allocate(r(1:n,1:2),v(1:n,1:2),t_col(1:n),partner(1:n)) !We use 'n' to describe the length of our matrixex
    
    !! Initail conditions taken from another external Archive
    
    open(20,file='initial_conditions.dat')

    do i=1,n
        read(20,*) r(i,1), r(i,2), v(i,1), v(i,2)
    end do

    close(20)

    !!      COLLISIONS      !!

    do k=1,m !Repeat for m collisions 
        
        !Reset the variables for each collision
        
        t_col=10000
        partner=0
        t_bound=10000
        boundary=[1,1,1]

        !! Boundary Conditions

        do i=1,n !Repeat for n particles

            do j=1,2 !1 = x-direction ; 2 = y-direction

                if(v(i,j).gt.0) then !Collision will take place at +x (or +y)

                    t=(L-sigma/2-r(i,j))/v(i,j)
                    updown=1
                
                else if(v(i,j).lt.0) then !Collision will take place at -x (or -y)

                    t=(sigma/2-r(i,j))/v(i,j)
                    updown=-1
                    
                end if

                if(t.lt.t_bound) then !Compare to the last particle taken into account

                    t_bound=t
                    boundary(1)=i               !We set the particle reaching the boundary
                    boundary(2)=j               !We set if an x-type or y-type boundary
                    boundary(3)=updown          !We set if it is an + or - type boundary
                
                 end if
            end do
        end do

        !! Interactions between particles

        do i=1,n !Repeat for n particles

            do j=1,n !Collision with he n-1 particles
                
                if (i.eq.j) then !Exclude itself from the calculus
                    cycle
                end if
                
                !Data needed to calculate the time at wich the collision takes place

                Rij=r(i,:)-r(j,:)
                Vij=v(i,:)-v(j,:)
                Mod_Rij=norm2(Rij)
                Mod_Vij=norm2(Vij)
                b=dot_product(Rij,Vij)
                

                if(b.gt.0) then !Exclude particles going away from each other
                    cycle
                else if(b**2-Mod_Vij**2*(Mod_Rij**2-sigma**2).le.0) then !Exclude particles not going to collide
                    cycle
                end if

                !Time needed for this collision to take place

                t=(-b-sqrt(b**2-Mod_Vij**2*(Mod_Rij**2-sigma**2)))/(Mod_Vij**2)     !!!! Tiemes t=NaN, don't know where the error comes from.

                
                if(t.lt.t_col(i)) then !t is less than the last t_col calculated
                
                    t_col(i)=t              !We set the time of the i-particle collision to t
                    partner(i)=j            !We set the i-particle partner to j
                
                 end if
            end do
        end do

        !! Find the real time needed to go to the next step

        !Minimal value of the collision

        t_real_col=minval(t_col) !The value of the time
        
        part_col=minval(minloc(t_col)) !The i-particle we will take into acount 
        !(We use minval here to turn the rank-1 value to rank-0 (Not sure if it can be done otherwise))

        real_partner=partner(part_col) !We find the partner of the 'part_col'-particle

        !Find if in this iteration a collision between particles or a collision with the boundaries takes place

        if (t_real_col.lt.t_bound) then !Collision between particles
            
            t_step=t_real_col        !The next step will take place at t_real_col
            r=r+v*t_step            !We move our particles to their positions at t_step
            
            b=dot_product(r(part_col,:)-r(real_partner,:),v(part_col,:)-v(real_partner,:))                      !Collision parameter

            v(part_col,:)=v(part_col,:) - (r(part_col,:)-r(real_partner,:))*b/(sigma**2)                        !Collision results of v(i,:)
            v(real_partner,:)=v(real_partner,:) + (r(part_col,:)-r(real_partner,:))*b/(sigma**2)                !Collision results of v(j,:)
 
        else if (t_real_col.gt.t_bound) then !Collision with a boundary
            
            t_step=t_bound
            r=r+v*t_step

            v(boundary(1),boundary(2))=-v(boundary(1),boundary(2))

        end if


        t_sym=t_sym+t_step
    end do

    !! We save the results from the last collision in an external file

    open(30,file="results.dat")

    do i=1,n
        write(30,*) r(i,1), r(i,2), v(i,1), v(i,2)
    end do

    close(30)

    open(40,file="velocities.dat")

    do i=1,n
        write(40,*) norm2(v(i,:))
    end do

    close(40)

    deallocate(r,v,t_col,partner) !We deallocate the variables from the ram

end program HSS