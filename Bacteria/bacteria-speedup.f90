!Leon Wu 01190736
!MATH 96012 Project 3
!This module contains four module variables and two subroutines;
!one of these routines must be developed for this assignment.
!Module variables--
! bm_l, bm_s0, bm_r0, bm_a: the parameters l, s0, r0, and a in the particle dynamics model
! numthreads: The number of threads that should be used in parallel regions within simulate2_omp
!
!Module routines---
! simulate2_f90: Simulate particle dynamics over m trials. Return: all x-y positions at final time
! and alpha at nt+1 times averaged across the m trials.
! simulate2_omp: Same input/output functionality as simulate2.f90 but parallelized with OpenMP

module bmotion
  implicit none
  integer :: numthreads
  real(kind=8) :: bm_l,bm_s0,bm_r0,bm_a
  real(kind=8), parameter :: pi = acos(-1.d0)
  complex(kind=8), parameter :: ii = complex(0.d0,1.d0)
contains

!Compute m particle dynamics simulations using the parameters,bm_l,bm_s0,bm_r0,bm_a.
!Input:
!m: number of simulations
!n: number of particles
!nt: number of time steps
!Output:
! x: x-positions at final time step for all m trials
! y: y-positions at final time step for all m trials
! alpha_ave: alignment parameter at each time step (including initial condition)
! averaged across the m trials
subroutine simulate2_f90(m,n,nt,x,y,alpha_ave)
  implicit none
  integer, intent(in) :: m,n,nt
  real(kind=8), dimension(m,n), intent(out) :: x,y
  real(kind=8), dimension(nt+1), intent(out) :: alpha_ave
  integer :: i1,j1,k1
  real(kind=8), dimension(m,n) :: nn !neighbors
  real(kind=8) :: r0sq !r0^2
  real(kind=8), dimension(m,n) :: phi_init,r_init,theta !used for initial conditions
  real(kind=8), dimension(m,n,n) :: dist2 !distance squared
  real(kind=8), allocatable, dimension(:,:,:) :: temp
  complex(kind=8), dimension(m,n) :: phase
  complex(kind=8), dimension(m,n,nt+1) :: exp_theta,AexpR

  !-------------------------
  !-------------------------

!---Set initial condition and initialize variables----
  allocate(temp(m,n,nt+1))
  call random_number(phi_init)
  call random_number(r_init)
  call random_number(theta)
  call random_number(temp)
  !------------------
  phi_init = phi_init*(2.d0*pi)
  r_init = sqrt(r_init)
  theta = theta*(2.d0*pi) !initial direction of motion
  x = r_init*cos(phi_init)+0.5d0*bm_l !initial positions
  y = r_init*sin(phi_init)+0.5d0*bm_l

  alpha_ave=0.d0
  r0sq = bm_r0*bm_r0
  exp_theta(:,:,1) = exp(ii*theta)

  AexpR = bm_a*exp(ii*temp*2.d0*pi) !noise term
  deallocate(temp)
  nn = 0

!----------------------------------------------
  !Time marching
  do i1 = 2,nt+1

    phase=0.d0
    dist2 = 0.d0

    !Compute distances
    do j1=1,n-1
      do k1 = j1+1,n
        dist2(:,j1,k1) = (x(:,j1)-x(:,k1))**2 + (y(:,j1)-y(:,k1))**2
        where (dist2(:,j1,k1)>r0sq)
          dist2(:,j1,k1)=0
        elsewhere
          dist2(:,j1,k1)=1
        end where
       dist2(:,k1,j1) = dist2(:,j1,k1)
      end do
    end do


    nn = sum(dist2,dim=3)+1

    !Update phase
    phase =  exp_theta(:,:,i1-1) +nn*AexpR(:,:,i1-1)

    do j1=1,m
      phase(j1,:) = phase(j1,:) + matmul(dist2(j1,:,:),exp_theta(j1,:,i1-1))
    end do


    !Update Theta
    theta = atan2(aimag(phase),real(phase))

    !Update X,Y
    x = x + bm_s0*cos(theta)
    y = y + bm_s0*sin(theta)

    x = modulo(x,bm_l)
    y = modulo(y,bm_l)

    exp_theta(:,:,i1) = exp(ii*theta)

  end do

  alpha_ave = (1.d0/dble(m*n))*sum(abs(sum(exp_theta,dim=2)),dim=1)

end subroutine simulate2_f90


!Same functionality as simulate2_f90, but parallelized with OpenMP
!Parallel regions should use numthreads threads.
!Compute m particle dynamics simulations using the parameters,bm_l,bm_s0,bm_r0,bm_a.
!Same functionality as simulate2_f90, but parallelized with OpenMP
!Parallel regions should use numthreads threads.
!Input:
!m: number of simulations
!n: number of particles
!nt: number of time steps
!Output:
! x: x-positions at final time step for all m trials
! y: y-positions at final time step for all m trials
! alpha_ave: alignment parameter at each time step (including initial condition)
! averaged across the m trials
subroutine simulate2_omp(m,n,nt,x,y,alpha_ave)
! I used one parallel region for a do loop across m, which required me to refactor
! some aspects of the code such as the distance calculations and indexing of variables.
! I put this loop outside the time marching loop to avoid repeatedly opening and closing
! parallel regions. I also removed the dimension corresponding to the m simulations and
! Nt time steps from most of the variables, greatly decreasing memory usage. These were declared
! private along with i1, j1 and k1 for inside loops, so that they were independant across simulations
! and therefore threads. I also reduced alpha_ave so that is was updated after every simulation
! and time step.
  use omp_lib
  implicit none
  integer, intent(in) :: m,n,nt
  real(kind=8), dimension(m,n), intent(out) :: x,y
  real(kind=8), dimension(nt+1), intent(out) :: alpha_ave
  integer :: i1,j1,k1,l1
  real(kind=8), dimension(m,n) :: phi_init,r_init !used for initial conditions
  real(kind=8), dimension(n) :: theta, nn, temp
  real(kind=8) :: r0sq !r0^2
  real(kind=8), dimension(n,n) :: dist2 !distance squared
  complex(kind=8), dimension(n) :: phase, exp_theta, AexpR

  call omp_set_num_threads(numthreads) ! Set number of threads for OMP to use
!---Set initial condition and initialize variables (does not need to be parallelized)----

  call random_number(phi_init)
  call random_number(r_init)
  !---------------------
  phi_init = phi_init*(2.d0*pi)
  r_init = sqrt(r_init)
  x = r_init*cos(phi_init)+0.5d0*bm_l !initial positions
  y = r_init*sin(phi_init)+0.5d0*bm_l
!-------------------------------------------------
  alpha_ave=0.d0
  r0sq = bm_r0*bm_r0
  nn = 0
    !----------------------------------------------
    !$OMP parallel do private(dist2, nn, phase, theta, exp_theta, AexpR, temp, i1, j1, k1) reduction(+:alpha_ave)
    !Parallel do loop across simulations
    do l1 = 1, m
        call random_number(theta)
        theta = theta*(2.d0*pi) !initial direction of motion

        exp_theta = exp(ii*theta)
        alpha_ave(1) = alpha_ave(1) + abs(sum(exp_theta)) !Compute first value of alpha_ave

    !Time marching
        do i1 = 2, nt+1

            call random_number(temp)
            AexpR = bm_a*exp(ii*temp*2.d0*pi)

            phase = 0.d0
            dist2 = 0

            do j1 = 1, n-1
                !Compute distances
                  do k1 = j1+1,n
                       dist2(j1,k1) = (x(l1,j1)-x(l1,k1))**2 + (y(l1,j1)-y(l1,k1))**2
                       if (dist2(j1,k1)>r0sq) then
                             dist2(j1,k1)=0
                       else
                             dist2(j1,k1)=1
                       end if
                             dist2(k1,j1) = dist2(j1,k1)
                    end do
             end do

             nn = sum(dist2,dim=2)+1
             phase = exp_theta +nn*AexpR

             phase = phase + matmul(dist2,exp_theta)

             !Update Theta
             theta = atan2(aimag(phase),real(phase))

             !Update X,Y
             x(l1, :) = x(l1, :) + bm_s0*cos(theta)
             y(l1, :) = y(l1, :) + bm_s0*sin(theta)

             x(l1, :) = modulo(x(l1, :),bm_l)
             y(l1, :) = modulo(y(l1, :),bm_l)

             exp_theta = exp(ii*theta)

             !alpha_ave reduction
             alpha_ave(i1) = alpha_ave(i1) + abs(sum(exp_theta))

        end do

    end do
    !$OMP end parallel do

    !Scale alpha_ave by m*n only at the end
    alpha_ave = alpha_ave / dble(m*n)

end subroutine simulate2_omp



subroutine simulate2_omp_full(m,n,nt,x,y,alpha_ave)
! Same as simulate2_omp but outputs x and y for all Nt+1 time steps
  use omp_lib
  implicit none
  integer, intent(in) :: m,n,nt
  real(kind=8), dimension(m,n,Nt+1), intent(out) :: x,y
  real(kind=8), dimension(m,Nt+1), intent(out) :: alpha_ave
  integer :: i1,j1,k1,l1
  real(kind=8), dimension(m,n) :: phi_init,r_init !used for initial conditions
  real(kind=8), dimension(n) :: theta, nn, temp
  !Note: use allocatable arrays where possible to avoid OpenMP memory issues

  real(kind=8) :: r0sq !r0^2
  real(kind=8), dimension(n,n) :: dist2 !distance squared
  complex(kind=8), dimension(n) :: phase, exp_theta, AexpR

  call omp_set_num_threads(numthreads) ! Set number of threads for OMP to use
!---Set initial condition and initialize variables (does not need to be parallelized)----

  call random_number(phi_init)
  call random_number(r_init)
  phi_init = phi_init*(2.d0*pi)
  r_init = sqrt(r_init)
  x(:, :, 1) = r_init*cos(phi_init)+0.5d0*bm_l !initial positions
  y(:, :, 1) = r_init*sin(phi_init)+0.5d0*bm_l
!-------------------------------------------------
  alpha_ave=0.d0
  r0sq = bm_r0*bm_r0
  nn = 0
!---------------------------------------------
    !$OMP parallel do private(dist2, nn, phase, theta, exp_theta, AexpR, temp, i1, j1, k1) reduction(+:alpha_ave)
    do l1 = 1, m
        call random_number(theta)
        theta = theta*(2.d0*pi) !initial direction of motion

        exp_theta = exp(ii*theta)
        alpha_ave(l1, 1) = alpha_ave(l1, 1) + abs(sum(exp_theta)) !Compute first value of alpha_ave

        !Time marching
        do i1 = 2, nt+1

            call random_number(temp)
            AexpR = bm_a*exp(ii*temp*2.d0*pi)

            phase = 0.d0
            dist2 = 0

           do j1 = 1, n-1
                !Compute distances
                do k1 = j1+1,n
                     dist2(j1,k1) = (x(l1,j1,i1-1)-x(l1,k1,i1-1))**2 + (y(l1,j1,i1-1)-y(l1,k1,i1-1))**2
                     if (dist2(j1,k1)>r0sq) then
                         dist2(j1,k1)=0
                     else
                         dist2(j1,k1)=1
                     end if
                       dist2(k1,j1) = dist2(j1,k1)
                    end do
            end do

            nn = sum(dist2,dim=2)+1
            phase = exp_theta +nn*AexpR

            phase = phase + matmul(dist2,exp_theta)

            !Update Theta
            theta = atan2(aimag(phase),real(phase))

            !Update X,Y
              x(l1, :, i1) = x(l1, :, i1-1) + bm_s0*cos(theta)
              y(l1, :, i1) = y(l1, :, i1-1) + bm_s0*sin(theta)

              x(l1, :, i1) = modulo(x(l1, :, i1),bm_l)
              y(l1, :, i1) = modulo(y(l1, :, i1),bm_l)

              exp_theta = exp(ii*theta)

              !alpha_ave reduction
              alpha_ave(l1, i1) = alpha_ave(l1, i1) + abs(sum(exp_theta))

        end do

    end do
    !$OMP end parallel do


    alpha_ave = alpha_ave / dble(n)


end subroutine simulate2_omp_full


end module bmotion
