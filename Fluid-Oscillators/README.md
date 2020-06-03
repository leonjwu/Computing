# Simulating Fluid Flow PDEs with Finite Differences and Simulating Weakly-Coupled Oscillators on a network

---
## Overview
- Simulated blood flow through a deformed artery
- Designed fast sparse matrix solvers in parallel Fortran code (OMP) to speedup the bottleneck in the finite-difference schemes by 100x
- Simulated Weakly-Coupled Osicalltors on a network using 1D and 2D decompositions for MPI parallel Computing in Fortran.
- Setup communication effectively between multiple processes on a grid using appropriate MPI directives
- Implemented finite difference methods in Python
 

Fortran vs Python Speedup  |  Number of Threads Speedup
:-------------------------:|:-------------------------:
![](https://github.com/leonwu4951/Computing/blob/master/Fluid-Oscillators/speedup.png)  |  ![](https://github.com/leonwu4951/Computing/blob/master/Fluid-Oscillators/threads.png)
---

---
# Relevant Code (Click Below to View Code)

<details><summary><big><big><big><b>Customised Sparse Matrix Solvers in Fortran</b></big></big></big></summary>
<p>

```fortran
!Module for flow simulations of liquid through tube
!This module contains a few module variables (see comments below)
!and four subroutines:
!jacobi: Uses jacobi iteration to compute solution
! to flow through tube
!sgisolve: To be completed. Use sgi method to
! compute flow through tube
!mvec: To be completed; matrix-vector multiplication z = Ay
!mtvec: To be completed; matrix-vector multiplication z = A^T y
module flow
    implicit none
    real(kind=8), parameter :: pi = acos(-1.d0)
    integer :: numthreads !number of threads used in parallel regions
    integer :: fl_kmax=10000 !max number of iterations
    real(kind=8) :: fl_tol=0.00000001d0 !convergence criterion
    real(kind=8), allocatable, dimension(:) :: fl_deltaw !|max change in w| each iteration
    real(kind=8) :: fl_s0=0.1d0 !deformation magnitude
    integer :: iters !Iterations required for convergence (for analysis)

contains
!-----------------------------------------------------
!Solve 2-d tube flow problem with Jacobi iteration
subroutine jacobi(n,w)
    !input  n: number of grid points (n+2 x n+2) grid
    !output w: final velocity field
    !Should also compute fl_deltaw(k): max(|w^k - w^k-1|)
    !A number of module variables can be set in advance.

    integer, intent(in) :: n
    real(kind=8), dimension(0:n+1,0:n+1), intent(out) :: w
    integer :: i1,j1,k1
    real(kind=8) :: del_r,del_t,del_r2,del_t2
    real(kind=8), dimension(0:n+1) :: s_bc,fac_bc
    real(kind=8), dimension(0:n+1,0:n+1) :: r,r2,t,RHS,w0,wnew,fac,fac2,facp,facm

    if (allocated(fl_deltaw)) then
      deallocate(fl_deltaw)
    end if
    allocate(fl_deltaw(fl_kmax))


    !grid--------------
    del_t = 0.5d0*pi/dble(n+1)
    del_r = 1.d0/dble(n+1)
    del_r2 = del_r**2
    del_t2 = del_t**2


    do i1=0,n+1
        r(i1,:) = i1*del_r
    end do

    do j1=0,n+1
        t(:,j1) = j1*del_t
    end do
    !-------------------

    !Update-equation factors------
    r2 = r**2
    fac = 0.5d0/(r2*del_t2 + del_r2)
    facp = r2*del_t2*fac*(1.d0+0.5d0*del_r/r) !alpha_p/gamma
    facm = r2*del_t2*fac*(1.d0-0.5d0*del_r/r) !alpha_m/gamma
    fac2 = del_r2 * fac !beta/gamma
    RHS = fac*(r2*del_r2*del_t2) !1/gamma
    !----------------------------

    !set initial condition/boundary deformation
    w0 = (1.d0-r2)/4.d0
    w = w0
    wnew = w0
    s_bc = fl_s0*exp(-10.d0*((t(0,:)-pi/2.d0)**2))/del_r
    fac_bc = s_bc/(1.d0+s_bc)

    !Jacobi iteration
    do k1=1,fl_kmax
        wnew(1:n,1:n) = RHS(1:n,1:n) + w(2:n+1,1:n)*facp(1:n,1:n) + w(0:n-1,1:n)*facm(1:n,1:n) + &
                                         (w(1:n,0:n-1) + w(1:n,2:n+1))*fac2(1:n,1:n)

        !Apply boundary conditions
        wnew(:,0) = wnew(:,1) !theta=0
        wnew(:,n+1) = wnew(:,n) !theta=pi/2
        wnew(0,:) = wnew(1,:) !r=0
        wnew(n+1,:) = wnew(n,:)*fac_bc !r=1s

        fl_deltaw(k1) = maxval(abs(wnew-w)) !compute relative error

        w=wnew    !update variable
        if (fl_deltaw(k1)<fl_tol) exit !check convergence criterion
        ! if (mod(k1,1000)==0) print *, k1,fl_deltaw(k1)
    end do

    print*, sum(w)
    ! print *, 'k,error=',k1,fl_deltaw(min(k1,fl_kmax))
    print*, 'Jacobi: k=', k1, 'dw_max=', fl_deltaw(min(k1, fl_kmax)), 'sum(theta)', sum(w)
    iters = k1

end subroutine jacobi
!-----------------------------------------------------


!Solve 2-d tube flow problem with sgi method
subroutine sgisolve(n,w)
    use omp_lib
    !input  n: number of grid points (n+2 x n+2) grid
    !output w: final velocity field stored in a column vector
    !Should also compute fl_deltaw(k): max(|w^k - w^k-1|)
    !A number of module variables can be set in advance.
    integer, intent(in) :: n
    real(kind=8), dimension((n+2)*(n+2)), intent(out) :: w
    real(kind=8) :: del_t,del_r,del_r2,del_t2
    !add other variables as needed
    real(kind=8), dimension((n+2)*(n+2)) :: b,d,Ad,ATAd,res,res_new,w_new,res2,Ad2
    real(kind=8), dimension(n+2) :: r,r2,t,fac,fac2,facp,facm,fac_bc,s_bc
    integer :: k1,i1
    real(kind=8) :: kappa,mu,top,top2,bot,temp

    call omp_set_num_threads(numthreads) ! Set number of threads for OMP to use

    if (allocated(fl_deltaw)) then
      deallocate(fl_deltaw)
    end if
    allocate(fl_deltaw(fl_kmax))

    !grid spacings------------
    del_t = 0.5d0*pi/dble(n+1)
    del_r = 1.d0/dble(n+1)
    del_r2 = del_r**2
    del_t2 = del_t**2

    do k1=1, n+2
        r(k1) = (k1-1)*del_r
        t(k1) = (k1-1)*del_t
    end do

    !-------------------
    !Update-equation factors------
    r2 = r**2
    fac = 0.5d0/(r2*del_t2 + del_r2)
    facp = r2*del_t2*fac*(1.d0+0.5d0*del_r/r) !alpha_p/gamma
    facm = r2*del_t2*fac*(1.d0-0.5d0*del_r/r) !alpha_m/gamma
    fac2 = del_r2 * fac !beta/gamma
    s_bc = fl_s0*exp(-10.d0*((t-pi/2.d0)**2))/del_r
    fac_bc = s_bc/(1.d0+s_bc)

    ! Initialize variables for iteration
    b = 0.d0
    do k1=2, n+1
        b((k1-1)*(n+2)+2:k1*(n+2)-1) = -fac(k1)*(r2(k1)*del_r2*del_t2) !-1/gamma
    end do

    ! Initial direction
    call mtvec(n,fac,fac2,facp,facm,fac_bc,b,d)

    ! Initial residual
    res = d

    ! Initial guess
    w = 0.d0

    ! Shifted-gradients iteration
    do k1=1, fl_kmax
        ! Initialize temp variables
        top = 0.d0
        top2 = 0.d0
        bot = 0.d0
        temp = 0.d0

        !OMP versions of mvec and mtvec were tested to be slower (see pdf)
        call mvec(n,fac,fac2,facp,facm,fac_bc,d,Ad)
        call mtvec(n,fac,fac2,facp,facm,fac_bc,Ad,ATAd)

        ! Begin parallel region across each iteration
        !$OMP parallel

        !$OMP do reduction(+:top) reduction(+:bot)
        do i1=1, (n+2)*(n+2)
            top = top + res(i1)**2
            bot = bot + Ad(i1)**2
        end do
        !$OMP end do

        ! kappa (and other variables outside do loops) is set numthreads times
        ! This is ignored since adding a condition to check or similar would cost
        ! time, and this method doesn't lose any time compared to a serial version.
        ! The speedup would also be very small if a method was found to do this.
        ! Ending the parallel region and using !$OMP single to do this was tested
        ! to be slower, so the code is left as is.
        kappa = top/bot

        !$OMP do reduction(+:top2)
        do i1=1, (n+2)*(n+2)
            w_new(i1) = w(i1) + kappa*d(i1)
            res_new(i1) = res(i1) - kappa*ATAd(i1)
            top2 = top2 + res_new(i1)**2
        end do
        !$OMP end do

        mu = top2/top

        !$OMP do
        do i1=1, (n+2)**2
            d(i1) = res_new(i1) + mu*d(i1)
        end do
        !$OMP end do

        !$OMP do reduction(max:temp)
        do i1=1, (n+2)**2
            temp = max(temp, abs(w_new(i1)-w(i1)))
        end do
        !$OMP end do
        fl_deltaw(k1) = temp

        !$OMP do
        do i1=1, (n+2)**2
            w(i1) = w_new(i1)
            res(i1) = res_new(i1)
        end do
        !$OMP end do

        !$OMP end parallel

        ! Check convergence outside parallel region, avoiding having to exit the
        ! parallel region in order to exit, slightly improving performance.
        ! Note that this outputs the updated values of w after converging within
        ! tol (as opposed to the old values of w returned by the non-omp version
        ! of simulate. Thus the algorithm is not perfectly consistant with the old
        ! one although this change is defendable since the solutions are within tolerance.)
        if (fl_deltaw(k1)<fl_tol) exit !check convergence criterion

    end do

    ! print *, 'k,error=',k1,fl_deltaw(min(k1,fl_kmax))
    print*, 'SGI OMP: k=', k1, 'dw_max=', fl_deltaw(min(k1, fl_kmax)), 'sum(theta)', sum(w), 'numthreads', numthreads
    iters = k1
end subroutine sgisolve



!Compute matrix-vector multiplication, z = Ay
subroutine mvec_omp(n,fac,fac2,facp,facm,fac_bc,y,z)
    !input n: grid is (n+2) x (n+2)
    ! fac,fac2,facp,facm,fac_bc: arrays that appear in
    !   discretized equations
    ! y: vector multipled by A
    !output z: result of multiplication Ay
    implicit none
    integer, intent(in) :: n
    real(kind=8), dimension(n+2), intent(in) :: fac,fac2,facp,facm,fac_bc
    real(kind=8), dimension((n+2)*(n+2)), intent(in) :: y
    real(kind=8), dimension((n+2)*(n+2)), intent(out) :: z
    integer :: m, i1, j1, k
    !add other variables as needed

    !Initialize output
    z = 0.d0

    ! Boundary conditions for r=0
    z(1) = -y(1) + y(n+3)
    z(n+2) = -y(n+2) + y(2*(n+2))

    ! Boundary conditions for r=1
    z((n+1)*(n+2)+1) = -fac_bc(1)*y(n*(n+2)+1) + y((n+1)*(n+2)+1)
    z((n+2)*(n+2)) = -fac_bc(n+2)*y((n+1)*(n+2)) + y((n+2)*(n+2))

    !$OMP parallel do private(j1)
    do i1=2, n+1
        ! Boundary conditions for t=0, t=2pi
        z((i1-1)*(n+2) + 1) = -y((i1-1)*(n+2) + 1) + y((i1-1)*(n+2) + 2)
        z(i1*(n+2)) = -y(i1*(n+2)-1) + y(i1*(n+2))

        ! Internal calculations
        do j1 = 2, n+1
            ! Boundary conditions for r=0
            z(j1) = -y(j1) + y(n+2+j1)
            ! Boundary conditions for r=1
            z((n+1)*(n+2)+j1) = -fac_bc(j1)*y(n*(n+2)+j1) + y((n+1)*(n+2)+j1)

            ! Calculate index for m-dim arrays
            k = (i1-1)*(n+2) + j1
            ! Set non-boundary elements
            z(k) = facm(i1)*y(k-(n+2)) + fac2(i1)*(y(k-1) + y(k+1)) - y(k) + facp(i1)*y(k+(n+2))
        end do
    end do
    !$OMP end parallel do


end subroutine mvec_omp


!Compute matrix-vector multiplication, z = A^T y
subroutine mtvec_omp(n,fac,fac2,facp,facm,fac_bc,y,z)
    !input n: grid is (n+2) x (n+2)
    ! fac,fac2,facp,facm,fac_bc: arrays that appear in
    !   discretized equations
    ! y: vector multipled by A^T
    !output z: result of multiplication A^T y
    implicit none
    integer, intent(in) :: n
    real(kind=8), dimension(n+2), intent(in) :: fac,fac2,facp,facm,fac_bc
    real(kind=8), dimension((n+2)*(n+2)), intent(in) :: y
    real(kind=8), dimension((n+2)*(n+2)), intent(out) :: z
    integer :: m, i1, j1, k
    !add other variables as needed
    !Initialize output
    z = 0.d0
    ! k = (i1-1)*(n+2) + j1

    ! First block
    z(1) = -y(1)
    z(2) = -y(2) + facm(2)*y(n+4)
    z(n+1) = -y(n+1) + facm(2)*y(2*n+3)
    z(n+2) = -y(n+2)

    ! Second block
    z(n+3) = y(1) - y(n+3) + fac2(2)*y(n+4)
    z(n+4) = y(2) + y(n+3) - y(n+4) + fac2(2)*y(n+5) + facm(3)*y(2*n+6)
    z(2*n+3) = y(n+1) + fac2(2)*y(2*n+2) - y(2*n+3)- y(2*n+4) + facm(3)*y(3*n+5)
    z(2*n+4) = y(n+2) + fac2(2)*y(2*n+3) + y(2*n+4)

    ! Penultimate block
    z(n*(n+2)+1) = -y(n*(n+2)+1) + fac2(n+1)*y(n*(n+2)+2) - fac_bc(1)*y((n+1)*(n+2)+1)
    z(n*(n+2)+2) = facp(n)*y((n-1)*(n+2)+2) + y(n*(n+2)+1) - y(n*(n+2)+2) + fac2(n+1)*y(n*(n+2)+3) - fac_bc(2)*y((n+1)*(n+2)+2)
    z((n+1)*(n+2)-1) = facp(n)*y(n*(n+2)-1) + fac2(n+1)*y((n+1)*(n+2)-2) - y((n+1)*(n+2)-1) &
        - y((n+1)*(n+2)) - fac_bc(n+1)*y((n+2)*(n+2)-1)
    z((n+1)*(n+2)) = fac2(n+1)*y((n+1)*(n+2)-1) + y((n+1)*(n+2)) - fac_bc(n+2)*y((n+2)*(n+2))

    ! Last block
    z((n+1)*(n+2)+1) = y((n+1)*(n+2)+1)
    z((n+1)*(n+2)+2) = facp(n+1)*y(n*(n+2)+2) + y((n+1)*(n+2)+2)
    z((n+2)*(n+2)-1) = facp(n+1)*y((n+1)*(n+2)-1) + y((n+2)*(n+2)-1)
    z((n+2)*(n+2)) = y((n+2)*(n+2))

    !$OMP parallel

    ! Main blocks
    !$OMP do
    do i1=3, n
        z((i1-1)*(n+2)+1) = - y((i1-1)*(n+2)+1) + fac2(i1)*y((i1-1)*(n+2)+2)
        z((i1-1)*(n+2)+2) = facp(i1-1)*y((i1-2)*(n+2)+2) + facm(i1+1)*y((i1)*(n+2)+2) - y((i1-1)*(n+2)+2) &
            + y((i1-1)*(n+2)+1) + fac2(i1)*y((i1-1)*(n+2)+3)
        z((i1)*(n+2)-1) = facp(i1-1)*y((i1-1)*(n+2)-1) + facm(i1+1)*y((i1+1)*(n+2)-1) + fac2(i1)*y((i1)*(n+2)-2) - y((i1)*(n+2)-1) &
            - y((i1)*(n+2))
        z((i1)*(n+2)) = fac2(i1)*y((i1)*(n+2)-1) + y((i1)*(n+2))
    end do
    !$OMP end do

    ! Parallelize all large calculation blocks
    !$OMP do private(i1)
    do j1=1, n-2
        ! First block
        z(j1+2) = -y(j1+2) + facm(2)*y(j1+n+4)
        ! Second block
        z(n+4+j1) = y(j1+2) + -y(n+4+j1) + fac2(2)*(y(n+3+j1)+ y(n+5+j1)) + facm(3)*y(2*n+6+j1)

        ! Main blocks
        do i1=3, n
            z((i1-1)*(n+2)+2+j1) = facp(i1-1)*y((i1-2)*(n+2)+2+j1) + facm(i1+1)*y((i1)*(n+2)+2+j1) &
                - y((i1-1)*(n+2)+2+j1) + fac2(i1)*(y((i1-1)*(n+2)+1+j1) + y((i1-1)*(n+2)+3+j1))
        end do

        ! Penultimate block
        z(n*(n+2)+2+j1) = facp(n)*y((n-1)*(n+2)+2+j1) - y(n*(n+2)+2+j1) + &
            + fac2(n+1)*(y(n*(n+2)+3+j1) + y(n*(n+2)+1+j1)) - fac_bc(2+j1)*y((n+1)*(n+2)+2+j1)
        ! Last block
        z((n+1)*(n+2)+2+j1) = facp(n+1)*y(n*(n+2)+2+j1) + y((n+1)*(n+2)+2+j1)
    end do
    !$OMP end do
    !$OMP end parallel


end subroutine mtvec_omp



end module flow
```

</p>
</details>


<details><summary><big><big><big><b>MPI Fortran Implementation for Osciallators</b></big></big></big></summary>
<p>

```fortran
!This file contains 1 module, 1 main program and 4 subroutines:
! params: module contain problem parameters and useful constants
! crickets: main program which reads in parameters from data.in
! 	calls simulate, and writes the computed results to output files
! simulate: subroutine for simulating coupled oscillator model
! 	using explicit-Euler time-marching and distributed-memory
! 	parallelization
! RHS: subroutine called by simulate, generates right-hand side
!		of oscillator model equations
! MPE_DECOMP1D: subroutine called by simulate and used to assign
!		oscillators to processes
! random_normal: subroutine called by main program and used to generate
!		natural frequencies, w

!-------------------------------------------------------------
module params
	implicit none
	real(kind=8), parameter :: pi = acos(-1.d0)
	complex(kind=8), parameter :: ii=cmplx(0.0,1.0) !ii = sqrt(-1)
    integer :: ntotal !total number of oscillators,
	real(kind=8) :: c,mu,sigma !coupling coefficient, mean, std for computing omega
	integer :: nlocal_min
	integer, allocatable, dimension(:) :: ai_mod !Coupling ranges for the process
	integer :: rec_upper, rec_lower ! Number of values to receive from upper and lower boundaries
	save
end module params
!-------------------------------

program crickets
    use mpi
    use params
    implicit none
    integer :: i1,j1
    integer :: nt !number of time steps
    real(kind=8) :: dt!time step
    integer :: myid, numprocs, ierr, istart, iend
    real(kind=8), allocatable, dimension(:) :: f0,w,f ! initial phases, frequencies, final phases
    real(kind=8), allocatable, dimension(:) :: r !synchronization parameter

 ! Initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

!gather input
    open(unit=10,file='data.in')
        read(10,*) ntotal !total number of oscillators
        read(10,*) nt !number of time steps
        read(10,*) dt !size of time step
        read(10,*) c ! coupling parameter
        read(10,*) sigma !standard deviation for omega calculation
    close(10)

    allocate(f0(ntotal),f(ntotal),w(ntotal),r(nt))

!generate initial phases
    call random_number(f0)
    f0 = f0*2.d0*pi


!generate frequencies
    mu = 1.d0
    call random_normal(ntotal,w)
    w = sigma*w+mu

!compute min(nlocal)
		nlocal_min = ntotal
		do i1 = 0,numprocs-1
			call mpe_decomp1d(ntotal,numprocs,i1,istart,iend)
			nlocal_min = min(iend-istart+1,nlocal_min)
		end do


!compute solution
    call simulate(MPI_COMM_WORLD,numprocs,ntotal,0.d0,f0,w,dt,nt,f,r)


!output solution (after collecting solution onto process 0 in simulate)
     call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
     if (myid==0) then
        open(unit=11,file='theta.dat')
        do i1=1,ntotal
            write(11,*) f(i1)
        end do
        close(11)

        open(unit=12,file='R.dat')
        do i1=1,nt
	    		write(12,*) r(i1)
				end do
				close(12)
    	end if
    !can be loaded in python, e.g. theta=np.loadtxt('theta.dat')

    call MPI_FINALIZE(ierr)
end program crickets



subroutine simulate(comm,numprocs,n,t0,y0,w,dt,nt,y,r)
    !explicit Euler method, parallelized with mpi
    !input:
    !comm: MPI communicator
    !numprocs: total number of processes
    !n: number of oscillators
    !t0: initial time
    !y0: initial phases of oscillators
    !w: array of frequencies, omega_i
    !dt: time step
    !nt: number of time steps
    !output: y, final solution
    !r: synchronization parameter at each time step
    use mpi
    use params
    implicit none
    integer, intent (in) :: n,nt
    real(kind=8), dimension(n), intent(in) :: y0,w
    real(kind=8), intent(in) :: t0,dt
    real(kind=8), dimension(n), intent(out) :: y
	real(kind=8), dimension(nt), intent(out) :: r
    real(kind=8) :: t
    integer :: i1,k,istart,iend
    integer :: comm,myid,ierr,numprocs
		integer, allocatable, dimension(:) :: seed,ai
		real(kind=8), allocatable, dimension(:) ::  temp
		integer :: nseed,time
		!add other variables as needed
		real(kind=8), allocatable, dimension(:) :: ylocal, Rpart
		integer :: send_upper,send_lower,proc_upper,proc_lower,request,istart_temp,&
			iend_temp,displs(numprocs),recvcounts(numprocs)
		integer, dimension(MPI_STATUS_SIZE) :: status
		complex(kind=8) :: r_partial, r_temp

    call MPI_COMM_RANK(comm, myid, ierr)
    print *, 'start simulate, myid=',myid

    !set initial conditions
    y = y0
    t = t0
    !generate decomposition and allocate sub-domain variables
    call mpe_decomp1d(size(y),numprocs,myid,istart,iend)
    print *, 'istart,iend,threadID=',istart,iend,myid

		!Set coupling ranges, ai
		allocate(ai(iend-istart+1),temp(iend-istart+1))
		call random_seed(size=nseed)
		call system_clock(time)
		allocate(seed(nseed))
		seed = myid+time !remove the "+time" to generate same ai each run
		call random_seed(put=seed)
		call random_number(temp)
		ai = 1 + FLOOR((nlocal_min-1)*temp)

		! Allocate Rpart
		allocate(Rpart(iend-istart+1))

		! Set ai module variable
		allocate(ai_mod(iend-istart+1))
		ai_mod = ai

		! Prepare senders and receivers
		if (myid<numprocs-1) then
			proc_upper = myid+1
		else
			proc_upper = 0
		end if

		if (myid>0) then
			proc_lower = myid-1
		else
			proc_lower = numprocs - 1
		end if

		! Calculate number of values to receive from top and bottom boundary
		do i1=1, iend-istart+1
			rec_upper = max(ai(i1) + i1 - (iend-istart+1), 0, rec_upper)
			rec_lower = max(-(i1-ai(i1)-1), 0, rec_lower)
		end do

		! Send and receive number of values to send from top and bottom boundary
		call MPI_ISEND(rec_upper,1,MPI_INTEGER,proc_upper,0,MPI_COMM_WORLD,request,ierr)
		call MPI_ISEND(rec_lower,1,MPI_INTEGER,proc_lower,1,MPI_COMM_WORLD,request,ierr)
		call MPI_RECV(send_lower,1,MPI_INTEGER,proc_lower,0,MPI_COMM_WORLD,status,ierr)
		call MPI_RECV(send_upper,1,MPI_INTEGER,proc_upper,1,MPI_COMM_WORLD,status,ierr)

		! Allocate ylocal, also making room to store neighbouring values
		allocate(ylocal(iend-istart+1 + rec_upper + rec_lower))

		! Set ylocal, which includes some number of neighbours on each side
		ylocal = y(istart-rec_lower:iend+rec_upper)

		! Send initial data across boundaries
		! Send data at bottom boundary down to next processor
		call MPI_ISEND(ylocal(rec_lower+1:send_lower+rec_lower),send_lower,MPI_DOUBLE_PRECISION,proc_lower,0,MPI_COMM_WORLD,request,ierr)
		! Send data at top boundary up to next processor
		call MPI_ISEND(ylocal(rec_lower+iend-istart+1-send_upper+1:rec_lower+iend-istart+1),&
			send_upper,MPI_DOUBLE_PRECISION,proc_upper,1,MPI_COMM_WORLD,request,ierr)

	!time marching
	do k = 1,nt

		! Receive data from top boundary
		call MPI_RECV(ylocal(rec_lower+iend-istart+2:),rec_upper,MPI_DOUBLE_PRECISION,proc_upper,0,MPI_COMM_WORLD,status,ierr)
		! Receive data from bottom boundary
		call MPI_RECV(ylocal(:rec_lower),rec_lower,MPI_DOUBLE_PRECISION,proc_lower,1,MPI_COMM_WORLD,status,ierr)

		! Calculate Rpart
		call RHS(iend-istart+1,0.d0,w(istart:iend),ylocal,Rpart)

		! Calculate ylocal
		ylocal(rec_lower+1:rec_lower+iend-istart+1) = ylocal(rec_lower+1:rec_lower+iend-istart+1) + dt*Rpart

		! Send data to other processors immediately to minimize waiting times for receiving data
		! Send data at bottom boundary down to next processor
		call MPI_ISEND(ylocal(rec_lower+1:send_lower+rec_lower),send_lower,MPI_DOUBLE_PRECISION,proc_lower,0,MPI_COMM_WORLD,request,ierr)
		! Send data at top boundary up to next processor
		call MPI_ISEND(ylocal(rec_lower+iend-istart+1-send_upper+1:rec_lower+iend-istart+1),&
			send_upper,MPI_DOUBLE_PRECISION,proc_upper,1,MPI_COMM_WORLD,request,ierr)

		! Calculate part of R on each processor
		r_partial = sum(exp(ii*ylocal(rec_lower+1:rec_lower+iend-istart+1)))

		! Reduce r+partial (sum) onto process 0
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		call MPI_REDUCE(r_partial,r_temp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)

		r(k) = 1/dble(n) * abs(r_temp)

    end do

    print *, 'before collection',myid, maxval(abs(ylocal(rec_lower+1:rec_lower+iend-istart+1)))

	! Set displs and recvcounts for MPI_GATHERV
	! This method is only performed once so a more efficient computation for displs
	! and recvcounts was not considered.
	if (myid==0) then
		displs(1) = 0
		do i1 = 0,numprocs-1
			call mpe_decomp1d(ntotal,numprocs,i1,istart_temp,iend_temp)
			displs(i1+2) = displs(i1+1) + iend_temp-istart_temp+1
			recvcounts(i1+1) = iend_temp-istart_temp+1
		end do
	end if

	call MPI_GATHERV(ylocal(rec_lower+1:rec_lower+iend-istart+1),iend-istart+1,MPI_DOUBLE_PRECISION,&
		y,recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    if (myid==0) print *, 'finished',maxval(abs(y))


end subroutine simulate
!-------------------------
subroutine RHS(nn,t,w,f,Rpart)
    !called by simulate
    !Rpart = (1/dt)*(f(t+dt)-f(t))
    use params
    implicit none
    integer, intent(in) :: nn
    real(kind=8), intent(in) :: t
!dimensions of variables below must be added
    real(kind=8), dimension(nn), intent(in) :: w
    real(kind=8), dimension(nn+rec_lower+rec_upper), intent(in) :: f
    real(kind=8), dimension(nn), intent(out) :: Rpart
	integer :: i1, j1
	real(kind=8) :: sum_term

	do i1=1,nn
		sum_term = 0.d0
		do j1 = 1, ai_mod(i1)
			 sum_term = sum_term + sin(f(rec_lower+i1)-f(rec_lower+i1-j1)) &
			  	+ sin(f(rec_lower+i1)-f(rec_lower+i1+j1))
		end do
		Rpart(i1) = w(i1) - dble(c)/dble(ntotal) * sum_term
	end do

end subroutine RHS


!--------------------------------------------------------------------
!  (C) 2001 by Argonne National Laboratory.
!      See COPYRIGHT in online MPE documentation.
!  This file contains a routine for producing a decomposition of a 1-d array
!  when given a number of processors.  It may be used in "direct" product
!  decomposition.  The values returned assume a "global" domain in [1:n]
!
subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )
    implicit none
    integer :: n, numprocs, myid, s, e
    integer :: nlocal
    integer :: deficit

    nlocal  = n / numprocs
    s       = myid * nlocal + 1
    deficit = mod(n,numprocs)
    s       = s + min(myid,deficit)
    if (myid .lt. deficit) then
        nlocal = nlocal + 1
    endif
    e = s + nlocal - 1
    if (e .gt. n .or. myid .eq. numprocs-1) e = n

end subroutine MPE_DECOMP1D

!--------------------------------------------------------------------

subroutine random_normal(n,rn)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

IMPLICIT NONE
integer, intent(in) :: n
real(kind=8), intent(out) :: rn(n)
!     Local variables
integer :: i1
REAL(kind=8)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region
do i1=1,n

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156d0 * (v - 0.5d0)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.d0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
rn(i1) = v/u
end do
RETURN


END subroutine random_normal
```
</p>
</details>

<details><summary><big><big><big><b>Finite Difference Scheme in Python using Fortran Solvers</b></big></big></big></summary>
<p>
	
```python
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import curve_fit
from m1 import flow as fl #assumes p41.f90 has been compiled with: f2py --f90flags='-fopenmp' -c p41.f90 -m m1 -lgomp


def jacobi(n,kmax=10000,tol=1.0e-8,s0=0.1,display=False):
    """ Solve liquid flow model equations with
        jacobi iteration.
        Input:
            n: number of grid points in r and theta
            kmax: max number of iterations
            tol: convergence test parameter
            s0: amplitude of cylinder deformation
            display: if True, plots showing the velocity field and boundary deformation
            are generated
        Output:
            w,deltaw: Final velocity field and |max change in w| each iteration
    """

    #-------------------------------------------
    #Set Numerical parameters and generate grid
    Del_t = 0.5*np.pi/(n+1)
    Del_r = 1.0/(n+1)
    Del_r2 = Del_r**2
    Del_t2 = Del_t**2
    r = np.linspace(0,1,n+2)
    t = np.linspace(0,np.pi/2,n+2) #theta
    tg,rg = np.meshgrid(t,r) # r-theta grid

    #Factors used in update equation (after dividing by gamma)
    rg2 = rg*rg
    fac = 0.5/(rg2*Del_t2 + Del_r2)
    facp = rg2*Del_t2*fac*(1+0.5*Del_r/rg) #alpha_p/gamma
    facm = rg2*Del_t2*fac*(1-0.5*Del_r/rg) #alpha_m/gamma
    fac2 = Del_r2*fac #beta/gamma
    RHS = fac*(rg2*Del_r2*Del_t2) #1/gamma

    #set initial condition/boundary deformation
    w0 = (1-rg**2)/4 #Exact solution when s0=0
    s_bc = s0*np.exp(-10.*((t-np.pi/2)**2))/Del_r
    fac_bc = s_bc/(1+s_bc)

    deltaw = []
    w = w0.copy()
    wnew = w0.copy()

    #Jacobi iteration
    for k in range(kmax):
        #Compute wnew
        wnew[1:-1,1:-1] = RHS[1:-1,1:-1] + w[2:,1:-1]*facp[1:-1,1:-1] + w[:-2,1:-1]*facm[1:-1,1:-1] + (w[1:-1,:-2] + w[1:-1,2:])*fac2[1:-1,1:-1] #Jacobi update

        #Apply boundary conditions
        wnew[:,0] = wnew[:,1] #theta=0
        wnew[:,-1] = wnew[:,-2] #theta=pi/2
        wnew[0,:] = wnew[1,:] #r=0
        wnew[-1,:] = wnew[-2,:]*fac_bc #r=1s

        #Compute delta_p
        deltaw += [np.max(np.abs(w-wnew))]
        w = wnew.copy()
        # if k%1000==0: print("k,dwmax:",k,deltaw[k])
        #check for convergence
        if deltaw[k]<tol:
            # print("Python converged,k=%d,dw_max=%28.16f " %(k,deltaw[k]))
            break

    print(f"Python: k={k}, dw_max={round(deltaw[k])}")
    deltaw = deltaw[:k+1]

    return w,deltaw

if __name__=='__main__':
    #performance()

```
</p>
</details>

