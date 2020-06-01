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
