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
