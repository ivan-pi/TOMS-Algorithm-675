module flight_of_a_ball

  use toms675, only: wp => dp, prmt
  implicit none

  abstract interface
    subroutine rkf45_f(t,y,yp)
      double precision, intent(in) :: t
      double precision, intent(in) :: y(*)
      double precision, intent(out) :: yp(*)
    end subroutine
  end interface

  interface 
    subroutine rkf45(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
      procedure(rkf45_f) :: f
      integer, intent(in) :: neqn
      double precision, intent(inout) :: y(neqn)
      double precision, intent(inout) :: t
      double precision, intent(in) :: tout
      double precision, intent(inout) :: relerr
      double precision, intent(in) :: abserr
      integer, intent(inout) :: iflag
      double precision, intent(inout) :: work(3+6*neqn)
      integer, intent(inout) :: iwork(5)
    end subroutine
  end interface

contains

  subroutine df_dt(t,y,yp)
    real(wp), intent(in) :: t
    real(wp), intent(in) :: y(*)
    real(wp), intent(out) :: yp(*)

    yp(1) = y(3)   ! dx/dt = v_x
    yp(2) = y(4)   ! dy/dt = v_y
    yp(3) = 0.0_wp  ! dv_x/dt = 0
    yp(4) = -9.81_wp ! dv_y/dt = gravity
  end subroutine

  subroutine generate_trajectory(t0,y0,tend,nsteps,t,y)
    real(wp), intent(in) :: t0, y0(4), tend
    integer, intent(in) :: nsteps
    real(wp), intent(out) :: t(nsteps), y(4,nsteps)

    real(wp) :: work(3+6*4), y_(4), t_, dt, tout
    integer :: iwork(5), iflag, i
    real(wp) :: relerr, abserr

    relerr = 1.e-6_wp
    abserr = 1.e-8_wp

    t(1) = t0
    t_ = t0
    y(:,1) = y0
    y_ = y0

    dt = (tend-t0)/real(nsteps-1,wp)

    ! First call
    iflag = 1
    i = 1
    tout = t(1) + dt

    do i = 2, nsteps
      call rkf45(df_dt,4,y_,t_,tout,relerr,abserr,iflag,work,iwork)
      if (iflag == 2 .or. iflag == 3 .or. iflag == 4) then
        tout = t_ + dt
        y(:,i) = y_
        t(i) = t_
      else
        ! if (iflag == 7) then
        !   iflag = 2
        !   cycle
        ! end if
        write(*,*) "Received iflag = ", iflag
        error stop 1
      end if
    end do
  end subroutine

  function eye(n) result(a)
    integer, intent(in) :: n
    real(wp) :: a(n,n)
    integer :: i
    a = 0
    do i = 1, n
      a(i,i) = 1.0_wp
    enddo
  end function

  subroutine filter_trajectory(t,ym,yf)
    real(wp), intent(in) :: t(:), ym(:,:)
    real(wp), intent(out) :: yf(:,:)

    real(wp) :: A(4,4), B(4,4), C(2,4), u(4), ytemp(4)
    real(wp) :: evolution_noise, observation_noise
    real(wp) :: Q(4,4), P(4,4), R(2,2), K(4,2), Pt(4,4)

    real(wp) :: wrk(6,10), tol, ts

    integer :: i, ld, info
    logical :: withk, multbq

    multbq = .false.
    tol = epsilon(tol)

    ld = 4

    ts = (t(size(t)) - t(1))/real(size(t)-1,wp)

    A = eye(4)
    B = eye(4)

    A(1:2,3:4) = A(1:2,3:4) + ts*eye(2)
    B(4,4) = ts

    do i = 1,4
      print *, (A(i,:))
    end do

    do i = 1,4
      print *, (B(i,:))
    end do

    u = 0
    u(4) = -9.81_wp

    C = 0
    C(1:2,1:2) = eye(2)

    P = matmul(A,matmul(eye(4),transpose(A)))

    call dpotrf('L',4,P,ld,info)
    if (info > 0) then
      write(*,*) "P is not positive definite"
      error stop
    end if

    evolution_noise = 0.001_wp
    Q = evolution_noise*eye(4)

    call dpotrf('L',4,Q,ld,info)
    if (info > 0) then
      write(*,*) "Q is not positive definite"
      error stop
    end if

    ! observation_noise = 0.5
    R = eye(2)
    R(1,1) = R(1,1)*0.1
    R(2,2) = R(2,2)*0.5
    
    call dpotrf('L',2,R,2,info)
    if (info > 0) then
      write(*,*) "R is not positive definite"
      error stop
    end if

    ! initial guess for position and velocity
    yf(:,1) = [0._wp,10._wp,20._wp,20._wp]

    do i = 2, size(t)

      withk = .true.
      call srcf(P,ld,A,ld,B,ld,Q,ld,C,2,R,2,4,4,2,K,4,wrk,6,multbq,withk)

      ytemp = matmul(A,yf(:,i-1)) + matmul(B,u)
      yf(:,i) = ytemp - matmul(K,matmul(C,ytemp)-ym(1:2,i-1))
    end do

    ! Calculate final state covariance matrix
      Do i = 1, 4
        pt(1:4,i) = p(i,1:i)
        Call dtrmv('L','N','N',i,p,4,pt(1,i),1)
        pt(i,1:i-1) = pt(1:i-1,i)
      End Do


    write(*,*)
    do i = 1,4
      print *, (pt(i,:))
    end do

  end subroutine

  subroutine ekf_filter_trajectory(t,ym,yf,zf,px,py)
    real(wp), intent(in) :: t(:), ym(:,:), px, py
    real(wp), intent(out) :: yf(:,:), zf(:,:)

    real(wp) :: A(4,4), B(4,4), C(2,4), u(4), ytemp(4), zm(2)
    real(wp) :: evolution_noise, observation_noise
    real(wp) :: Q(4,4), P(4,4), R(2,2), K(4,2), Pt(4,4)

    real(wp) :: wrk(6,10), tol, ts, xx, yy, sqr

    integer :: i, ld, info
    logical :: withk, multbq

    multbq = .false.
    tol = epsilon(tol)

    ld = 4

    ts = (t(size(t)) - t(1))/real(size(t)-1,wp)

    A = eye(4)
    B = eye(4)

    A(1:2,3:4) = A(1:2,3:4) + ts*eye(2)
    B(4,4) = ts

    do i = 1,4
      print *, (A(i,:))
    end do

    do i = 1,4
      print *, (B(i,:))
    end do

    u = 0
    u(4) = -9.81_wp

    P = matmul(A,matmul(eye(4),transpose(A)))

    call dpotrf('L',4,P,ld,info)
    if (info > 0) then
      write(*,*) "P is not positive definite"
      error stop
    end if

    evolution_noise = 0.001_wp
    Q = evolution_noise*eye(4)

    call dpotrf('L',4,Q,ld,info)
    if (info > 0) then
      write(*,*) "Q is not positive definite"
      error stop
    end if

    ! observation_noise = 0.5
    R = eye(2)
    R(1,1) = R(1,1)*0.5
    R(2,2) = R(2,2)*0.001
    
    call dpotrf('L',2,R,2,info)
    if (info > 0) then
      write(*,*) "R is not positive definite"
      error stop
    end if

    ! initial guess for position and velocity
    xx = px + ym(1,1)*cos(ym(2,1))
    yy = py + ym(1,1)*sin(ym(2,1))
    print *, xx, yy
    yf(:,1) = [xx,yy,20*cos(0.25_wp*3.14_wp),20*sin(0.25_wp*3.14_wp)]
    zf(:,1) = ym(:,1)

    do i = 2, size(t)

      ytemp = matmul(A,yf(:,i-1)) + matmul(B,u)
      
      xx = ytemp(1) - px
      yy = ytemp(2) - py

      sqr = hypot(xx,yy)
      C = 0
      C(1,1) = xx/sqr
      C(1,2) = yy/sqr
      C(2,1) = -yy/(xx**2 + yy**2)
      C(2,2) = xx/(xx**2 + yy**2)

      withk = .true.
      call srcf(P,ld,A,ld,B,ld,Q,ld,C,2,R,2,4,4,2,K,4,wrk,6,multbq,withk)

      zm(1) = sqr
      zm(2) = atan(yy,xx)
      zf(:,i) = zm
      yf(:,i) = ytemp - matmul(K,zm-ym(1:2,i-1))
    end do

    ! Calculate final state covariance matrix
      Do i = 1, 4
        pt(1:4,i) = p(i,1:i)
        Call dtrmv('L','N','N',i,p,4,pt(1,i),1)
        pt(i,1:i-1) = pt(1:i-1,i)
      End Do


    write(*,*)
    do i = 1,4
      print *, (pt(i,:))
    end do

  end subroutine

end module

program kf_example

  use toms675, only: srcf, wp => dp
  use flight_of_a_ball, only: generate_trajectory, filter_trajectory, &
    ekf_filter_trajectory
  implicit none

  interface
    function random_normal() result(ran_norm)
      import wp
      real(wp) :: ran_norm
    end function
  end interface

  real(wp) :: t0, tend, y0(4)
  real(wp) :: v0, phi0, px, py
  real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)

  integer, parameter :: nsteps = 101
  real(wp) :: t(nsteps), y(4,nsteps), yf(4,nsteps), la(2,nsteps), laf(2,nsteps)

  integer :: funit, i

  t0 = 0.0_wp
  tend = 3.0_wp

  v0 = 20.0_wp
  phi0 = 0.25_wp*pi

  y0 = [0.0_wp,0.0_wp,v0*cos(phi0),v0*sin(phi0)]

  call generate_trajectory(t0,y0,tend,nsteps,t,y)

  ! Calculate angular distance and angle for EKF example later
  ! and perturb them with Gaussian noise
  px = 20._wp
  py = 0._wp
  do i = 1, nsteps
    la(1,i) = hypot(y(1,i)-px,y(2,i)-py) + 0.5_wp*random_normal()
    la(2,i) = atan(y(2,i)-py,y(1,i)-px)! + 0.1_wp*random_normal()
  end do

  open(newunit=funit,file="ball_path.txt")
  do i = 1, nsteps
    write(funit,*) t(i), y(:,i)
  end do
  close(funit)

  ! Now the path is generated, let's perturb it with some Gaussian noise
  do i = 1, nsteps
    y(1,i) = y(1,i) + 0.1_wp*random_normal()
    y(2,i) = y(2,i) + 0.5_wp*random_normal()
  end do

  open(newunit=funit,file="perturbed_ball_path.txt")
  do i = 1, nsteps
    write(funit,*) t(i), y(1:2,i)
  end do
  close(funit)

  call filter_trajectory(t,y,yf)

  open(newunit=funit,file="filtered_ball_path.txt")
  do i = 1, nsteps
    write(funit,*) t(i), yf(1:2,i)
  end do
  close(funit)

  ! Now let's repeat but this time using the extended Kalman filter algorithm
  ! since the measurement is nonlinear.
  call ekf_filter_trajectory(t,la,yf,laf,px,py)

  open(newunit=funit,file="ekf_filtered_ball_path.txt")
  do i = 1, nsteps
    write(funit,*) t(i), yf(1:2,i), la(1:2,i), laf(1:2,i)
  end do
  close(funit)
end program

! gfortran -o kf_example1 toms675.f90 kf_example1.f90 -L../src/ -ltoms675 -lblas -llapack
