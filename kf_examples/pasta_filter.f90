module pasta_filter

    use toms675, only: srcf, prmt
    use iso_fortran_env, only: stdout => output_unit

    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: wp = dp

    real(wp), parameter :: ponsart(4) = [0.17096_wp,-0.001309_wp,2.0536e-1_wp,0.004821_wp]

    real(wp), parameter :: temperature = 70.0_wp
    real(wp), parameter :: temperatureK = temperature + 273.15_wp

  interface
    function random_normal() result(ran_norm)
      import wp
      real(wp) :: ran_norm
    end function
  end interface

contains

  function eye(n) result(a)
    integer, intent(in) :: n
    real(wp) :: a(n,n)
    integer :: i
    a = 0
    do i = 1, n
      a(i,i) = 1.0_wp
    enddo
  end function
    
    !> Saturation vapor pressure
    pure function buck(T)
        real(wp), intent(in) :: T   !! Temperature in °C.
        real(wp) :: buck            !! Water vapor pressure in Pascals.

        buck = 0.61121_wp*exp((18.678_wp-T/234.5_wp)*(T/(257.14+T)))*1.e3_wp
    end function

    !> Sorption isotherm
    pure function oswin_sorption_isotherm(X,T,p) result(water_activity)
        
        ! PARAMETERS
        real(wp), intent(in) :: X   !! Moisture content on a dry basis.
        real(wp), intent(in) :: T   !! Temperature in degrees Celsius.
        real(wp), intent(in) :: p(4)   !! A set of parameters, \(k_0,k_1,n_0,n_1\).
        
        ! RETURNS
        real(wp) :: water_activity  !! The water activity (a value between 0 and 1)-

        real(wp) :: a, b

        a = p(1) + p(2)*T
        b = p(3) + p(4)*T
        water_activity = (X/a)**(1.0_wp/b)
        water_activity = water_activity/(1.0_wp + water_activity)
    end function

    subroutine oswin_der(X,T,p,aw,daw)
        real(wp), intent(in) :: X
        real(wp), intent(in) :: T
        real(wp), intent(in) :: p(4)
        real(wp), intent(out) :: aw, daw
        
        real(wp) :: a, b

        a = p(1) + p(2)*T
        b = p(3) + p(4)*T

        daw = (X/a)**(1.0_wp/b)
        aw = daw/(1 + daw)
        daw = daw/(b*X*(1 + daw)**2)
    end subroutine

    pure real(wp) function conc_(psat,RH,t)
        real(wp), intent(in) :: psat    !! Saturation vapor pressure (Pa).
        real(wp), intent(in) :: RH      !! Relative humidity [/] or water activity.
        real(wp), intent(in) :: t       !! Temperature (°K).

        real(wp), parameter :: M_H20 = 0.01801528_wp 
            !! Molecular weight of water (kg/mol).
        real(wp), parameter :: R = 8.314459848_wp    
            !! Universal gas constant (m^3 Pa K^-1 mol^-1).

        conc_ = psat*RH*M_H20/(R*t) ! kg/m^3
    end function

    subroutine trf_prog(t,rh,temp)
        real(wp), intent(in) :: t !! Time in seconds
        real(wp), intent(out) :: rh !! Relative humidity
        real(wp), intent(out) :: temp !! Temperature in Celsius

        real(wp), parameter :: tb(0:10) = [0.0_wp,&
                                     0.23809524_wp,& 
                                     0.33333333_wp,& 
                                     0.44444444_wp,& 
                                     0.53968254_wp,& 
                                     1.25396825_wp,&
                                     1.98412698_wp,& 
                                     2.15873016_wp,& 
                                     4.04761905_wp,& 
                                     4.22222222_wp,& 
                                     4.33333333_wp]
        real(wp), parameter :: rhb(10) = [59,64,64,71,72,73,94,77,94,90]
        real(wp), parameter :: tempb(10) = [49,57,66,74,79,88,81,81,81,47]
        integer :: i
        real(wp) :: thours

        thours = t/3600._wp 

        do i = 1, 10
            if (tb(i-1) <= thours .and. thours < tb(i)) then
                rh = rhb(i)/100._wp
                temp = tempb(i)
                return
            end if
        end do
        ! ambient conditions (one step is missing!)
        rh = 30/100._wp
        temp = 20
    end subroutine

    subroutine perform_measurements(u0,tend,nsteps,nout,D,k,A,rinf,nmeas)
        real(wp), intent(inout) :: u0(:)
        real(wp), intent(in) :: tend
        integer, intent(in) :: nsteps, nout
        real(wp), intent(in) :: D, k, A
        real(wp), intent(in), optional :: rinf
        integer, intent(out) :: nmeas
        
        integer :: n, step, i
        real(wp), pointer :: uold(:)=>null(), unew(:) => null(), utmp(:) => null()
        real(wp) :: dx, dt, fct, aw, rs, psat, tflux, rinf_, rh, ltemp
        real(wp) :: trapz
        character(len=20) :: fname1, fname2, num
        integer :: funit1, funit2

        fname1 = "mprof"

        fname2 = "mflux"
        open(newunit=funit2,file=trim(fname2)//".out")

        n = size(u0)
        dx = 0.0005_wp/real(n-1,wp) ! 1 mm thick pasta
        dt = tend/real(nsteps,wp)
        write(*,*) "dt = ",dt
        write(*,*) "dx = ",dx

        fct = dt*D/dx**2
        if (fct > 0.5_wp) then
            write(*,*) "dt = ",dt
            write(*,*) "dx = ",dx
            write(*,*) "D = ",D
            write(*,*) "Explicit Euler is unstable, dt*D/dx**2 = ", fct
            stop 1
        end if

        allocate(uold(n),unew(n))
        uold = u0
        nmeas = 0

        do step = 0, nsteps
            if (present(rinf)) then
                aw = oswin_sorption_isotherm(uold(n),temperature,ponsart)
                psat = buck(temperature)
                rs = conc_(psat,aw,temperatureK)
                rinf_ = rinf
                rh = 0.2_wp
            else
                call trf_prog(step*dt,rh,ltemp)
                aw = oswin_sorption_isotherm(uold(n),ltemp,ponsart)
                psat = buck(ltemp)
                rs = conc_(psat,aw,ltemp+273.15)
                rinf_ = conc_(psat,rh,ltemp+273.15)
            end if

            call predict(n,unew,uold,fct,dt,dx,rs,rinf_,k)

            if (mod(step,nout)==0) then
                write(num,'(I0.8)') step
                open(newunit=funit1,file=trim(fname1)//trim(num)//".out")
                do i = 1, n
                    write(funit1,*) (i-1)*dx, uold(i)
                end do
                close(funit1)
            ! trapezoidal integration
            trapz = 0
            do i = 1, n-1
                trapz = trapz + dx*0.5*(uold(i+1) + uold(i))
            end do
            trapz = 2*trapz/0.001

            tflux = k*A*(rinf_ - rs)
            write(funit2,*) step, step*dt, uold(n), rs, aw, tflux, rh, ltemp, trapz
            nmeas = nmeas + 1

            end if




            utmp => unew
            unew => uold
            uold => utmp

        end do

        u0 = unew

        deallocate(unew,utmp,uold)
        nullify(unew,utmp,uold)

        close(funit2)

    end subroutine

    subroutine predict(n,unew,uold,fct,dt,dx,rs,rinf,k)
        integer, intent(in) :: n
        real(wp), intent(out) :: unew(n)
        real(wp), intent(in) :: uold(n)
        real(wp), intent(in) :: fct,dt,dx,rs,rinf,k
        integer :: i
        unew(1) = uold(1) + fct*2*(uold(2)-uold(1))
        do i = 2, n-1
            unew(i) = uold(i) + fct*(uold(i+1)-2*uold(i)+uold(i-1))
        end do
        unew(n) = uold(n) + fct*2*(uold(n-1)-uold(n)) + dt*2*k*(rinf - rs)/dx
    end subroutine

    subroutine perform_filtering(nmeas,u0,tend,nsteps,nout,D,k,area,rinf)
        integer, intent(in) :: nmeas
        real(wp), intent(inout) :: u0(:)
        real(wp), intent(in) :: tend
        integer, intent(in) :: nsteps, nout
        real(wp), intent(in) :: D, k, area
        real(wp), intent(in), optional :: rinf

        real(wp), allocatable :: xf(:), xhat(:)
        real(wp), allocatable :: J(:,:), H(:,:), B(:,:), P(:,:), Q(:,:), R(:,:)
        real(wp), allocatable :: zmeas(:,:), zm(:), Kgain(:,:)
        real(wp), allocatable :: wrk(:,:)
        integer :: n, np1, lwrk, nw

        real(wp) :: fct,dt,dx, aw, daw
        real(wp) :: step_dt, uold, rs, tflux, rh, ltemp, trapz, psat, tol, rinf_
        integer :: funit1, funit2, step, i, info
        character(len=20) :: fname1 = "mprof", fname2, num
        logical :: withk, multbq
        real(wp) :: local_temp, local_temp_K


        !
        ! Retrieve measurements
        !
        fname2 = "mflux"
        open(newunit=funit2,file=trim(fname2)//".out",status="old")

        allocate(zmeas(nmeas,2))
        do i = 1, nmeas
            read(funit2,*) step, step_dt, uold, rs, aw, tflux, rh, ltemp, trapz
            zmeas(i,1) = uold + 0.02*random_normal()
            zmeas(i,2) = tflux! + 2.e-6*random_normal()
        end do

        close(funit2)

        print *, zmeas(nmeas,:)

        fname2 = "mflux"
        open(newunit=funit2,file=trim(fname2)//"_filtered.out")


        !
        ! Allocate matrices, state vectors etc
        !
        n = size(u0)


        np1 = n + 1
        nw = np1 + 2 + 3
        lwrk = np1 + 2
        allocate(xf(np1),xhat(np1),P(np1,np1),J(np1,np1),H(2,np1),&
            B(np1,2),Q(2,2),R(2,2),Kgain(np1,2),wrk(lwrk,nw),zm(2))

        !
        ! Kalman filter settings
        !
        multbq = .false.
        tol = epsilon(tol)

        !
        ! Timestep and spatial step
        !
        dx = 0.0005_wp/real(n-1,wp) ! 1 mm thick pasta
        dt = tend/real(nsteps,wp)
        write(*,*) "dt = ",dt
        write(*,*) "dx = ",dx

        fct = dt*D/dx**2
        if (fct > 0.5_wp) then
            write(*,*) "dt = ",dt
            write(*,*) "dx = ",dx
            write(*,*) "D = ",D
            write(*,*) "Explicit Euler is unstable, dt*D/dx**2 = ", fct
            stop 1
        end if

        !
        ! Initialize state vector
        !
        xf(1:n) = u0
        do i = 1, n
            xf(i) = xf(i)! + 0.005_wp*random_normal()
        end do

        ! Estimate initial surface humidity
        if (present(rinf)) then
            call oswin_der(u0(n),temperature,ponsart,aw,daw)
            psat = buck(temperature)
            rs = conc_(psat,aw,temperatureK)
            rinf_ = rinf
            rh = 0.2_wp
        else
            call trf_prog(0.0_wp,rh,ltemp)
            call oswin_der(u0(n),ltemp,ponsart,aw,daw)
            psat = buck(ltemp)
            rs = conc_(psat,aw,ltemp+273.15)
            rinf_ = conc_(psat,rh,ltemp+273.15)
        end if

        ! Initial convection coefficient (supply a known or estimated value)
        xf(np1) = k

        !
        ! Prepare the Jacobian matrix
        !
        J = 0

        ! First comes the identity part
        do i = 1, n
            J(i,i) = 1.0_wp
        end do

        ! Add the finite-difference terms
        do i = 1, n
            if (i == 1) then
                J(i,1) = J(i,1) - fct*2
                J(i,2) = J(i,2) + fct*2
            else if (i == n) then
                J(i,n)   = J(i,n)   - fct*2
                J(i,n-1) = J(i,n-1) + fct*2
            else
                J(i,i-1) = J(i,i-1) + fct
                J(i,i)   = J(i,i)   - fct*2
                J(i,i+1) = J(i,i+1) + fct
            end if
        end do
        

        ! Calculate surface humidity
        if (present(rinf)) then
            call oswin_der(u0(n),temperature,ponsart,aw,daw)
            psat = buck(temperature)
            rs = conc_(psat,aw,temperatureK)
            daw = conc_(psat,daw,temperatureK)
            rinf_ = rinf
            rh = 0.2_wp
        else
            call trf_prog(0.0_wp,rh,ltemp)
            call oswin_der(u0(n),ltemp,ponsart,aw,daw)
            psat = buck(ltemp)
            rs = conc_(psat,aw,ltemp+273.15)
            daw = conc_(psat,daw,ltemp+273.15)
            rinf_ = conc_(psat,rh,ltemp+273.15)
        end if

        J(n,np1) =  dt*2*D*(rinf_ - rs)/dx
        
        ! And last, the equation for the convective coefficient \dot{k} = 0
        J(np1,np1) = 1.0_wp


        !===================
        ! Covariance matrix
        !===================
        ! P = matmul(J,matmul(eye(np1),transpose(J)))
        P = 0.001*eye(np1)

        call dpotrf('L',np1,P,np1,info)
        if (info > 0) then
          write(*,*) "P is not positive definite"
          error stop
        end if

        !====================
        ! State noise matrix
        !====================

        B = 0
        B(1:n,1) = 1
        B(np1,2) = 1

        Q = 0
        Q(1,1) = 0.0000001 ! diffusion noise
        Q(2,2) = 0.00001 ! k noise
        Q(1,2) = 0.00000003
        Q(2,1) = 0.00000003

        call PRMT(Q, 2, 2, 2," Q = ", stdout, 2)
        call dpotrf('L',2,Q,2,info)
        if (info > 0) then
          write(*,*) "Q is not positive definite"
          error stop
        end if

        Q(1,2) = 0.0_wp

        call PRMT(Q, 2, 2, 2," Q = ", stdout, 2)

        Q = matmul(Q,transpose(Q))
        call PRMT(Q, 2, 2, 2," Q = ", stdout, 2)
        stop

        !====================
        ! Observation noise matrix
        !====================
        R = 0
        R(1,1) = 0.04
        R(2,2) = 0.000001
        
        call dpotrf('L',2,R,2,info)
        if (info > 0) then
          write(*,*) "R is not positive definite"
          error stop
        end if

        !
        ! Measurement matrix
        !
        H = 0
        H(1,n) = 1

        do step = 0, nsteps-1

            ! Calculate surface humidity
            if (present(rinf)) then
                call oswin_der(xf(n),temperature,ponsart,aw,daw)
                local_temp = temperature
                local_temp_K = temperatureK
                psat = buck(temperature)
                rs = conc_(psat,aw,temperatureK)
                daw = conc_(psat,daw,temperatureK)
                rinf_ = rinf
                rh = 0.2_wp
            else
                call trf_prog(step*dt,rh,ltemp)
                call oswin_der(xf(n),ltemp,ponsart,aw,daw)
                local_temp = ltemp
                local_temp_K = ltemp+273.15
                psat = buck(ltemp)
                rs = conc_(psat,aw,local_temp_K)
                daw = conc_(psat,daw,local_temp_K)
                rinf_ = conc_(psat,rh,local_temp_K)
            end if

            ! Update Jacobian matrix
            J(n,n) = 1.0_wp - fct*2 + dt*D*2*xf(np1)*daw/dx
            J(n,np1) = dt*D*2*(rinf_ - rs)/dx
            J(np1,np1) = 1.0

            ! Prepare A and Bu
            call predict(n,xhat(1:n),xf(1:n),fct,dt,dx,rs,rinf_,xf(np1))
            xhat(np1) = xf(np1)

            ! Use prediction to evaluate measurement matrix
            ! (Jacobian of measurement function)
            !H(2,n) = daw
            call oswin_der(xhat(n),local_temp,ponsart,aw,daw)
            psat = buck(local_temp)
            rs = conc_(psat,aw,local_temp_K)
            daw = conc_(psat,daw,local_temp_K)
            rinf_ = conc_(psat,rh,local_temp_K)

            H(1,n) = 1.0_wp
            H(2,n) = xhat(np1)*area*daw
            H(2,np1) = area*(rinf_ - rs)

            ! Calculate Kalman gain
            withk = .true.
            call srcf(P,np1,J,np1,B,np1,Q,2,H,2,R,2,np1,2,2,Kgain,np1,wrk,lwrk,multbq,withk,tol)
            if (.not. withk) then
                print *, step, "Nearly singular matrix in srcf"
            end if

            ! Predicted measurements z = h(x)
            zm(1) = xhat(n)                            ! surface concentration
            zm(2) = xhat(np1)*area*(rinf_ - rs) ! total surface flux

            ! Updated state from Kalman gain
            xf = xhat - matmul(Kgain,zm-zmeas(step+1,:))

            ! if (maxval(abs(Kgain)) > 0) then
            !     print *, maxval(Kgain), minval(Kgain)
            ! end if
            ! print *, Kgain(np1,:), dot_product(Kgain(np1,:),zm-zmeas(step+1,:))
            ! print *, step, zm-zmeas(step+1,:), Kgain(n,:), maxval(Kgain(:,1))
            ! xf = xhat

            if (mod(step,nout)==0) then
                write(num,'(I0.8)') step
                open(newunit=funit1,file=trim(fname1)//trim(num)//"_filtered.out")
                do i = 1, n
                    write(funit1,*) (i-1)*dx, xf(i)
                end do
                close(funit1)
            end if

            write(funit2,*) step, step*dt, zm(1), zm(2), zmeas(step+1,1), zmeas(step+1,2), xf(np1)

        end do

        close(funit2)
    end subroutine

end module

program run_pasta_filter

    use pasta_filter
    implicit none

    integer, parameter :: n = 81
    real(wp), allocatable :: u0(:), u(:)

    real(wp) :: tend, D, k, A, rinf
    integer :: nsteps, nout, nmeas

    allocate(u(n))

    u = 0.45_wp
    u0 = u

    tend = 5*60*60._wp
    nsteps = 1000000
    nout = 1000
    D = 1.1e-11_wp
    k = 0.0001_wp
    A = 100._wp

    rinf = conc_(buck(temperature),&
                 oswin_sorption_isotherm(0.2_wp,temperature,ponsart),&
                 temperatureK)

    call perform_measurements(u,tend,nsteps,nout,D,k,A,nmeas=nmeas)

    ! call perform_filtering(nmeas,u0,tend,nsteps,nout,D,k,A)


end program