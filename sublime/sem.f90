!-------------------------------------------------------------------------------
!   module for Synthetic Eddy Method.
!-------------------------------------------------------------------------------
    module var_sem
        implicit none
        private
        public:: sem_ini_eddy,sem_MU2008_update,sem_LD_update
        public:: sem_MU2016_ini,sem_MU2016_update

        integer, parameter:: dpL=  kind(.true.)
        integer, parameter:: dpI=  kind(1)
        integer, parameter:: dpR=  kind(1.0d0)
        real(dpR),parameter:: pi=  4.0d0*atan(1.0d0)
        real(dpR),parameter:: cE=  1.452762113d0

        contains
!       ------------------------------------------------------------------------
!       set the box.
!       ------------------------------------------------------------------------
        subroutine sem_set_box
        implicit none


        return
        end subroutine sem_set_box
!       ------------------------------------------------------------------------
!       initialize the eddies.
!       ------------------------------------------------------------------------
        subroutine sem_ini_eddy(n_eddy,LR,xyz_eddy,e_eddy)
        implicit none
        integer(dpI),intent(in):: n_eddy
        real   (dpR),intent(in):: LR(3,*)
        integer(dpI):: i
        real   (dpR):: e_eddy(3,*),xyz_eddy(3,*),v(6,n_eddy)
                    
        call random_seed()
        call random_number(v)
        do i=1,n_eddy
            xyz_eddy(1:3,i) =  LR(1:3,1)+(LR(1:3,2)-LR(1:3,1))*v(1:3,i)
            e_eddy  (1  ,i) =  sign(1.0d0, v(4,i)-0.5d0)
            e_eddy  (2  ,i) =  sign(1.0d0, v(5,i)-0.5d0)
            e_eddy  (3  ,i) =  sign(1.0d0, v(6,i)-0.5d0)
        end do

        return
        end subroutine sem_ini_eddy
!       ------------------------------------------------------------------------
!       get the velocity fluctuation, Jarrin method.
!       ------------------------------------------------------------------------
        subroutine sem_MU2008_update(dt,U,LR,n_eddy,e_eddy,xyz_eddy,n_f,xyz_f,length_f,stress_f,duvw)
        implicit none
        integer(dpI),intent(in):: n_eddy,n_f
        real   (dpR),intent(in):: dt,U(*),LR(3,*),length_f(*),stress_f(6,*),xyz_f(3,*)
        integer(dpI):: i,k
        real   (dpR):: VOL,e_eddy(3,*),xyz_eddy(3,*),duvw(3,*),xyz(3),uvw(3), &
                    &  dx,dy,dz,A11,A21,A31,A22,A32,A33,f,v(6),C
                    
        C   =  sqrt(1.5d0**3)
!       ------------------------------------------------------------------------
!       update the position of the eddies.
        do i=1,n_eddy
            xyz(1:3)=  xyz_eddy(1:3,i)+dt*U(1:3)

            if(xyz(1) .gt. LR(1,2)) then
                call random_seed()
                call random_number(v)
                xyz(1)      =  LR(1,1)+xyz(1)-LR(1,2)
                xyz(2)      =  LR(2,1)+(LR(2,2)-LR(2,1))*v(2)
                xyz(3)      =  LR(3,1)+(LR(3,2)-LR(3,1))*v(3)
                e_eddy(1,i) =  sign(1.0d0, v(4)-0.5d0)
                e_eddy(2,i) =  sign(1.0d0, v(5)-0.5d0)
                e_eddy(3,i) =  sign(1.0d0, v(6)-0.5d0)
            elseif(xyz(1) .lt. LR(1,1)) then
                call random_seed()
                call random_number(v)
                xyz(1)      =  LR(1,2)+xyz(1)-LR(1,1)
                xyz(2)      =  LR(2,1)+(LR(2,2)-LR(2,1))*v(2)
                xyz(3)      =  LR(3,1)+(LR(3,2)-LR(3,1))*v(3)
                e_eddy(1,i) =  sign(1.0d0, v(4)-0.5d0)
                e_eddy(2,i) =  sign(1.0d0, v(5)-0.5d0)
                e_eddy(3,i) =  sign(1.0d0, v(6)-0.5d0)
            end if

            if(xyz(2) .gt. LR(2,2)) then
                call random_seed()
                call random_number(v)
                xyz(1)      =  LR(1,1)+(LR(1,2)-LR(1,1))*v(1)
                xyz(2)      =  LR(2,1)+xyz(2)-LR(2,2)
                xyz(3)      =  LR(3,1)+(LR(3,2)-LR(3,1))*v(3)
                e_eddy(1,i) =  sign(1.0d0, v(4)-0.5d0)
                e_eddy(2,i) =  sign(1.0d0, v(5)-0.5d0)
                e_eddy(3,i) =  sign(1.0d0, v(6)-0.5d0)
            elseif(xyz(2) .lt. LR(2,1)) then
                call random_seed()
                call random_number(v)
                xyz(1)      =  LR(1,1)+(LR(1,2)-LR(1,1))*v(1)
                xyz(2)      =  LR(2,2)+xyz(2)-LR(2,1)
                xyz(3)      =  LR(3,1)+(LR(3,2)-LR(3,1))*v(3)
                e_eddy(1,i) =  sign(1.0d0, v(4)-0.5d0)
                e_eddy(2,i) =  sign(1.0d0, v(5)-0.5d0)
                e_eddy(3,i) =  sign(1.0d0, v(6)-0.5d0)
            end if

            if(xyz(3) .gt. LR(3,2)) then
                call random_seed()
                call random_number(v)
                xyz(1)      =  LR(1,1)+(LR(1,2)-LR(1,1))*v(1)
                xyz(2)      =  LR(2,1)+(LR(2,2)-LR(2,1))*v(2)
                xyz(3)      =  LR(3,1)+xyz(3)-LR(3,2)
                e_eddy(1,i) =  sign(1.0d0, v(4)-0.5d0)
                e_eddy(2,i) =  sign(1.0d0, v(5)-0.5d0)
                e_eddy(3,i) =  sign(1.0d0, v(6)-0.5d0)
            elseif(xyz(3) .lt. LR(3,1)) then
                call random_seed()
                call random_number(v)
                xyz(1)      =  LR(1,1)+(LR(1,2)-LR(1,1))*v(1)
                xyz(2)      =  LR(2,1)+(LR(2,2)-LR(2,1))*v(2)
                xyz(3)      =  LR(3,2)+xyz(3)-LR(3,1)
                e_eddy(1,i) =  sign(1.0d0, v(4)-0.5d0)
                e_eddy(2,i) =  sign(1.0d0, v(5)-0.5d0)
                e_eddy(3,i) =  sign(1.0d0, v(6)-0.5d0)
            end if

            xyz_eddy(1:3,i) =  xyz(1:3)
        end do
!       update the position of the eddies.
!       ------------------------------------------------------------------------

        VOL = (LR(1,2)-LR(1,1))*(LR(2,2)-LR(2,1))*(LR(3,2)-LR(3,1))

!       ------------------------------------------------------------------------
!       reconstruct the velocity fluctuation.
        do k=1,n_f
            uvw =  0.0d0
            do i=1,n_eddy
                dx  =  abs(xyz_f(1,k)-xyz_eddy(1,i))
                dy  =  abs(xyz_f(2,k)-xyz_eddy(2,i))
                dz  =  abs(xyz_f(3,k)-xyz_eddy(3,i))

                if((dx .le. length_f(k)) .and. (dy .le. length_f(k)) .and. &
                  &(dz .le. length_f(k))) then
                    f   = (1.0d0-dx/length_f(k))*(1.0d0-dy/length_f(k)) &
                        &*(1.0d0-dz/length_f(k))*C
                    uvw(1:3)=  uvw(1:3)+e_eddy(1:3,i)*f &
                            & *sqrt(VOL/(length_f(k)**3*real(n_eddy, dpR)))
                end if
            end do

            A11 =  sqrt(stress_f(1,k)) 
            A21 =  stress_f(2,k) / A11
            A22 =  sqrt(stress_f(4,k) - A21*A21) 
            A31 =  stress_f(3,k) / A11
            A32 = (stress_f(5,k) - A21*A31) / A22
            A33 =  sqrt(stress_f(6,k) - A31*A31 - A32*A32)

!           Reconstruction of the correct velocity fluctuations
            duvw(1,k)   =  A11*uvw(1)
            duvw(2,k)   =  A21*uvw(1) + A22*uvw(2)
            duvw(3,k)   =  A31*uvw(1) + A32*uvw(2) + A33*uvw(3)
        end do
!       reconstruct the velocity fluctuation.
!       ------------------------------------------------------------------------

        return
        end subroutine sem_MU2008_update
!       ------------------------------------------------------------------------
!       generate velocity fluctuations.
!       ------------------------------------------------------------------------
        subroutine sem_LD_update(n_modes,dx,length,eps,nu,velocity,n_vtx,xyz,uvw)
        implicit none
        integer(dpI),intent(in):: n_modes,n_vtx
        real   (dpR),intent(in):: dx,xyz(3,*),length,eps,nu,velocity
        integer(dpI):: ivtx,i,m
        real   (dpR):: uvw(3,*),wave_max,wave_e,wave_min,dw,e,utn,wave(n_modes), &
                    &  kx,ky,kz,v(n_modes,4),fi(n_modes),psi(n_modes),alp(n_modes), &
                    &  theta(n_modes),kxio(n_modes),kyio(n_modes),kzio(n_modes), &
                    &  sxio(n_modes),syio(n_modes),szio(n_modes),rk,arg,tfunk,r, &
                    &  wave_eta

!       highest wave number
        wave_max=  2.0d0*pi/(1.0d0*dx)

!       k_e related to peak energy wave number
        wave_e  =  9.0d0*pi*cE/(55.0d0*length)

!       lowest wave number
        wave_min=  wave_e/2.0d0

!       wavenumber used in the viscous expression in the von Karman spectrum
        wave_eta= (eps/nu**3)**0.25d0

!       wavenumber step
        dw  = (wave_max-wave_min)/real(n_modes, dpR)
        do i=1,n_modes
            wave(i) =  dw*(real(i, dpR)-0.5d0)
        end do

        call random_seed()
        call random_number(v)
        fi (1:n_modes)  =  v(1:n_modes,1)*2.0d0*pi
        psi(1:n_modes)  =  v(1:n_modes,2)*2.0d0*pi
        alp(1:n_modes)  =  v(1:n_modes,3)*2.0d0*pi
        do i=1,n_modes
            theta(i)=  acos(1.0d0-v(i,4)*2.0d0)
        end do

!       wavenumber vector from random angles
        do m=1,n_modes
            kxio(m) =  sin(theta(m))*cos(fi(m))
            kyio(m) =  sin(theta(m))*sin(fi(m))
            kzio(m) =  cos(theta(m))
            sxio(m) =  cos(fi(m))*cos(theta(m))*cos(alp(m))-sin(fi(m))*sin(alp(m))
            syio(m) =  sin(fi(m))*cos(theta(m))*cos(alp(m))+cos(fi(m))*sin(alp(m))
            szio(m) = -sin(theta(m))*cos(alp(m))
        end do

        do ivtx=1,n_vtx
!           initiate turbulent velocities to zero.
            uvw(1:3,ivtx)   =  0.0d0

            do m=1,n_modes
                kx  =  kxio(m)*wave(m)
                ky  =  kyio(m)*wave(m)
                kz  =  kzio(m)*wave(m)
                rk  =  sqrt(kx**2+ky**2+kz**2)
                if(rk .ge. wave_max)    cycle

                arg     =  kx*xyz(1,ivtx)+ky*xyz(2,ivtx)+kz*xyz(3,ivtx)+psi(m)
                tfunk   =  cos(arg)

                r   =  wave(m)/wave_e
                e   =  cE/wave_e*r**4/(1.0d0+r*r)**(17.0d0/6.0d0) &
                    & *exp(-2.0d0*(wave(m)/wave_eta)**2)*velocity**2
                utn =  sqrt(e*dw)

                uvw(1,ivtx) =  uvw(1,ivtx)+2.0d0*utn*tfunk*sxio(m)
                uvw(2,ivtx) =  uvw(2,ivtx)+2.0d0*utn*tfunk*syio(m)
                uvw(3,ivtx) =  uvw(3,ivtx)+2.0d0*utn*tfunk*szio(m)
            end do
        end do

        return
        end subroutine sem_LD_update
!       ------------------------------------------------------------------------
!       get the initial eddy concentration.
!       ------------------------------------------------------------------------
        subroutine sem_MU2016_ini(n_eddy,xyz_eddy,ec,n_f,xyz_f,length_f)
        implicit none
        integer(dpI),intent(in):: n_eddy,n_f
        real   (dpR),intent(in):: length_f(*),xyz_f(3,*)
        integer(dpI):: i,k
        real   (dpR):: xyz_eddy(3,*),ec(*),dx,dy,dz,f,C
                    
        C   =  sqrt(1.5d0**3)
        do k=1,n_f
            ec(k)   =  0.0d0
            do i=1,n_eddy
                dx  =  abs(xyz_f(1,k)-xyz_eddy(1,i))
                dy  =  abs(xyz_f(2,k)-xyz_eddy(2,i))
                dz  =  abs(xyz_f(3,k)-xyz_eddy(3,i))

                if((dx .le. length_f(k)) .and. (dy .le. length_f(k)) .and. &
                  &(dz .le. length_f(k))) then
                    f       = (1.0d0-dx/length_f(k))*(1.0d0-dy/length_f(k)) &
                            &*(1.0d0-dz/length_f(k))*C
                    ec(k)   =  ec(k)+f*f
                end if
            end do
        end do

        return
        end subroutine sem_MU2016_ini
!       ------------------------------------------------------------------------
!       get the velocity fluctuation, Jarrin-2016 method.
!       ------------------------------------------------------------------------
        subroutine sem_MU2016_update(dt,U,LR,n_eddy,e_eddy,xyz_eddy,n_f,xyz_f,length_f,stress_f,duvw,ec)
        implicit none
        integer(dpI),intent(in):: n_eddy,n_f
        real   (dpR),intent(in):: dt,U(*),LR(3,*),length_f(*),stress_f(6,*),xyz_f(3,*)
        integer(dpI):: i,k
        real   (dpR):: VOL,e_eddy(3,*),xyz_eddy(3,*),duvw(3,*),ec(*),xyz(3),uvw(3), &
                    &  dx,dy,dz,A11,A21,A31,A22,A32,A33,f,v(6),C,T,e,a
                    
        C   =  sqrt(1.5d0**3)
        T   =  max(5.0d0*dt, (LR(1,2)-LR(1,1))/U(1))

!       ------------------------------------------------------------------------
!       update the position of the eddies.
        do i=1,n_eddy
            xyz(1:3)=  xyz_eddy(1:3,i)+dt*U(1:3)

            if(xyz(1) .gt. LR(1,2)) then
                call random_seed()
                call random_number(v)
                xyz(1)      =  LR(1,1)+xyz(1)-LR(1,2)
                xyz(2)      =  LR(2,1)+(LR(2,2)-LR(2,1))*v(2)
                xyz(3)      =  LR(3,1)+(LR(3,2)-LR(3,1))*v(3)
                e_eddy(1,i) =  sign(1.0d0, v(4)-0.5d0)
                e_eddy(2,i) =  sign(1.0d0, v(5)-0.5d0)
                e_eddy(3,i) =  sign(1.0d0, v(6)-0.5d0)
            elseif(xyz(1) .lt. LR(1,1)) then
                call random_seed()
                call random_number(v)
                xyz(1)      =  LR(1,2)+xyz(1)-LR(1,1)
                xyz(2)      =  LR(2,1)+(LR(2,2)-LR(2,1))*v(2)
                xyz(3)      =  LR(3,1)+(LR(3,2)-LR(3,1))*v(3)
                e_eddy(1,i) =  sign(1.0d0, v(4)-0.5d0)
                e_eddy(2,i) =  sign(1.0d0, v(5)-0.5d0)
                e_eddy(3,i) =  sign(1.0d0, v(6)-0.5d0)
            end if

            if(xyz(2) .gt. LR(2,2)) then
                call random_seed()
                call random_number(v)
                xyz(1)      =  LR(1,1)+(LR(1,2)-LR(1,1))*v(1)
                xyz(2)      =  LR(2,1)+xyz(2)-LR(2,2)
                xyz(3)      =  LR(3,1)+(LR(3,2)-LR(3,1))*v(3)
                e_eddy(1,i) =  sign(1.0d0, v(4)-0.5d0)
                e_eddy(2,i) =  sign(1.0d0, v(5)-0.5d0)
                e_eddy(3,i) =  sign(1.0d0, v(6)-0.5d0)
            elseif(xyz(2) .lt. LR(2,1)) then
                call random_seed()
                call random_number(v)
                xyz(1)      =  LR(1,1)+(LR(1,2)-LR(1,1))*v(1)
                xyz(2)      =  LR(2,2)+xyz(2)-LR(2,1)
                xyz(3)      =  LR(3,1)+(LR(3,2)-LR(3,1))*v(3)
                e_eddy(1,i) =  sign(1.0d0, v(4)-0.5d0)
                e_eddy(2,i) =  sign(1.0d0, v(5)-0.5d0)
                e_eddy(3,i) =  sign(1.0d0, v(6)-0.5d0)
            end if

            if(xyz(3) .gt. LR(3,2)) then
                call random_seed()
                call random_number(v)
                xyz(1)      =  LR(1,1)+(LR(1,2)-LR(1,1))*v(1)
                xyz(2)      =  LR(2,1)+(LR(2,2)-LR(2,1))*v(2)
                xyz(3)      =  LR(3,1)+xyz(3)-LR(3,2)
                e_eddy(1,i) =  sign(1.0d0, v(4)-0.5d0)
                e_eddy(2,i) =  sign(1.0d0, v(5)-0.5d0)
                e_eddy(3,i) =  sign(1.0d0, v(6)-0.5d0)
            elseif(xyz(3) .lt. LR(3,1)) then
                call random_seed()
                call random_number(v)
                xyz(1)      =  LR(1,1)+(LR(1,2)-LR(1,1))*v(1)
                xyz(2)      =  LR(2,1)+(LR(2,2)-LR(2,1))*v(2)
                xyz(3)      =  LR(3,2)+xyz(3)-LR(3,1)
                e_eddy(1,i) =  sign(1.0d0, v(4)-0.5d0)
                e_eddy(2,i) =  sign(1.0d0, v(5)-0.5d0)
                e_eddy(3,i) =  sign(1.0d0, v(6)-0.5d0)
            end if

            xyz_eddy(1:3,i) =  xyz(1:3)
        end do
!       update the position of the eddies.
!       ------------------------------------------------------------------------

        VOL = (LR(1,2)-LR(1,1))*(LR(2,2)-LR(2,1))*(LR(3,2)-LR(3,1))

!       ------------------------------------------------------------------------
!       reconstruct the velocity fluctuation.
        do k=1,n_f
            uvw =  0.0d0
            e   =  0.0d0
            do i=1,n_eddy
                dx  =  abs(xyz_f(1,k)-xyz_eddy(1,i))
                dy  =  abs(xyz_f(2,k)-xyz_eddy(2,i))
                dz  =  abs(xyz_f(3,k)-xyz_eddy(3,i))

                if((dx .le. length_f(k)) .and. (dy .le. length_f(k)) .and. &
                  &(dz .le. length_f(k))) then
                    f   = (1.0d0-dx/length_f(k))*(1.0d0-dy/length_f(k)) &
                        &*(1.0d0-dz/length_f(k))*C
                    uvw(1:3)=  uvw(1:3)+e_eddy(1:3,i)*f
                    e       =  e+f*f
                end if
            end do
            a       =  dt/T
            ec(k)   =  a*e+(1.0d0-a)*ec(k)
            uvw     =  uvw/sqrt(ec(k))

            A11 =  sqrt(stress_f(1,k)) 
            A21 =  stress_f(2,k) / A11
            A22 =  sqrt(stress_f(4,k) - A21*A21) 
            A31 =  stress_f(3,k) / A11
            A32 = (stress_f(5,k) - A21*A31) / A22
            A33 =  sqrt(stress_f(6,k) - A31*A31 - A32*A32)

!           Reconstruction of the correct velocity fluctuations
            duvw(1,k)   =  A11*uvw(1)
            duvw(2,k)   =  A21*uvw(1) + A22*uvw(2)
            duvw(3,k)   =  A31*uvw(1) + A32*uvw(2) + A33*uvw(3)
        end do
!       reconstruct the velocity fluctuation.
!       ------------------------------------------------------------------------

        return
        end subroutine sem_MU2016_update
    end module var_sem
