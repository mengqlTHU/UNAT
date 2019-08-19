!-------------------------------------------------------------------------------
!   properties of air, perfect gas.
!-------------------------------------------------------------------------------
    module var_air
        use var_kind_def
        implicit none

        real   (dpR),parameter:: gk     =  1.4d0
        real   (dpR),parameter:: gk1    =  gk-1.0d0
        real   (dpR),parameter:: rr     =  2.87058d2
        real   (dpR),parameter:: t_suth =  110.4d0
        real   (dpR),parameter:: cp     =  gk*rr/gk1
        real   (dpR),parameter:: cv     =     rr/gk1
        real   (dpR),parameter:: prL    =  0.72d0
        real   (dpR),parameter:: prT    =  0.9d0
        contains
!       ------------------------------------------------------------------------
!       get total temperature.
!       ------------------------------------------------------------------------
        function air_get_tt(u) result(tt)
        implicit none
        real(dpR),intent(in):: u(*)
        real(dpR):: tt

        tt  =  u(5)/(rr*u(1))+0.5d0*(u(2)**2+u(3)**2+u(4)**2)/cp

        return
        end function air_get_tt
!       ------------------------------------------------------------------------
!       get total pressure.
!       ------------------------------------------------------------------------
        function air_get_pt(u) result(pt)
        implicit none
        real(dpR),intent(in):: u(*)
        real(dpR):: pt,tt,t

        t   =  u(5)/(rr*u(1))
        tt  =  t+0.5d0*(u(2)**2+u(3)**2+u(4)**2)/cp
        pt  =  u(5)*(t/tt)**(-gk/gk1)

        return
        end function air_get_pt
    end module var_air
!-------------------------------------------------------------------------------
!   Roe flux formulation.
!-------------------------------------------------------------------------------
    subroutine rhs_conv_roe(calD,addD,m)
    use var_kind_def
    use var_air, only: gk,gk1
    use var_temp, only: uLf,uRf,vgf,rhsl,rhsD,fac_1d
    implicit none
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: m
    integer(dpI):: i
    real   (dpR):: uL(5),uR(5),duc(5),aL,aR,hL,hR,unL,unR,d,u(5),h,a,ke,un,vn, &
                &  eig(3),e(3),area,uan,eiga,eigb,ad(5),tv(5),c1,c2,tmp1,tmp2

    vn  =  0.0d0
    do i=1,m
        uL(1:5) =  uLf(1:5,i)
        uR(1:5) =  uRf(1:5,i)

        tmp1    =  0.5d0*(uL(2)*uL(2)+uL(3)*uL(3)+uL(4)*uL(4))
        tmp2    =  0.5d0*(uR(2)*uR(2)+uR(3)*uR(3)+uR(4)*uR(4))
        duc(1)  =  uR(1)-uL(1)
        duc(2:4)=  uR(1)*uR(2:4)-uL(1)*uL(2:4)
        duc(5)  = (uR(5)-uL(5))/gk1+uR(1)*tmp2-uL(1)*tmp1
        aL      =  sqrt(gk*uL(5)/uL(1))
        hL      =  aL*aL/gk1+tmp1
        aR      =  sqrt(gk*uR(5)/uR(1))
        hR      =  aR*aR/gk1+tmp2

        u(1)    =  sqrt(uL(1)*uR(1))
        d       =  sqrt(uR(1)/uL(1))
        tmp1    =  1.0d0/(1.0d0+d)
        u(2:4)  = (uL(2:4)+d*uR(2:4))*tmp1
        h       = (hL     +d*hR     )*tmp1
        ke      =  0.5d0*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
        u(5)    = (h-ke)*gk1*u(1)/gk
        a       =  sqrt(abs(gk*u(5)/u(1)))

        e(1:3)  =  fac_1d(1:3,i)
        area    =  fac_1d(4  ,i)
        vn      =  vgf(1,i)*e(1)+vgf(2,i)*e(2)+vgf(3,i)*e(3)

        unL =  e(1)*uL(2)+e(2)*uL(3)+e(3)*uL(4)-vn
        unR =  e(1)*uR(2)+e(2)*uR(3)+e(3)*uR(4)-vn
        if(.true.) then
            tmp1    =  uL(1)*unL
            tmp2    =  uR(1)*unR
            tv(1)   =  tmp1        +tmp2
            tv(2:4) =  tmp1*uL(2:4)+tmp2*uR(2:4)+e(1:3)*(uL(5)+uR(5))
            tv(5)   =  tmp1*hL     +tmp2*hR     +vn    *(uL(5)+uR(5))
        else
            tv(1)   =  uL(1)*unL+uR(1)*unR
            tv(2:4) =  tv(1)*0.5d0*(uL(2:4)+uR(2:4))+e(1:3)*(uL(5)+uR(5))
            tv(5  ) =  tv(1)*0.5d0*(hL     +hR     )+vn    *(uL(5)+uR(5))
        end if
        if(.not. calD) then
            rhsl(1:5,i) =  tv(1:5)*area*0.5d0
            cycle
        end if

        uan     =  u(2)*e(1)+u(3)*e(2)+u(4)*e(3)
        un      =  uan-vn
        eig(1)  =  abs(un)
        eig(2)  =  abs(un+a)
        eig(3)  =  abs(un-a)
        tmp1    =  2.0d0*abs(unR+aR-unL-aL)
        if(eig(2) .lt. tmp1)    eig(2)  =  0.5d0*(eig(2)*eig(2)/tmp1+tmp1)
        tmp1    =  2.0d0*abs(unR-aR-unL+aL)
        if(eig(3) .lt. tmp1)    eig(3)  =  0.5d0*(eig(3)*eig(3)/tmp1+tmp1)

        eiga    =  0.5d0*(eig(2)+eig(3))-eig(1)
        eigb    =  0.5d0*(eig(2)-eig(3))
        tmp1    = -uan*duc(1)+e(1)*duc(2)+e(2)*duc(3)+e(3)*duc(4)
        tmp2    =  ke*duc(1)-u(2)*duc(2)-u(3)*duc(3)-u(4)*duc(4)+duc(5)
        tmp2    =  tmp2*gk1

        a   =  1.0d0/a
        c1  = (eiga*tmp2*a+eigb*tmp1)*a
        c2  =  eigb*tmp2*a+eiga*tmp1
        ad  =  eig(1)*duc
        ad(1  ) =  ad(1  )+c1
        ad(2:4) =  ad(2:4)+c1*u(2:4)+c2*e(1:3)
        ad(5  ) =  ad(5  )+c1*h     +c2*uan
!       ad(1:5) =  eig(1)*duc(1:5) &
!               &+(eiga*tmp2*a+eigb*tmp1)*a*(/1.0d0, u(2:4), h/) &
!               &+(eigb*tmp2*a+eiga*tmp1)  *(/0.0d0, e(1), e(2), e(3), uan/)

        if(addD) then
            rhsl(1:5,i) =  0.5d0*area*(tv(1:5)-ad(1:5))
        else
            rhsl(1:5,i) =  0.5d0*area*tv(1:5)
            rhsD(1:5,i) =  0.5d0*area*ad(1:5)
        end if
        ! uLf(1,i) = tmp3;
    end do

    return
    end subroutine rhs_conv_roe
!-------------------------------------------------------------------------------
!   Roe flux formulation, reduced dissipation.
!-------------------------------------------------------------------------------
    subroutine rhs_conv_roe_ld(calD,addD,m,sigma)
    use var_kind_def
    use var_air, only: gk,gk1
    use var_temp, only: uLf,uRf,vgf,rhsl,rhsD,fac_1d
    implicit none
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: m
    real   (dpR),intent(in):: sigma(*)
    integer(dpI):: i
    real   (dpR):: uL(5),uR(5),duc(5),aL,aR,hL,hR,unL,unR,d,u(5),h,a,ke,un,vn, &
                &  eig(3),e(3),area,uan,eiga,eigb,ad(5),tv(5),tmp1,tmp2

    vn  =  0.0d0
    do i=1,m
        uL(1:5) =  uLf(1:5,i)
        uR(1:5) =  uRf(1:5,i)

        tmp1    =  0.5d0*(uL(2)*uL(2)+uL(3)*uL(3)+uL(4)*uL(4))
        tmp2    =  0.5d0*(uR(2)*uR(2)+uR(3)*uR(3)+uR(4)*uR(4))
        duc(1)  =  uR(1)-uL(1)
        duc(2:4)=  uR(1)*uR(2:4)-uL(1)*uL(2:4)
        duc(5)  = (uR(5)-uL(5))/gk1+uR(1)*tmp2-uL(1)*tmp1
        hL      =  gk*uL(5)/(gk1*uL(1))+tmp1
        hR      =  gk*uR(5)/(gk1*uR(1))+tmp2

        u(1)    =  sqrt(uL(1)*uR(1))
        d       =  sqrt(uR(1)/uL(1))
        tmp1    =  1.0d0/(1.0d0+d)
        u(2:4)  = (uL(2:4)+d*uR(2:4))*tmp1
        h       = (hL     +d*hR     )*tmp1
        ke      =  0.5d0*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
        u(5)    = (h-ke)*gk1*u(1)/gk
        a       =  sqrt(abs(gk*u(5)/u(1)))

        e(1:3)  =  fac_1d(1:3,i)
        area    =  fac_1d(4  ,i)
        vn      =  vgf(1,i)*e(1)+vgf(2,i)*e(2)+vgf(3,i)*e(3)

        unL =  e(1)*uL(2)+e(2)*uL(3)+e(3)*uL(4)-vn
        unR =  e(1)*uR(2)+e(2)*uR(3)+e(3)*uR(4)-vn
        if(.false.) then
            tmp1    =  uL(1)*unL
            tmp2    =  uR(1)*unR
            tv(1)   =  tmp1        +tmp2
            tv(2:4) =  tmp1*uL(2:4)+tmp2*uR(2:4)+e(1:3)*(uL(5)+uR(5))
            tv(5)   =  tmp1*hL     +tmp2*hR     +vn    *(uL(5)+uR(5))
        else
!           tv(1)   =  uL(1)*unL+uR(1)*unR
            tv(1)   = (uL(1)+uR(1))*(unL+unR)*0.5d0
            tv(2:4) =  tv(1)*0.5d0*(uL(2:4)+uR(2:4))+e(1:3)*(uL(5)+uR(5))
            tv(5  ) =  tv(1)*0.5d0*(hL     +hR     )+vn    *(uL(5)+uR(5))
        end if
        if(.not. calD) then
            rhsl(1:5,i) =  tv(1:5)*area*0.5d0
            cycle
        end if

        aL  =  sqrt(gk*uL(5)/uL(1))
        aR  =  sqrt(gk*uR(5)/uR(1))
        eiga=  a
        if(.false.) then
!           low speed correction.
            tmp1=  sqrt((uL(2)-vgf(1,i))**2+(uL(3)-vgf(2,i))**2+(uL(4)-vgf(3,i))**2)
            tmp2=  sqrt((uR(2)-vgf(1,i))**2+(uR(3)-vgf(2,i))**2+(uR(4)-vgf(3,i))**2)
            d   =  sqrt((u (2)-vgf(1,i))**2+(u (3)-vgf(2,i))**2+(u (4)-vgf(3,i))**2)
            aL  =  min(aL  , tmp1)
            aR  =  min(aR  , tmp2)
            eiga=  min(eiga, d   )
        end if
        tmp1=  2.0d0*abs(unR+aR-unL-aL)
        tmp2=  2.0d0*abs(unR-aR-unL+aL)

        uan     =  u(2)*e(1)+u(3)*e(2)+u(4)*e(3)
        un      =  uan-vn
        eig(1)  =  abs(un)
        eig(2)  =  abs(un+eiga)
        eig(3)  =  abs(un-eiga)
        if(eig(2) .lt. tmp1)    eig(2)  =  0.5d0*(eig(2)*eig(2)/tmp1+tmp1)
        if(eig(3) .lt. tmp2)    eig(3)  =  0.5d0*(eig(3)*eig(3)/tmp2+tmp2)

        eiga    =  0.5d0*(eig(2)+eig(3))-eig(1)
        eigb    =  0.5d0*(eig(2)-eig(3))
        tmp1    = -uan*duc(1)+e(1)*duc(2)+e(2)*duc(3)+e(3)*duc(4)
        tmp2    =  ke*duc(1)-u(2)*duc(2)-u(3)*duc(3)-u(4)*duc(4)+duc(5)
        tmp2    =  tmp2*gk1

        a   =  1.0d0/a
        ad(1:5) =  eig(1)*duc(1:5) &
                &+(eiga*tmp2*a+eigb*tmp1)*a*(/1.0d0, u(2:4), h/) &
                &+(eigb*tmp2*a+eiga*tmp1)  *(/0.0d0, e(1), e(2), e(3), uan/)

        if(addD) then
            rhsl(1:5,i) =  0.5d0*area*(tv(1:5)-sigma(i)*ad(1:5))
        else
            rhsl(1:5,i) =  0.5d0*area*tv(1:5)
            rhsD(1:5,i) =  0.5d0*area*ad(1:5)*sigma(i)
        end if
    end do

    return
    end subroutine rhs_conv_roe_ld
!-------------------------------------------------------------------------------
!   HLLC flux formulation.
!-------------------------------------------------------------------------------
    subroutine rhs_conv_hllc(calD,addD,m)
    use var_kind_def
    use var_air, only: gk,gk1
    use var_temp
    implicit none
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: m
    integer(dpI):: i
    real   (dpR):: uL(5),uR(5),n(3),S,rh,u(3),p,d,h,a,ke,un,vn,sL,sR,sm,ps, &
                &  aL,aR,unL,unR,hL,hR,rheL,rheR,us(5),f(5),rtm1,rtm2

    vn  =  0.0d0
    do i=1,M
        n(1:3)  =  fac_1d(1:3,i)
        S       =  fac_1d(4  ,i)
        uL(1:5) =  uLf   (1:5,i)
        uR(1:5) =  uRf   (1:5,i)

        vn  =  vgf(1,i)*n(1)+vgf(2,i)*n(2)+vgf(3,i)*n(3)
        unL =  uL(2)*n(1)+uL(3)*n(2)+uL(4)*n(3)-vn
        unR =  uR(2)*n(1)+uR(3)*n(2)+uR(4)*n(3)-vn
        aL  =  sqrt(gk*uL(5)/uL(1))
        ke  =  0.5d0*(uL(2)*uL(2)+uL(3)*uL(3)+uL(4)*uL(4))
        hL  =  aL*aL/gk1+ke
        rheL=  uL(5)/gk1+uL(1)*ke
        aR  =  sqrt(gk*uR(5)/uR(1))
        ke  =  0.5d0*(uR(2)*uR(2)+uR(3)*uR(3)+uR(4)*uR(4))
        hR  =  aR*aR/gk1+ke
        rheR=  uR(5)/gk1+uR(1)*ke
        if(.not. calD) then
            sL      =  uL(1)*unL
            sR      =  uR(1)*unR
            ps      =  uL(5)+uR(5)
            f(1  )  =  sL        +sR
            f(2:4)  =  sL*uL(2:4)+sR*uR(2:4)+ps*n(1:3)
            f(5  )  =  sL*hL     +sR*hR     +ps*vn
            rhsl(1:5,i) =  f(1:5)*S*0.5d0
            cycle
        end if

        rh      =  sqrt(uL(1)*uR(1))
        d       =  sqrt(uR(1)/uL(1))
        rtm1    =  1.0d0/(1.0d0+d)
        u(1:3)  = (uL(2:4)+d*uR(2:4))*rtm1
        h       = (hL+d*hR)*rtm1
        ke      =  0.5d0*(u(1)**2+u(2)**2+u(3)**2)
        p       = (h-ke)*gk1*rh/gk
        a       =  sqrt(abs(gk*p/rh))

        un  =  u(1)*n(1)+u(2)*n(2)+u(3)*n(3)-vn
        sL  =  min(unL-aL, un-a)
        sR  =  max(unR+aR, un+a)
        rtm1=  sL-unL
        rtm2=  sR-unR
        sm  = (uR(1)*unR*rtm2-uL(1)*unL*rtm1+uL(5)-uR(5))/(uR(1)*rtm2-uL(1)*rtm1)

        if    (sL .gt. 0.0d0) then
            f(1  )  =  uL(1)*unL
            f(2:4)  =  f(1)*uL(2:4)+uL(5)*n(1:3)
            f(5  )  =  f(1)*hL     +uL(5)*vn
        elseif(sR .lt. 0.0d0) then
            f(1  )  =  uR(1)*unR
            f(2:4)  =  f(1)*uR(2:4)+uR(5)*n(1:3)
            f(5  )  =  f(1)*hR     +uR(5)*vn
        elseif(sm .gt. 0.0d0) then
            ps      =  uL(1)*(unL-sL)*(unL-sm)+uL(5)
            rtm2    =  ps-uL(5)
            us(1  ) =  rtm1*uL(1)
            us(2:4) =  rtm1*uL(1)*uL(2:4)+rtm2*n(1:3)
            us(5  ) =  rtm1*rheL-uL(5)*unL+ps*sm
            us      =  us/(sL-sm)
            f(1  )  =  sm*us(1)
            f(2:4)  =  sm*us(2:4)+ps*n(1:3)
            f(5  )  =  sm*us(5)  +ps*(sm+vn)
        else
            ps      =  uR(1)*(unR-sR)*(unR-sm)+uR(5)
            rtm1    =  ps-uR(5)
            us(1)   =  rtm2*uR(1)
            us(2:4) =  rtm2*uR(1)*uR(2:4)+rtm1*n(1:3)
            us(5)   =  rtm2*rheR-uR(5)*unR+ps*sm
            us      =  us/(sR-sm)
            f(1)    =  sm*us(1)
            f(2:4)  =  sm*us(2:4)+ps*n(1:3)
            f(5)    =  sm*us(5)+ps*(sm+vn)
        end if

        if(addD) then
            rhsl(1:5,i) =  f(1:5)*S
        else
            sL      =  uL(1)*unL*0.5d0
            sR      =  uR(1)*unR*0.5d0
            ps      = (uL(5)+uR(5))*0.5d0
            us(1)   =  sL+sR
            us(2:4) =  sL*uL(2:4)+sR*uR(2:4)+ps*n(1:3)
            us(5)   =  sL*hL     +sR*hR     +ps*vn
            rhsl(1:5,i) =  S* us(1:5)
            rhsD(1:5,i) =  S*(us(1:5)-f(1:5))
        end if
    end do

    return
    end subroutine rhs_conv_hllc
!-------------------------------------------------------------------------------
!   AUSM flux formulation.
!-------------------------------------------------------------------------------
    subroutine rhs_conv_ausm(calD,addD,m)
    use var_kind_def
    use var_air, only: gk,gk1
    use var_temp, only: uLf,uRf,vgf,rhsl,rhsD,fac_1d
    implicit none
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: m
    integer(dpI):: i
    real   (dpR):: e(3),S,aL,aR,unL,unR,vn,maL,maR,hL,hR,mL,pL,mR,pR,mf,p,pm, &
                &  uL(5),uR(5),asL,asR,am,fL(5),fR(5),f(5),ua

    vn  =  0.0d0
    do i=1,M
        e(1:3)  =  fac_1d(1:3,i)
        S       =  fac_1d(4  ,i)
        uL(1:5) =  uLf(1:5,i)
        uR(1:5) =  uRf(1:5,i)
        aL      =  sqrt(gk*uL(5)/uL(1))
        aR      =  sqrt(gk*uR(5)/uR(1))
        hL      =  aL*aL/gk1+0.5d0*(uL(2)*uL(2)+uL(3)*uL(3)+uL(4)*uL(4))
        hR      =  aR*aR/gk1+0.5d0*(uR(2)*uR(2)+uR(3)*uR(3)+uR(4)*uR(4))
        pm      =  0.5d0*(uL(5)+uR(5))

        vn  =  vgf(1,i)*e(1)+vgf(2,i)*e(2)+vgf(3,i)*e(3)
        unL =  uL(2)*e(1)+uL(3)*e(2)+uL(4)*e(3)-vn
        unR =  uR(2)*e(1)+uR(3)*e(2)+uR(4)*e(3)-vn
        if(.not. calD) then
            mL      =  0.5d0*uL(1)*unL
            mR      =  0.5d0*uR(1)*unR
            f(1  )  =  mL+mR
            f(2:4)  =  mL*uL(2:4)+mR*uR(2:4)+pm*e(1:3)
            f(5  )  =  mL*hL     +mR*hR     +pm*vn
            rhsl(1:5,i) =  S*f(1:5)
            cycle
        end if

        asL     =  sqrt(2.0d0*hL*gk1/(gk+1.0d0))
        asR     =  sqrt(2.0d0*hR*gk1/(gk+1.0d0))
        fL(1:5) = (/1.0d0, uL(2:4), hL/)*uL(1)
        fR(1:5) = (/1.0d0, uR(2:4), hR/)*uR(1)
        ua      = (uL(2)-vgf(1,i))**2+(uL(3)-vgf(2,i))**2+(uL(4)-vgf(3,i))**2 &
                &+(uR(2)-vgf(1,i))**2+(uR(3)-vgf(2,i))**2+(uR(4)-vgf(3,i))**2
        ua      =  sqrt(0.5d0*ua)

        am  =  min(asL*asL/max(asL, unL), asR*asR/max(asR, -unR))
        maL =  unL/am
        maR =  unR/am
        if    (maL .ge. 1.0d0) then
            mL  =  maL
            pL  =  1.0d0
        elseif(maL .le. -1.0d0) then
            mL  =  0.0d0
            pL  =  0.0d0
        else
            mL  =  0.25d0*(MaL+1.0d0)**2+0.125d0*(MaL*MaL-1.0d0)**2
!           pL  =  0.25d0*(MaL+1.0d0)**2*(2.0d0-MaL)+0.1875d0*MaL*(MaL**2-1.0d0)**2
            pL  =  0.25d0*(MaL+1.0d0)**2*(2.0d0-MaL)
        end if

        if    (maR .le. -1.0d0) then
            mR  =  maR
            pR  =  1.0d0
        elseif(maR .ge. 1.0d0) then
            mR  =  0.0d0
            pR  =  0.0d0
        else
            mR  = -0.25d0*(MaR-1.0d0)**2-0.125d0*(MaR*MaR-1.0d0)**2
!           pR  =  0.25d0*(MaR-1.0d0)**2*(2.0d0+MaR)-0.1875d0*MaR*(MaR**2-1.0d0)**2
            pR  =  0.25d0*(MaR-1.0d0)**2*(2.0d0+MaR)
        end if
        mf  = (mL+mR)*am

        if(.true.) then
            p   =  pm+0.5d0*(pL-pR)*(uL(5)-uR(5)) &
                &+(pL+pR-1.0d0)*ua*0.25d0*(uL(1)+uR(1))*(aL+aR)
        else
!           p   =  pL*uL(5)+pR*uR(5)
!           p   =  pL*uL(5)+pR*uR(5)-pL*pR*(uL(1)+uR(1))*am*(unR-unL)
            p   =  pL*uL(5)+pR*uR(5)-pL*pR*(uL(1)+uR(1))*ua*(unR-unL)
        end if
        mL      =  max(mf, 0.0d0)
        mR      =  min(mf ,0.0d0)
        f(1)    =  mL*fL(1)  +mR*fR(1)
        f(2:4)  =  mL*fL(2:4)+mR*fR(2:4)+p*e(1:3)
        f(5)    =  mL*fL(5)  +mR*fR(5)  +p*vn
        if(addD) then
            rhsl(1:5,i) =  S*f(1:5)
        else
            mL          =  0.5d0*uL(1)*unL
            mR          =  0.5d0*uR(1)*unR
            rhsl(1  ,i) =  S*(mL        +mR)
            rhsl(2:4,i) =  S*(mL*uL(2:4)+mR*uR(2:4)+pm*e(1:3))
            rhsl(5  ,i) =  S*(mL*hL     +mR*hR     +pm*vn)
            rhsD(1:5,i) =  rhsl(1:5,i)-S*f(1:5)
        end if
    end do

    return
    end subroutine rhs_conv_ausm
!-------------------------------------------------------------------------------
!   LF flux formulation.
!-------------------------------------------------------------------------------
    subroutine rhs_conv_lf(calD,addD,m)
    use var_kind_def
    use var_air, only: gk,gk1
    use var_temp, only: uLf,uRf,vgf,rhsl,rhsD,fac_1d
    implicit none
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: m
    integer(dpI):: i
    real   (dpR):: uL(5),uR(5),duc(5),aL,aR,hL,hR,unL,unR,vn,e(3),area,ad(5),tv(5),tmp1,tmp2

    do i=1,m
        uL(1:5) =  uLf(1:5,i)
        uR(1:5) =  uRf(1:5,i)

        tmp1    =  0.5d0*(uL(2)*uL(2)+uL(3)*uL(3)+uL(4)*uL(4))
        tmp2    =  0.5d0*(uR(2)*uR(2)+uR(3)*uR(3)+uR(4)*uR(4))
        duc(1)  =  uR(1)-uL(1)
        duc(2:4)=  uR(1)*uR(2:4)-uL(1)*uL(2:4)
        duc(5)  = (uR(5)-uL(5))/gk1+uR(1)*tmp2-uL(1)*tmp1
        aL      =  sqrt(gk*uL(5)/uL(1))
        hL      =  aL*aL/gk1+tmp1
        aR      =  sqrt(gk*uR(5)/uR(1))
        hR      =  aR*aR/gk1+tmp2

        e(1:3)  =  fac_1d(1:3,i)
        area    =  fac_1d(4  ,i)
        vn      =  vgf(1,i)*e(1)+vgf(2,i)*e(2)+vgf(3,i)*e(3)

        unL =  e(1)*uL(2)+e(2)*uL(3)+e(3)*uL(4)-vn
        unR =  e(1)*uR(2)+e(2)*uR(3)+e(3)*uR(4)-vn
        if(.true.) then
            tmp1    =  uL(1)*unL
            tmp2    =  uR(1)*unR
            tv(1)   =  tmp1        +tmp2
            tv(2:4) =  tmp1*uL(2:4)+tmp2*uR(2:4)+e(1:3)*(uL(5)+uR(5))
            tv(5)   =  tmp1*hL     +tmp2*hR     +vn    *(uL(5)+uR(5))
        else
            tv(1)   =  uL(1)*unL+uR(1)*unR
            tv(2:4) =  tv(1)*0.5d0*(uL(2:4)+uR(2:4))+e(1:3)*(uL(5)+uR(5))
            tv(5  ) =  tv(1)*0.5d0*(hL     +hR     )+vn    *(uL(5)+uR(5))
        end if
        if(.not. calD) then
            rhsl(1:5,i) =  tv(1:5)*area*0.5d0
            cycle
        end if

        ad  =  max(abs(unL)+aL, abs(unR)+aR)*duc
        if(addD) then
            rhsl(1:5,i) =  0.5d0*area*(tv(1:5)-ad(1:5))
        else
            rhsl(1:5,i) =  0.5d0*area*tv(1:5)
            rhsD(1:5,i) =  0.5d0*area*ad(1:5)
        end if
    end do

    return
    end subroutine rhs_conv_lf
!-------------------------------------------------------------------------------
!   primitive variables to conservative variables.
!-------------------------------------------------------------------------------
    subroutine u_to_akh(N,u,akh)
    use var_kind_def
    use var_air, only: gk,gk1
    implicit none
    integer(dpI),intent(in):: N
    integer(dpI):: i
    real   (dpR):: u(5,*),akh(3,*)

    do i=1,N
        akh(1,i)=  sqrt(gk*u(5,i)/u(1,i))
        akh(2,i)=  0.5d0*(u(2,i)*u(2,i)+u(3,i)*u(3,i)+u(4,i)*u(4,i))
        akh(3,i)=  akh(1,i)*akh(1,i)/gk1+akh(2,i)
    end do

    return
    end subroutine u_to_akh
!-------------------------------------------------------------------------------
!   primitive variables to mu.
!-------------------------------------------------------------------------------
    subroutine u_to_mu(N,u,mu)
    use var_kind_def
    use var_air, only: rr,t_Suth
    implicit none
    integer(dpI),intent(in):: N
    integer(dpI):: i
    real   (dpR):: u(5,*),mu(*),t

    do i=1,N
        t       =  u(5,i)/(rr*u(1,i))
        mu(i)   =  1.461d-6*sqrt(t**3)/(t+t_suth)
    end do

    return
    end subroutine u_to_mu
!-------------------------------------------------------------------------------
!   u-->uc.
!-------------------------------------------------------------------------------
    subroutine u_to_uc(N,u_in,uc_out,t_out,mu_out)
    use var_kind_def
    use var_air, only: rr,gk,gk1,t_Suth
    use var_turb, only: is_tur_cal
    implicit none
    integer(dpI),intent(in):: N
    real   (dpR),intent(in):: u_in(5,*)
    integer(dpI):: i
    real   (dpR):: uc_out(5,*),t_out(*),mu_out(*),a,ke,h

    if(N .le. 0)    return
    do i=1,N
        a               =  sqrt(abs(gk*u_in(5,i)/u_in(1,i)))
        ke              =  0.5d0*(u_in(2,i)**2+u_in(3,i)**2+u_in(4,i)**2)
        h               =  a*a/gk1+ke
        uc_out(1  ,i)   =  u_in(1,i)
        uc_out(2:4,i)   =  u_in(1,i)*u_in(2:4,i)
        uc_out(5  ,i)   =  u_in(5,i)/gk1+ke*u_in(1,i)
        t_out (    i)   =  u_in(5,i)/(rr*u_in(1,i))
        if(is_tur_cal) then
            mu_out(2*i-1)   =  1.461d-6*sqrt(t_out(i)**3)/(t_out(i)+t_suth)
        else
            mu_out(i    )   =  1.461d-6*sqrt(t_out(i)**3)/(t_out(i)+t_suth)
        end if
    end do

    return
    end subroutine u_to_uc
!-------------------------------------------------------------------------------
!   transfrom uc into u.
!-------------------------------------------------------------------------------
    subroutine uc_to_u(N,uc_in,u_out,t_out,mu_out)
    use var_kind_def
    use var_air, only: rr,gk1,t_Suth
    use var_turb, only: is_tur_cal
    implicit none
    integer(dpI),intent(in):: N
    real   (dpR),intent(in):: uc_in(5,*)
    integer(dpI):: i
    real   (dpR):: u_out(5,*),t_out(*),mu_out(*),ke,mu

    if(N .le. 0)    return
    do i=1,N
        u_out(1  ,i)=  uc_in(1  ,i)
        u_out(2:4,i)=  uc_in(2:4,i)/uc_in(1,i)
        ke          =  0.5d0*(u_out(2,i)**2+u_out(3,i)**2+u_out(4,i)**2)
        u_out(5  ,i)=  abs(uc_in(5,i)-uc_in(1,i)*ke)*gk1
        t_out(    i)=  u_out(5,i)/(rr*uc_in(1,i))
        mu  =  1.461d-6*sqrt(t_out(i)**3)/(t_out(i)+t_suth)
        if(is_tur_cal) then
            mu_out(2*i-1)   =  mu
        else
            mu_out(i)       =  mu
        end if
    end do

    return
    end subroutine uc_to_u
!-------------------------------------------------------------------------------
!   Jacobian, Flux over conservative variable.
!-------------------------------------------------------------------------------
    subroutine jac_Fn_U(n,u,vg,J)
    use var_kind_def
    use var_air, only: gk,gk1
    implicit none
    real(dpR):: n(*),u(*),vg(*),J(5,*),uan,vgn,un,ke,H

    ke  =  0.5d0*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
    H   =  gk*u(5)/(gk1*u(1))+ke
    uan =  u(2)*n(1)+u(3)*n(2)+u(4)*n(3)
    vgn =  vg(1)*n(1)+vg(2)*n(2)+vg(3)*n(3)
    un  =  uan-vgn
    J(1,1)  = -vgn
    J(2,1)  = -uan*u(2)+gk1*ke*n(1)
    J(3,1)  = -uan*u(3)+gk1*ke*n(2)
    J(4,1)  = -uan*u(4)+gk1*ke*n(3)
    J(5,1)  = (gk1*ke-H)*uan
    J(1,2)  =  n(1)
    J(2,2)  = (2.0d0-gk)*u(2)*n(1)+un
    J(3,2)  =  u(3)*n(1)-gk1*u(2)*n(2)
    J(4,2)  =  u(4)*n(1)-gk1*u(2)*n(3)
    J(5,2)  =  H   *n(1)-gk1*u(2)*uan
    J(1,3)  =  n(2)
    J(2,3)  =  u(2)*n(2)-gk1*u(3)*n(1)
    J(3,3)  = (2.0d0-gk)*u(3)*n(2)+un
    J(4,3)  =  u(4)*n(2)-gk1*u(3)*n(3)
    J(4,3)  =  H   *n(2)-gk1*u(3)*uan
    J(1,4)  =  n(3)
    J(2,4)  =  u(2)*n(3)-gk1*u(4)*n(1)
    J(3,4)  =  u(3)*n(3)-gk1*u(4)*n(2)
    J(4,4)  = (2.0d0-gk)*u(4)*n(3)+un
    J(5,4)  =  H*n(3)-gk1*u(4)*uan
    J(1,5)  =  0.0d0
    J(2:4,5)=  gk1*n(1:3)
    J(5,5)  =  un+gk1*uan
    J(1:5,1:5)  =  J(1:5,1:5)*n(4)

    return
    end subroutine jac_Fn_U
!-------------------------------------------------------------------------------
!   $\fpp{F_n}{W}$.
!-------------------------------------------------------------------------------
    subroutine Fn_over_W(n,u,vg,J,sp)
    use var_kind_def
    use var_air, only: gk,gk1
    implicit none
    real(dpR):: n(*),u(*),vg(*),J(5,*),sp,un,vn,H,p_h_u(5)

    un      =  u(2)*n(1)+u(3)*n(2)+u(4)*n(3)
    vn      =  vg(1)*n(1)+vg(2)*n(2)+vg(3)*n(3)
    H       =  gk*u(5)/(u(1)*gk1)+0.5d0*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
    p_h_u   = (/-gk*u(5)/(gk1*u(1)**2), u(2:4), gk/(gk1*u(1))/)
    J(1,1:5)= (/un-vn, u(1)*n(1:3), 0.0d0/)
    J(2,1:5)=  u(2)*J(1,1:5)+(/0.0d0, u(1)*(un-vn), 0.0d0, 0.0d0, n(1)/)
    J(3,1:5)=  u(3)*J(1,1:5)+(/0.0d0, 0.0d0, u(1)*(un-vn), 0.0d0, n(2)/)
    J(4,1:5)=  u(4)*J(1,1:5)+(/0.0d0, 0.0d0, 0.0d0, u(1)*(un-vn), n(3)/)
    J(5,1:5)=  H   *J(1,1:5)+u(1)*(un-vn)*p_h_u(1:5)
    J(1:5,1:5)  =  J(1:5,1:5)*n(4)
    sp      = (abs(un-vn)+sqrt(gk*u(5)/u(1)))*n(4)

    return
    end subroutine Fn_over_W
!-------------------------------------------------------------------------------
!   Jacobian, Fn over n.
!-------------------------------------------------------------------------------
    subroutine jac_Fn_n(vg,u,J)
    use var_kind_def
    use var_air, only: gk,gk1
    implicit none
    real   (dpR),intent(in):: vg(*),u(*)
    real   (dpR):: J(5,*),H

    H       =  gk*u(5)/(gk1*u(1))+0.5d0*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
    J(1,1:3)=  u(1)*(u(2:4)-vg(1:3))
    J(2,1:3)=  J(1,1:3)*u(2)+(/u(5) , 0.0d0, 0.0d0/)
    J(3,1:3)=  J(1,1:3)*u(3)+(/0.0d0, u(5) , 0.0d0/)
    J(4,1:3)=  J(1,1:3)*u(4)+(/0.0d0, 0.0d0, u(5) /)
    J(5,1:3)=  J(1,1:3)*H   +u(5)*vg(1:3)

    return
    end subroutine jac_Fn_n
!-------------------------------------------------------------------------------
!   calculate Jacobians for the SW solver.
!-------------------------------------------------------------------------------
    subroutine jac_sw(is_trans,n,uL,uR,JL,JR)
    use var_kind_def
    implicit none
    logical(dpR),intent(in):: is_trans
    real   (dpR):: n(*),uL(*),uR(*),JL(5,*),JR(5,*),F_U(5,5),U_u(5,5),spra,akh(3), &
                &  eig(3),rtmm(5,5)

    call u_to_akh(1, uL, akh)
    call reigl(uL,(/0.0d0, 0.0d0, 0.0d0/),akh,n, 1.0d0,1.0d0,0.0d0,F_U,spra,eig)
    call U_over_W(uL, U_u, rtmm)
    JL(1:5,1:5) =  matmul(F_U, U_u)

    call u_to_akh(1, uR, akh)
    call reigl(uR,(/0.0d0, 0.0d0, 0.0d0/),akh,n,-1.0d0,1.0d0,0.0d0,F_U,spra,eig)
    call U_over_W(uR, U_u, rtmm)
    JR(1:5,1:5) =  matmul(F_U, U_u)

    if(is_trans) then
        call matrix_transpose(5, JL)
        call matrix_transpose(5, JR)
    end if

    return
    end subroutine jac_sw
!-------------------------------------------------------------------------------
!   transformation between conservative and primitive variables.
!-------------------------------------------------------------------------------
    subroutine W_over_U(u,WU)
    use var_kind_def
    use var_air
    implicit none
    real(dpR):: u(*),WU(*)

    WU(1:25)=  0.0d0
    WU(1   )=  1.0d0
    WU(2:4 )= -u(2:4)/u(1)
    WU(5   )=  0.5d0*gk1*(u(2)**2+u(3)**2+u(4)**2)
    WU(7   )=  1.0d0/u(1)
    WU(10  )= -gk1*u(2)
    WU(13  )=  1.0d0/u(1)
    WU(15  )= -gk1*u(3)
    WU(19  )=  1.0d0/u(1)
    WU(20  )= -gk1*u(4)
    WU(25  )=  gk1

    return
    end subroutine W_over_U
!-------------------------------------------------------------------------------
!   transformation between conservative and primitive variables.
!-------------------------------------------------------------------------------
    subroutine U_over_W(u,UW,WU)
    use var_kind_def
    use var_air
    implicit none
    real(dpR):: u(*),UW(*),WU(*)

    UW(1:25)=  0.0d0
    UW(1   )=  1.0d0
    UW(2:4 )=  u(2:4)
    UW(5   )=  0.5d0*(u(2)**2+u(3)**2+u(4)**2)
    UW(7   )=  u(1)
    UW(10  )=  u(1)*u(2)
    UW(13  )=  u(1)
    UW(15  )=  u(1)*u(3)
    UW(19  )=  u(1)
    UW(20  )=  u(1)*u(4)
    UW(25  )=  1.0d0/gk1

    WU(1:25)=  0.0d0
    WU(1   )=  1.0d0
    WU(2:4 )= -u(2:4)/u(1)
    WU(5   )=  0.5d0*gk1*(u(2)**2+u(3)**2+u(4)**2)
    WU(7   )=  1.0d0/u(1)
    WU(10  )= -gk1*u(2)
    WU(13  )=  1.0d0/u(1)
    WU(15  )= -gk1*u(3)
    WU(19  )=  1.0d0/u(1)
    WU(20  )= -gk1*u(4)
    WU(25  )=  gk1

    return
    end subroutine U_over_W
!-------------------------------------------------------------------------------
!   isentropic vortex.
!-------------------------------------------------------------------------------
    subroutine iv2d(c,r,u)
    use var_kind_def
    use var_air, only: gk,gk1,rr
    use var_global, only: pi
    use var_bndv, only: u_fs,t_fs,a_fs
    implicit none
    real(dpR):: c(*),r(*),u(*),cv,rv,d(2),r2,alp,phi,t,e

!   u(1)    =  u_fs(1)*(1.0d0+5.0d-2*sin(0.2d0*pi*r(1)))
!   u(2:4)  =  u_fs(2:4)
!   u(5  )  =  u_fs(5)
!   return

!   the Lodato implementation.
    cv      =  5.0d0
    rv      =  5.0d0
    d(1:2)  =  r(1:2)-c(1:2)

    cv      =  5.0d0
    rv      =  1.0d0
    d(1:2)  =  r(1:2)-c(1:2)
    r2      =  d(1)**2+d(2)**2
    u(2)    = -cv*d(2)*dexp(-r2*0.5d0/(rv*rv))/(rv*rv)+u_fs(2)
    u(3)    =  cv*d(1)*dexp(-r2*0.5d0/(rv*rv))/(rv*rv)
    u(4)    =  0.0d0
    u(5)    =  u_fs(5)*dexp(-0.5d0*gk*(cv/(a_fs*rv))**2*dexp(-r2/(rv*rv)))
    u(1)    =  u(5)/(rr*t_fs)
    return

!   the Mavriplis implementation.
    alp =  4.0d0
    phi =  1.0d0
    e   =  dexp(phi*(1.0d0-r2))
    u(2)= -0.5d0*alp*d(2)*e/pi*a_fs+u_fs(2)
    u(3)=  0.5d0*alp*d(1)*e/pi*a_fs+u_fs(3)
    u(4)=  0.0d0
    t   = (1.0d0-alp*alp*gk1/(1.6d1*phi*gk*pi*pi)*e*e)*t_fs
    u(5)=  u_fs(5)*(t_fs/t)**(-gk/gk1)
    u(1)=  u(5)/(rr*t)

    return
    end subroutine iv2d
!-------------------------------------------------------------------------------
!   Taylor-Green vortex.
!-------------------------------------------------------------------------------
    subroutine taylor_green_vortex(xyz,u)
    use var_kind_def
    use var_global, only: pi
    use var_bndv, only: u_fs,ua_fs
    implicit none
    real(dpR):: xyz(*),u(*),x,y

    u(1)=  1.0d0
    x   = (xyz(1)+5.0d0)*0.1d0
    y   = (xyz(2)+5.0d0)*0.1d0
    u(2)=  ua_fs*sin(2.0d0*pi*x)*cos(2.0d0*pi*y)
    u(3)= -ua_fs*cos(2.0d0*pi*x)*sin(2.0d0*pi*y)
    u(4)=  0.0d0
    u(5)=  u_fs(5)+0.25d0*ua_fs*ua_fs*(cos(4.0d0*pi*x)+cos(4.0d0*pi*y))

    return
    end subroutine taylor_green_vortex
!-------------------------------------------------------------------------------
!   compute flux.
!-------------------------------------------------------------------------------
    subroutine u_to_F(M,n,u,vg,f)
    use var_kind_def
    use var_air
    implicit none
    integer(dpI),intent(in):: M
    real   (dpR),intent(in):: n(3,*),u(5,*),vg(3,*)
    integer(dpI):: i
    real   (dpR):: f(5,*),H,un,vn

    do i=1,M
        un      =  u (2,i)*n(1,i)+u (3,i)*n(2,i)+u (4,i)*n(3,i)
        vn      =  vg(1,i)*n(1,i)+vg(2,i)*n(2,i)+vg(3,i)*n(3,i)
        f(1  ,i)=  u(1,i)*(un-vn)
        f(2:4,i)=  f(1,i)*u(2:4,i)+u(5,i)*n(1:3,i)
        H       =  gk*u(5,i)/(gk1*u(1,i))+0.5d0*(u(2,i)**2+u(3,i)**2+u(4,i)**2)
        f(5  ,i)=  f(1,i)*H+u(5,i)*vn
    end do

    return
    end subroutine u_to_F
!-------------------------------------------------------------------------------
!   compute viscous flux.
!-------------------------------------------------------------------------------
    subroutine get_vis_fac(M,e,u,mu,g,k,rhsD)
    use var_kind_def
    use var_air, only: cp,prL,prT
    use var_turb, only: is_tur_cal,is_Spalart_QCR,is_KO
    implicit none
    integer(dpI),intent(in):: M
    real   (dpR),intent(in):: e(4,*),u(5,*),mu(2,*),g(12,*),k(*)
    integer(dpI):: i
    real   (dpR):: rhsD(5,*),mul,mut,kq,dU3,S(3,3),Tlam(3,3),Ttur(3,3),f(2:5),gra(12),turk

    do i=1,M
        mul     =  mu(1,i)
        gra(1:12)   =  g(1:12,i)
        dU3     = (gra(1)+gra(5)+gra(9))/3.0d0
        S(1,1)  =  gra(1)-dU3
        S(2,1)  =  0.5d0*(gra(2)+gra(4))
        S(3,1)  =  0.5d0*(gra(3)+gra(7))
        S(1,2)  =  S(2,1)
        S(2,2)  =  gra(5)-dU3
        S(3,2)  =  0.5d0*(gra(6)+gra(8))
        S(1,3)  =  S(3,1)
        S(2,3)  =  S(3,2)
        S(3,3)  =  gra(9)-dU3
        Tlam    =  2.0d0*mul*S
        kq      =  cp*mul/prL

        if(is_tur_cal) then
            mut =  mu(2,i)
            Ttur=  2.0d0*mut*S
            if(is_KO) then
                turk        = -2.0d0/3.0d0*u(1,i)*k(i)
                Ttur(1,1)   =  Ttur(1,1)+turk
                Ttur(2,2)   =  Ttur(2,2)+turk
                Ttur(3,3)   =  Ttur(3,3)+turk
            end if
            if(is_Spalart_QCR)  call get_QCR_stress(mut, gra, S, Ttur)
            kq  =  kq+cp*mut/prT
            Tlam=  Tlam+Ttur
        end if

        f(2)=  e(1,i)*Tlam(1,1)+e(2,i)*Tlam(1,2)+e(3,i)*Tlam(1,3)
        f(3)=  e(1,i)*Tlam(2,1)+e(2,i)*Tlam(2,2)+e(3,i)*Tlam(2,3)
        f(4)=  e(1,i)*Tlam(3,1)+e(2,i)*Tlam(3,2)+e(3,i)*Tlam(3,3)
        f(5)=  u(2,i)*f(2)+u(3,i)*f(3)+u(4,i)*f(4) &
            & +kq*(e(1,i)*gra(10)+e(2,i)*gra(11)+e(3,i)*gra(12))
        rhsD(1  ,i) =  0.0d0
        rhsD(2:5,i) =  f(2:5)*e(4,i)
    end do

    return
    end subroutine get_vis_fac
!-------------------------------------------------------------------------------
!   compute viscous flux.
!-------------------------------------------------------------------------------
    subroutine get_vis_fac_test(M,e,u,mu,g,k,rhsD)
    use var_kind_def
    use var_air, only: cp,prL,prT
    use var_turb, only: is_tur_cal,is_Spalart_QCR,is_KO
    implicit none
    integer(dpI),intent(in):: M
    real   (dpR),intent(in):: e(4,*),u(5,*),mu(2,*),g(12,*),k(*)
    integer(dpI):: i
    real   (dpR):: rhsD(5,*),mul,mut,kq,dU3,S(3,3),Tlam(3,3),Ttur(3,3),f(2:5),gra(12),turk

    do i=1,M
        mul     =  mu(1,i)
        gra(1:12)   =  g(1:12,i)
        dU3     = (gra(1)+gra(5)+gra(9))/3.0d0
        S(1,1)  =  gra(1)-dU3
        S(2,1)  =  0.5d0*(gra(2)+gra(4))
        S(3,1)  =  0.5d0*(gra(3)+gra(7))
        S(1,2)  =  S(2,1)
        S(2,2)  =  gra(5)-dU3
        S(3,2)  =  0.5d0*(gra(6)+gra(8))
        S(1,3)  =  S(3,1)
        S(2,3)  =  S(3,2)
        S(3,3)  =  gra(9)-dU3
        Tlam    =  2.0d0*mul*S
        kq      =  cp*mul/prL

        if(is_tur_cal) then
            mut =  mu(2,i)
            Ttur=  2.0d0*mut*S
            ! if(is_KO) then
            !     turk        = -2.0d0/3.0d0*u(1,i)*k(i)
            !     Ttur(1,1)   =  Ttur(1,1)+turk
            !     Ttur(2,2)   =  Ttur(2,2)+turk
            !     Ttur(3,3)   =  Ttur(3,3)+turk
            ! end if
            if(is_Spalart_QCR)  call get_QCR_stress(mut, gra, S, Ttur)
            kq  =  kq+cp*mut/prT
            Tlam=  Tlam+Ttur
        end if

        f(2)=  e(1,i)*Tlam(1,1)+e(2,i)*Tlam(1,2)+e(3,i)*Tlam(1,3)
        f(3)=  e(1,i)*Tlam(2,1)+e(2,i)*Tlam(2,2)+e(3,i)*Tlam(2,3)
        f(4)=  e(1,i)*Tlam(3,1)+e(2,i)*Tlam(3,2)+e(3,i)*Tlam(3,3)
        f(5)=  u(2,i)*f(2)+u(3,i)*f(3)+u(4,i)*f(4) &
            & +kq*(e(1,i)*gra(10)+e(2,i)*gra(11)+e(3,i)*gra(12))

        rhsD(1  ,i) =  0.0d0
        rhsD(2:5,i) =  f(2:5)*e(4,i)
    end do

    return
    end subroutine get_vis_fac_test
!-------------------------------------------------------------------------------
!   get the QCR Reynolds Stress.
!-------------------------------------------------------------------------------
    subroutine get_QCR_stress(mut,gra,S,Ttur)
    use var_kind_def
    use var_global, only: uref,L_ref
    use var_turb, only: RANS_model,RANS_SA,QCR_version
    implicit none
    real(dpR),intent(in):: mut,gra(*),S(3,*)
    real(dpR):: Ttur(3,*),O(3,3),T(3,3),D(3,3),rtmp

    rtmp=  gra(1)*gra(1)+gra(2)*gra(2)+gra(3)*gra(3) &
        & +gra(4)*gra(4)+gra(5)*gra(5)+gra(6)*gra(6) &
        & +gra(7)*gra(7)+gra(8)*gra(8)+gra(9)*gra(9)
    rtmp=  max(sqrt(rtmp), 1.0d-6*uref/L_ref)
    O(1,1)  =  0.0d0
    O(2,1)  =  gra(4)-gra(2)
    O(3,1)  =  gra(7)-gra(3)
    O(1,2)  = -O(2,1)
    O(2,2)  =  0.0d0
    O(3,2)  =  gra(8)-gra(6)
    O(1,3)  = -O(3,1)
    O(2,3)  = -O(3,2)
    O(3,3)  =  0.0d0
    O       =  O/rtmp

    T(1:3,1:3)  =  Ttur(1:3,1:3)
    D(1,1)  = (O(1,2)*T(1,2)+O(1,3)*T(1,3))*2.0d0
    D(2,1)  =  O(2,1)*T(1,1)+O(2,2)*T(1,2)+O(2,3)*T(1,3) &
            & +O(1,1)*T(2,1)+O(1,2)*T(2,2)+O(1,3)*T(2,3)
    D(3,1)  =  O(3,1)*T(1,1)+O(3,2)*T(1,2)+O(3,3)*T(1,3) &
            & +O(1,1)*T(3,1)+O(1,2)*T(3,2)+O(1,3)*T(3,3)
    D(1,2)  =  D(2,1)
    D(2,2)  = (O(2,1)*T(2,1)+O(2,3)*T(2,3))*2.0d0
    D(3,2)  =  O(3,1)*T(2,1)+O(3,2)*T(2,2)+O(3,3)*T(2,3) &
            & +O(2,1)*T(3,1)+O(2,2)*T(3,2)+O(2,3)*T(3,3)
    D(1,3)  =  D(3,1)
    D(2,3)  =  D(3,2)
    D(3,3)  = (O(3,1)*T(3,1)+O(3,2)*T(3,2))*2.0d0

    Ttur(1:3,1:3)   =  T-0.3d0*D

    if(QCR_version .eq. 2013) then
        if(RANS_model .eq. RANS_SA) then
            rtmp=  S(1,1)**2+S(2,1)**2+S(3,1)**2+S(1,2)**2+S(2,2)**2+S(3,2)**2 &
                & +S(1,3)**2+S(2,3)**2+S(3,3)**2
            rtmp=  2.5d0*mut*dsqrt(2.0d0*rtmp)
            Ttur(1,1)   =  Ttur(1,1)-rtmp
            Ttur(2,2)   =  Ttur(2,2)-rtmp
            Ttur(3,3)   =  Ttur(3,3)-rtmp
        end if
    end if

    return
    end subroutine get_QCR_stress
!-------------------------------------------------------------------------------
!   A*d, the invisicd part.
!-------------------------------------------------------------------------------
    subroutine AprdV(N,e,vg,u,akh,sgn,efix,d,v)
    use var_kind_def
    use var_air, only: gk1
    implicit none
    integer(dpI),intent(in):: N
    real   (dpR),intent(in):: e(5,*),vg(3,*),u(5,*),akh(3,*),sgn(*),efix(*),d(5,*)
    integer(dpI):: i
    real   (dpR):: v(5,*),uan,un,eig(3),eiga,eigb,tmp1,tmp2,c1,c2

    do i=1,N
        uan     =  u(2,i)*e(1,i)+u(3,i)*e(2,i)+u(4,i)*e(3,i)
        un      =  uan-(vg(1,i)*e(1,i)+vg(2,i)*e(2,i)+vg(3,i)*e(3,i))
        eig(1)  =  un
        eig(2)  =  eig(1)+akh(1,i)
        eig(3)  =  eig(1)-akh(1,i)
        eig     =  0.5d0*e(4,i)*(eig+sgn(i)*(abs(eig)+efix(i)*(abs(un)+akh(1,i))))

        eiga=  0.5d0*(eig(2)+eig(3))-eig(1)
        eigb=  0.5d0*(eig(2)-eig(3))
        tmp1= -uan     *d(1,i)+e(1,i)*d(2,i)+e(2,i)*d(3,i)+e(3,i)*d(4,i)
        tmp2= (akh(2,i)*d(1,i)-u(2,i)*d(2,i)-u(3,i)*d(3,i)-u(4,i)*d(4,i)+d(5,i))*gk1

        un      =  1.0d0/akh(1,i)
        c1      = (eiga*tmp2*un+eigb*tmp1)*un
        c2      =  eigb*tmp2*un+eiga*tmp1 
        v(1:5,i)=  eig(1)*d(1:5,i)
        v(1  ,i)=  v(1  ,i)+c1
        v(2:4,i)=  v(2:4,i)+c1*u  (2:4,i)+c2*e(1:3,i)
        v(5  ,i)=  v(5  ,i)+c1*akh(3  ,i)+c2*uan
!       v(1:5,i)=  eig(1)*d(1:5,i) &
!               &+(eiga*tmp2*un+eigb*tmp1)*un*(/1.0d0, u(2:4,i), akh(3,i)/) &
!               &+(eigb*tmp2*un+eiga*tmp1)*   (/0.0d0, e(1:3,i), uan     /)
    end do

    return
    end subroutine AprdV
!-------------------------------------------------------------------------------
!   A*d, the viscous part.
!-------------------------------------------------------------------------------
    subroutine AprdV_vis(N,e,uL,uR,mu_L,mu_R,d,v)
    use var_kind_def
    use var_air, only: gk1,rr,cp,prL,prT
    use var_global, only: R13
    use var_turb, only: is_tur_cal
    implicit none
    integer(dpI),intent(in):: N
    real   (dpR),intent(in):: uL(*),uR(5,*),mu_L(*),mu_R(2,*),e(5,*),d(5,*)
    integer(dpI):: i
    real   (dpR):: v(5,*),ke,un,mu,kq,muL,muT,s1,s2,s3,u(5),rh1,rtmp

    muT =  0.0d0
    do i=1,N
        u(1:5)  =  0.5d0*(uL(1:5)+uR(1:5,i))
        rh1     =  1.0d0/u(1)
        ke      =  0.5d0*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
        un      =  u(2)*e(1,i)+u(3)*e(2,i)+u(4)*e(3,i)
        mul     =  0.5d0*(mu_L(1)+mu_R(1,i))
        if(is_tur_cal)  muT =  0.5d0*(mu_L(2)+mu_R(2,i))

        mu  =  muL+muT
        kq  =  cp*(muL/prL+muT/prT)
        rtmp=  kq/(mu*rr)
        s1  =  R13*(-un*d(1,i)+e(1,i)*d(2,i)+e(2,i)*d(3,i)+e(3,i)*d(4,i))
        s2  =  u(2)*d(2,i)+u(3)*d(3,i)+u(4)*d(4,i)
        s3  =  gk1*rtmp*(ke*d(1,i)-s2+d(5,i))-u(5)*d(1,i)*rtmp*rh1
        v(2:4,i)=  s1*e(1:3,i)-u(2:4)*d(1,i)+d(2:4,i)
        v(5  ,i)=  s1*un-ke*d(1,i)*2.0d0+s2+s3
        v(2:5,i)=  v(2:5,i)*mu*e(4,i)*e(5,i)*rh1
    end do

    return
    end subroutine AprdV_vis
!-------------------------------------------------------------------------------
!   R*D*L.
!-------------------------------------------------------------------------------
    subroutine reigl(u,vg,akh,e,sgn_eig,sgn_mat,eps,RDL,spra,eig)
    use var_kind_def
    use var_air
    implicit none
    integer(dpI):: j
    real   (dpR):: u(5),vg(3),akh(*),e(*),sgn_eig,sgn_mat,eps,spra,eig(*),rdl(5,*), &
                &  an(3),qn,vn,qt1,qt2,qt3,phi,R(5,5),L(5,5),tv(5)

    an(1:3) =  akh(1)*e(1:3)
    qn      =  e(1)*u (2)+e(2)*u (3)+e(3)*u (4)
    vn      =  e(1)*vg(1)+e(2)*vg(2)+e(3)*vg(3)
    qt1     =  u(3)*e(3)-u(4)*e(2)
    qt2     =  u(4)*e(1)-u(2)*e(3)
    qt3     =  u(2)*e(2)-u(3)*e(1)
    phi     =  akh(1)*akh(1)-gk1*akh(2)

    tv(1)   =  qn-vn
    tv(2)   =  tv(1)+akh(1)
    tv(3)   =  tv(1)-akh(1)
    spra    =  abs(tv(1))+akh(1)
    eig(1:3)=  0.5d0*sgn_mat*e(4)*(tv(1:3)+sgn_eig*(abs(tv(1:3))+eps*spra))

    tv(1)   =       e(1)
    tv(2)   =  u(2)*e(1)
    tv(3)   =  u(3)*e(1)+an(3)
    tv(4)   =  u(4)*e(1)-an(2)
    tv(5)   =  akh(2)*e(1)+akh(1)*qt1
    R(:,1)  =  tv(:)*eig(1)
    tv(1)   =       e(2)
    tv(2)   =  u(2)*e(2)-an(3)
    tv(3)   =  u(3)*e(2)
    tv(4)   =  u(4)*e(2)+an(1)
    tv(5)   =  akh(2)*e(2)+akh(1)*qt2
    R(:,2)  =  tv(:)*eig(1)
    tv(1)   =       e(3)
    tv(2)   =  u(2)*e(3)+an(2)
    tv(3)   =  u(3)*e(3)-an(1)
    tv(4)   =  u(4)*e(3)
    tv(5)   =  akh(2)*e(3)+akh(1)*qt3
    R(:,3)  =  tv(:)*eig(1)
    tv(1)   =  1.0d0
    tv(2:4) =  u(2:4)+an(1:3)
    tv(5)   =  akh(3)+akh(1)*qn
    R(:,4)  =  tv(:)*eig(2)
    tv(1)   =  1.0d0
    tv(2:4) =  u(2:4)-an(1:3)
    tv(5)   =  akh(3)-akh(1)*qn
    R(:,5)  =  tv(:)*eig(3)

    L(1,1)  =  phi*e(1)-akh(1)*qt1
    L(2,1)  =  phi*e(2)-akh(1)*qt2
    L(3,1)  =  phi*e(3)-akh(1)*qt3
    L(4,1)  =  0.5d0*(gk1*akh(2)-akh(1)*qn)
    L(5,1)  =  0.5d0*(gk1*akh(2)+akh(1)*qn)
    L(1,2)  =  gk1*u(2)*e(1)
    L(2,2)  =  gk1*u(2)*e(2)-an(3)
    L(3,2)  =  gk1*u(2)*e(3)+an(2)
    L(4,2)  = -0.5d0*(gk1*u(2)-an(1))
    L(5,2)  = -0.5d0*(gk1*u(2)+an(1))
    L(1,3)  =  gk1*u(3)*e(1)+an(3)
    L(2,3)  =  gk1*u(3)*e(2)
    L(3,3)  =  gk1*u(3)*e(3)-an(1)
    L(4,3)  = -0.5d0*(gk1*u(3)-an(2))
    L(5,3)  = -0.5d0*(gk1*u(3)+an(2))
    L(1,4)  =  gk1*u(4)*e(1)-an(2)
    L(2,4)  =  gk1*u(4)*e(2)+an(1)
    L(3,4)  =  gk1*u(4)*e(3)
    L(4,4)  = -0.5d0*(gk1*u(4)-an(3))
    L(5,4)  = -0.5d0*(gk1*u(4)+an(3))
    L(1:3,5)= -gk1*e(1:3)
    L(4:5,5)=  0.5d0*gk1
    L       =  L/(akh(1)*akh(1))

    do j=1,5
        RDL(1:5,j)  =  R(1:5,1)*L(1,j)+R(1:5,2)*L(2,j)+R(1:5,3)*L(3,j) &
                    & +R(1:5,4)*L(4,j)+R(1:5,5)*L(5,j)
    end do
    spra=  spra*e(4)

    return
    end subroutine reigl
!-------------------------------------------------------------------------------
!   R*D*L.
!-------------------------------------------------------------------------------
    subroutine reigl_new(u,vg,e,sgn_eig,sgn_mat,eps,RDL,spra,eig)
    use var_kind_def
    use var_air
    implicit none
    real   (dpR),intent(in):: u(*),vg(*),e(*),sgn_eig,sgn_mat,eps
    real   (dpR):: a,spra,eig(*),rdl(5,*),qn,vn,eiga,eigb,m(5),n(5),ke,H,a1,a2

    a       =  sqrt(gk*u(5)/u(1))
    qn      =  e(1)*u (2)+e(2)*u (3)+e(3)*u (4)
    vn      =  e(1)*vg(1)+e(2)*vg(2)+e(3)*vg(3)
    m(1)    =  qn-vn
    m(2)    =  m(1)+a
    m(3)    =  m(1)-a
    spra    =  abs(m(1))+a
    eig(1:3)=  0.5d0*sgn_mat*e(4)*(m(1:3)+sgn_eig*(abs(m(1:3))+eps*spra))
    eiga    =  0.5d0*(eig(2)-eig(3))
    eigb    =  0.5d0*(eig(2)+eig(3))-eig(1)

    m(1)    = -qn
    m(2:4)  =  e(1:3)
    m(5)    =  0.0d0
    ke      =  0.5d0*(u(2)**2+u(3)**2+u(4)**2)
    n(1)    =  gk1*ke
    n(2:4)  = -gk1*u(2:4)
    n(5)    =  gk1
    H       =  gk*u(5)/(gk1*u(1))+ke

    a1          =  1.0d0/a
    a2          =  a1*a1
    RDL(1,1:5)  =  eiga*a1*      m(1:5)             +eigb                  *a2*n(1:5)
    RDL(2,1:5)  =  eiga*a1*(u(2)*m(1:5)+e(1)*n(1:5))+eigb*(e(1)*m(1:5)+u(2)*a2*n(1:5))
    RDL(3,1:5)  =  eiga*a1*(u(3)*m(1:5)+e(2)*n(1:5))+eigb*(e(2)*m(1:5)+u(3)*a2*n(1:5))
    RDL(4,1:5)  =  eiga*a1*(u(4)*m(1:5)+e(3)*n(1:5))+eigb*(e(3)*m(1:5)+u(4)*a2*n(1:5))
    RDL(5,1:5)  =  eiga*a1*(H   *m(1:5)+qn  *n(1:5))+eigb*(qn  *m(1:5)+H   *a2*n(1:5))
    RDL(1,1)    =  RDL(1,1)+eig(1)
    RDL(2,2)    =  RDL(2,2)+eig(1)
    RDL(3,3)    =  RDL(3,3)+eig(1)
    RDL(4,4)    =  RDL(4,4)+eig(1)
    RDL(5,5)    =  RDL(5,5)+eig(1)
!   do i=1,5
!       RDL(i,i)=  RDL(i,i)+eig(1)
!   end do
    spra=  spra*e(4)

    return
    end subroutine reigl_new
!-------------------------------------------------------------------------------
!   R*D*L.
!-------------------------------------------------------------------------------
    subroutine reigl_new_test(u,vg,e,sgn_eig,sgn_mat,eps,RDL,spra,eig)
    use var_kind_def
    use var_air
    implicit none
    real   (dpR),intent(in):: u(*),vg(*),e(*),sgn_eig,sgn_mat,eps
    real   (dpR):: a,spra,eig(*),rdl(5,*),qn,vn,eiga,eigb,m(5),n(5),ke,H,a1,a2

    a       =  sqrt(gk*u(5)/u(1))
    qn      =  e(1)*u (2)+e(2)*u (3)+e(3)*u (4)
    vn      =  e(1)*vg(1)+e(2)*vg(2)+e(3)*vg(3)
    m(1)    =  qn-vn
    m(2)    =  m(1)+a
    m(3)    =  m(1)-a
    spra    =  abs(m(1))+a
    eig(1:3)=  0.5d0*sgn_mat*e(4)*(m(1:3)+sgn_eig*(abs(m(1:3))+eps*spra))
    eiga    =  0.5d0*(eig(2)-eig(3))
    eigb    =  0.5d0*(eig(2)+eig(3))-eig(1)

    m(1)    = -qn
    m(2:4)  =  e(1:3)
    m(5)    =  0.0d0
    ke      =  0.5d0*(u(2)**2+u(3)**2+u(4)**2)
    n(1)    =  gk1*ke
    n(2:4)  = -gk1*u(2:4)
    n(5)    =  gk1
    H       =  gk*u(5)/(gk1*u(1))+ke

    a1          =  1.0d0/a
    a2          =  a1*a1
    RDL(1,1:5)  =  eiga*a1*      m(1:5)             +eigb                  *a2*n(1:5)
    RDL(2,1:5)  =  eiga*a1*(u(2)*m(1:5)+e(1)*n(1:5))+eigb*(e(1)*m(1:5)+u(2)*a2*n(1:5))
    RDL(3,1:5)  =  eiga*a1*(u(3)*m(1:5)+e(2)*n(1:5))+eigb*(e(2)*m(1:5)+u(3)*a2*n(1:5))
    RDL(4,1:5)  =  eiga*a1*(u(4)*m(1:5)+e(3)*n(1:5))+eigb*(e(3)*m(1:5)+u(4)*a2*n(1:5))
    RDL(5,1:5)  =  eiga*a1*(H   *m(1:5)+qn  *n(1:5))+eigb*(qn  *m(1:5)+H   *a2*n(1:5))
    RDL(1,1)    =  RDL(1,1)+eig(1)
    RDL(2,2)    =  RDL(2,2)+eig(1)
    RDL(3,3)    =  RDL(3,3)+eig(1)
    RDL(4,4)    =  RDL(4,4)+eig(1)
    RDL(5,5)    =  RDL(5,5)+eig(1)
!   do i=1,5
!       RDL(i,i)=  RDL(i,i)+eig(1)
!   end do
    spra=  spra*e(4)

    return
    end subroutine reigl_new_test
!-------------------------------------------------------------------------------
!   R*D*L for viscous flux.
!-------------------------------------------------------------------------------
    subroutine reigl_vis(u,mul,mut,e,tmpa,spra)
    use var_kind_def
    use var_air
    implicit none
    real(dpR):: u(5),mul,mut,e(5),tmpa(*),spra,un,ke,mu,mupr,kq,rtmp

    mupr    =  mul/prl+mut/prt
    un      =  e(1)*u(2)+e(2)*u(3)+e(3)*u(4)
    ke      =  0.5d0*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
    mu      =  mul+mut
    kq      =  cp*mupr
    rtmp    =  kq/(mu*rr)

    tmpa(1 )=  0.0d0
    tmpa(2 )= -u(2)-un*e(1)/3.0d0
    tmpa(3 )= -u(3)-un*e(2)/3.0d0
    tmpa(4 )= -u(4)-un*e(3)/3.0d0
    tmpa(5 )= -rtmp*u(5)/u(1)-un*un/3.0d0-ke*2.0d0+gk1*rtmp*ke
    tmpa(6 )=  0.0d0
    tmpa(7 )=  e(1)*e(1)/3.0d0+1.0d0
    tmpa(8 )=  e(1)*e(2)/3.0d0
    tmpa(9 )=  e(1)*e(3)/3.0d0
    tmpa(10)=  u(2)+un*e(1)/3.0d0-gk1*rtmp*u(2)
    tmpa(11)=  0.0d0
    tmpa(12)=  e(1)*e(2)/3.0d0
    tmpa(13)=  e(2)*e(2)/3.0d0+1.0d0
    tmpa(14)=  e(2)*e(3)/3.0d0
    tmpa(15)=  u(3)+un*e(2)/3.0d0-gk1*rtmp*u(3)
    tmpa(16)=  0.0d0
    tmpa(17)=  e(1)*e(3)/3.0d0
    tmpa(18)=  e(2)*e(3)/3.0d0
    tmpa(19)=  e(3)*e(3)/3.0d0+1.0d0
    tmpa(20)=  u(4)+un*e(3)/3.0d0-gk1*rtmp*u(4)
    tmpa(21:24) =  0.0d0
    tmpa(25)=  gk1*rtmp
    tmpa(1:25)  =  tmpa(1:25)*e(4)*mu*e(5)/u(1)
    spra=  e(4)*e(5)*gk*mupr/u(1)

    return
    end subroutine reigl_vis
!-------------------------------------------------------------------------------
!   R*D*L for viscous flux.
!-------------------------------------------------------------------------------
    subroutine reigl_vis_test(u,mul,mut,e,tmpa,spra)
    use var_kind_def
    use var_air
    implicit none
    real(dpR):: u(5),mul,mut,e(5),tmpa(*),spra,un,ke,mu,mupr,kq,rtmp

    mupr    =  mul/prl+mut/prt
    un      =  e(1)*u(2)+e(2)*u(3)+e(3)*u(4)
    ke      =  0.5d0*(u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
    mu      =  mul+mut
    kq      =  cp*mupr
    rtmp    =  kq/(mu*rr)

    tmpa(1 )=  0.0d0
    tmpa(2 )= -u(2)-un*e(1)/3.0d0
    tmpa(3 )= -u(3)-un*e(2)/3.0d0
    tmpa(4 )= -u(4)-un*e(3)/3.0d0
    tmpa(5 )= -rtmp*u(5)/u(1)-un*un/3.0d0-ke*2.0d0+gk1*rtmp*ke
    tmpa(6 )=  0.0d0
    tmpa(7 )=  e(1)*e(1)/3.0d0+1.0d0
    tmpa(8 )=  e(1)*e(2)/3.0d0
    tmpa(9 )=  e(1)*e(3)/3.0d0
    tmpa(10)=  u(2)+un*e(1)/3.0d0-gk1*rtmp*u(2)
    tmpa(11)=  0.0d0
    tmpa(12)=  e(1)*e(2)/3.0d0
    tmpa(13)=  e(2)*e(2)/3.0d0+1.0d0
    tmpa(14)=  e(2)*e(3)/3.0d0
    tmpa(15)=  u(3)+un*e(2)/3.0d0-gk1*rtmp*u(3)
    tmpa(16)=  0.0d0
    tmpa(17)=  e(1)*e(3)/3.0d0
    tmpa(18)=  e(2)*e(3)/3.0d0
    tmpa(19)=  e(3)*e(3)/3.0d0+1.0d0
    tmpa(20)=  u(4)+un*e(3)/3.0d0-gk1*rtmp*u(4)
    tmpa(21:24) =  0.0d0
    tmpa(25)=  gk1*rtmp
    tmpa(1:25)  =  tmpa(1:25)*e(4)*mu*e(5)/u(1)
    spra=  e(4)*e(5)*gk*mupr/u(1)

    spra = tmpa(5 )
    return
    end subroutine reigl_vis_test
!-------------------------------------------------------------------------------
!   $\fpp{a, ke, H}{W}$.
!-------------------------------------------------------------------------------
    subroutine akh_over_W(u,J)
    use var_kind_def
    use var_air, only: gk,gk1
    implicit none
    real(dpR),intent(in):: u(*)
    real(dpR):: J(3,*),a

    a       =  sqrt(gk*u(5)/u(1))
    J(1,1:5)=  0.5d0*(/-a/u(1), 0.0d0, 0.0d0, 0.0d0, sqrt(gk/(u(1)*u(5)))/)
    J(2,1:5)= (/0.0d0, u(2), u(3), u(4), 0.0d0/)
    J(3,1:5)= (/-gk*u(5)/(gk1*u(1)**2), u(2:4), gk/(gk1*u(1))/)

    return
    end subroutine akh_over_W
!-------------------------------------------------------------------------------
!   mu over u.
!-------------------------------------------------------------------------------
    subroutine mu_over_u(u,J)
    use var_kind_def
    use var_air
    implicit none
    real(dpR),intent(in):: u(*)
    real(dpR):: J(*),t

    t       =  u(5)/(rr*u(1))
    J(1:5)  =  1.461d-6*(1.5d0*sqrt(t)/(t+t_suth)-t**1.5d0/(t+t_suth)**2) &
            & *(/-t/u(1), 0.0d0, 0.0d0, 0.0d0, t/u(5)/)

    return
    end subroutine mu_over_u
!-------------------------------------------------------------------------------
!   t over u.
!-------------------------------------------------------------------------------
    subroutine t_over_u(u,J)
    use var_kind_def
    use var_air
    implicit none
    real(dpR),intent(in):: u(*)
    real(dpR):: J(*),t

    t       =  u(5)/(rr*u(1))
    J(1:5)  =(/-t/u(1), 0.0d0, 0.0d0, 0.0d0, t/u(5)/)

    return
    end subroutine t_over_u
!-------------------------------------------------------------------------------
!   AUSM flux formulation.
!-------------------------------------------------------------------------------
    subroutine rhs_conv_ausm_LD(calD,addD,m,sigma)
    use var_kind_def
    use var_air, only: gk,gk1
    use var_temp, only: uLf,uRf,vgf,rhsl,rhsD,fac_1d
    implicit none
    logical(dpL),intent(in):: calD,addD
    integer(dpI),intent(in):: m
    real   (dpR),intent(in):: sigma(*)
    integer(dpI):: i
    real   (dpR):: e(3),S,aL,aR,unL,unR,vn,maL,maR,hL,hR,mL,pL,mR,pR,mf,p,pm, &
                &  uL(5),uR(5),asL,asR,am,fL(5),fR(5),f(5),ua

    vn  =  0.0d0
    do i=1,M
        e(1:3)  =  fac_1d(1:3,i)
        S       =  fac_1d(4  ,i)
        uL(1:5) =  uLf(1:5,i)
        uR(1:5) =  uRf(1:5,i)
        aL      =  sqrt(gk*uL(5)/uL(1))
        aR      =  sqrt(gk*uR(5)/uR(1))
        hL      =  aL*aL/gk1+0.5d0*(uL(2)*uL(2)+uL(3)*uL(3)+uL(4)*uL(4))
        hR      =  aR*aR/gk1+0.5d0*(uR(2)*uR(2)+uR(3)*uR(3)+uR(4)*uR(4))
        pm      =  0.5d0*(uL(5)+uR(5))

        vn  =  vgf(1,i)*e(1)+vgf(2,i)*e(2)+vgf(3,i)*e(3)
        unL =  uL(2)*e(1)+uL(3)*e(2)+uL(4)*e(3)-vn
        unR =  uR(2)*e(1)+uR(3)*e(2)+uR(4)*e(3)-vn
        if(.not. calD) then
            mL      =  0.5d0*uL(1)*unL
            mR      =  0.5d0*uR(1)*unR
            f(1  )  =  mL+mR
            f(2:4)  =  mL*uL(2:4)+mR*uR(2:4)+pm*e(1:3)
            f(5  )  =  mL*hL     +mR*hR     +pm*vn
            rhsl(1:5,i) =  S*f(1:5)
            cycle
        end if

        asL     =  sqrt(2.0d0*hL*gk1/(gk+1.0d0))
        asR     =  sqrt(2.0d0*hR*gk1/(gk+1.0d0))
        fL(1:5) = (/1.0d0, uL(2:4), hL/)*uL(1)
        fR(1:5) = (/1.0d0, uR(2:4), hR/)*uR(1)
        ua      = (uL(2)-vgf(1,i))**2+(uL(3)-vgf(2,i))**2+(uL(4)-vgf(3,i))**2 &
                &+(uR(2)-vgf(1,i))**2+(uR(3)-vgf(2,i))**2+(uR(4)-vgf(3,i))**2
        ua      =  sqrt(0.5d0*ua)

        am  =  min(asL*asL/max(asL, unL), asR*asR/max(asR, -unR))
        maL =  unL/am
        maR =  unR/am
        if    (maL .ge. 1.0d0) then
            mL  =  maL
            pL  =  1.0d0
        elseif(maL .le. -1.0d0) then
            mL  =  0.0d0
            pL  =  0.0d0
        else
            mL  =  0.25d0*(MaL+1.0d0)**2+0.125d0*(MaL*MaL-1.0d0)**2
!           pL  =  0.25d0*(MaL+1.0d0)**2*(2.0d0-MaL)+0.1875d0*MaL*(MaL**2-1.0d0)**2
            pL  =  0.25d0*(MaL+1.0d0)**2*(2.0d0-MaL)
        end if

        if    (maR .le. -1.0d0) then
            mR  =  maR
            pR  =  1.0d0
        elseif(maR .ge. 1.0d0) then
            mR  =  0.0d0
            pR  =  0.0d0
        else
            mR  = -0.25d0*(MaR-1.0d0)**2-0.125d0*(MaR*MaR-1.0d0)**2
!           pR  =  0.25d0*(MaR-1.0d0)**2*(2.0d0+MaR)-0.1875d0*MaR*(MaR**2-1.0d0)**2
            pR  =  0.25d0*(MaR-1.0d0)**2*(2.0d0+MaR)
        end if
        mf  = (mL+mR)*am

        if(.true.) then
            p   =  pm+0.5d0*(pL-pR)*(uL(5)-uR(5)) &
                &+(pL+pR-1.0d0)*ua*0.25d0*(uL(1)+uR(1))*(aL+aR)*sigma(i)
        else
!           p   =  pL*uL(5)+pR*uR(5)
!           p   =  pL*uL(5)+pR*uR(5)-pL*pR*(uL(1)+uR(1))*am*(unR-unL)
            p   =  pL*uL(5)+pR*uR(5)-pL*pR*(uL(1)+uR(1))*ua*(unR-unL)
        end if
        mL      =  max(mf, 0.0d0)
        mR      =  min(mf ,0.0d0)
        f(1)    =  mL*fL(1)  +mR*fR(1)
        f(2:4)  =  mL*fL(2:4)+mR*fR(2:4)+p*e(1:3)
        f(5)    =  mL*fL(5)  +mR*fR(5)  +p*vn
        if(addD) then
            rhsl(1:5,i) =  S*f(1:5)
        else
            mL          =  0.5d0*uL(1)*unL
            mR          =  0.5d0*uR(1)*unR
            rhsl(1  ,i) =  S*(mL        +mR)
            rhsl(2:4,i) =  S*(mL*uL(2:4)+mR*uR(2:4)+pm*e(1:3))
            rhsl(5  ,i) =  S*(mL*hL     +mR*hR     +pm*vn)
            rhsD(1:5,i) =  rhsl(1:5,i)-S*f(1:5)
        end if
    end do

    return
    end subroutine rhs_conv_ausm_LD
!-------------------------------------------------------------------------------
!   get the vortex tilting parameter.
!-------------------------------------------------------------------------------
    subroutine get_vortex_tilting(g,nu,nuT,tilt)
    use var_kind_def
    implicit none
    real(dpR),intent(in):: g(*),nu,nuT
    real(dpR):: S(3,3),o(3),oo,SDo(3),v(3),t0,t1,tilt

    S(1,1)  =  g(1)
    S(2,1)  =  0.5d0*(g(2)+g(4))
    S(3,1)  =  0.5d0*(g(3)+g(7))
    S(1,2)  =  S(2,1)
    S(2,2)  =  g(5)
    S(3,2)  =  0.5d0*(g(6)+g(8))
    S(1,3)  =  S(3,1)
    S(2,3)  =  S(3,2)
    S(3,3)  =  g(9)

    o(1)    =  g(8)-g(6)
    o(2)    =  g(3)-g(7)
    o(3)    =  g(4)-g(2)
    oo      =  sqrt(o(1)**2+o(2)**2+o(3)**2)
    SDo     =  S(1:3,1)*o(1)+S(1:3,2)*o(2)+S(1:3,3)*o(3)
    v(1)    =  SDo(2)*o(3)-SDO(3)*o(2)
    v(2)    =  SDo(3)*o(1)-SDO(1)*o(3)
    v(3)    =  SDo(1)*o(2)-SDO(2)*o(1)

    t0  =  3.0d0*(S(1,1)**2+S(2,2)**2+S(3,3)**2)
    t1  = (S(1,1)+S(2,2)+S(3,3))**2
    tilt=  sqrt(6.0d0)*sqrt(v(1)**2+v(2)**2+v(3)**2)/(oo*oo*sqrt(t1-t0)) &
        & *max(1.0d0, 0.2d0*nu/nuT)

    return
    end subroutine get_vortex_tilting
