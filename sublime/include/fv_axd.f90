    n_nb=  sec(isec)%iA_face_neighbour(iele+1)-sec(isec)%iA_face_neighbour(iele)
    do i_nb=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
        j_nb=  i_nb-sec(isec)%iA_face_neighbour(iele)+1
        im  =  sec(isec)%jA_face_neighbour(i_nb)
        if(im .gt. 0) then
            s_nb        =  mesh(lev)%mortar_LR(3  ,im)
            e_nb        =  mesh(lev)%mortar_LR(4  ,im)
            e (1:5,j_nb)=  mesh(lev)%mortar_n_vg(1:5,im)
            vg(1:3,j_nb)=  mesh(lev)%mortar_n_vg(6:8,im)
        else
            s_nb        =  mesh(lev)%mortar_LR(1  ,-im)
            e_nb        =  mesh(lev)%mortar_LR(2  ,-im)
            e (1:3,j_nb)= -mesh(lev)%mortar_n_vg(1:3,-im)
            e (4:5,j_nb)=  mesh(lev)%mortar_n_vg(4:5,-im)
            vg(1:3,j_nb)=  mesh(lev)%mortar_n_vg(6:8,-im)
        end if
        uR(1:5,j_nb)=  fv(s_nb)%u  (1:5,e_nb)
        d (1:5,j_nb)=  fv(s_nb)%duc(1:5,e_nb)
        if(is_vis_cal)  mu_R(1,j_nb)=  fv(s_nb)%mu(1,e_nb)
        if(is_tur_cal)  mu_R(2,j_nb)=  fv(s_nb)%mu(2,e_nb)
    end do
    call u_to_akh(n_nb, uR, akh_R)

    sgn (1:n_nb)= -1.0d0
    efix(1:n_nb)=  LHS_efix
    call AprdV(n_nb, e, vg, uR, akh_R, sgn, efix, d, v)
    r(1:5)  = -fv(isec)%rhs(1:5,iele)
    do j_nb=1,n_nb
        r(1:5)  =  r(1:5)-v(1:5,j_nb)
    end do

    if(is_vis_cal) then
        uL(1:5) =  fv(isec)%u (1:5,iele)
        mu_L(1) =  fv(isec)%mu(1  ,iele)
        if(is_tur_cal)  mu_L(2) =  fv(isec)%mu(2,iele)

        call AprdV_vis(n_nb, e, uL, uR, mu_L, mu_R, d, v)
        do j_nb=1,n_nb
            r(2:5)  =  r(2:5)+v(2:5,j_nb)
        end do
    end if

    n_nb=  iele+iA(isec)-1
    fv(isec)%duc(1:5,iele)  =  p%D(1:5,5*n_nb-4)*r(1) &
                            & +p%D(1:5,5*n_nb-3)*r(2) &
                            & +p%D(1:5,5*n_nb-2)*r(3) &
                            & +p%D(1:5,5*n_nb-1)*r(4) &
                            & +p%D(1:5,5*n_nb  )*r(5)
