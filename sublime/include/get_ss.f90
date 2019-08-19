    uc(1:5) =  fv(isec)%uc(1:5,iele)
    p       =  fv(isec)%u (  5,iele)

    lap     =  0.0d0
    lap_2   =  0.0d0
    ss      =  0.0d0
    ss_2    =  0.0d0
    do j=sec(isec)%iA_face_neighbour(iele),sec(isec)%iA_face_neighbour(iele+1)-1
        im  =  sec(isec)%jA_face_neighbour(j)
        if(im .gt. 0) then
            s   =  mesh(lev)%mortar_LR  (3, im)
            e   =  mesh(lev)%mortar_LR  (4, im)
            Ad  =  mesh(lev)%mortar_n_vg(4, im)*mesh(lev)%mortar_n_vg(5, im)
        else
            s   =  mesh(lev)%mortar_LR  (1,-im)
            e   =  mesh(lev)%mortar_LR  (2,-im)
            Ad  =  mesh(lev)%mortar_n_vg(4,-im)*mesh(lev)%mortar_n_vg(5,-im)
        end if
        if(sec(s)%is_bnd) then
            call u_to_uc(1, fv(s)%uR(1:5,e), uc_R, t_R, mu_R)
            pR  =  fv(s)%uR(5,e)
        else
            uc_R(1:5)   =  fv(s)%uc(1:5,e)
            pR          =  fv(s)%u (  5,e)
        end if

        lap     =  lap  +Ad*(uc_R(1:5)-uc(1:5))
        lap_2   =  lap_2+Ad
        ss      =  ss   +abs(pR-p)
        ss_2    =  ss_2 +(pR+p)
    end do
    fv(isec)%ss_lap(1:5,iele)   =  lap(1:5)/lap_2
    fv(isec)%ss_lap(6  ,iele)   =  abs(ss)/ss_2
