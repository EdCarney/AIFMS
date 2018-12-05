function [ inside,dist_outside_pct,pt_in_ellipse_tot ] = is_ellipse_inside_ellipse( a1,b1,h1,k1,t1,a2,b2,h2,k2,t2 )
%is_ellipse_inside_ellipse: checks if 1st ellipse is inside 2nd ellipse
    %Inputs
    %ep - 4x2 matrix of a,b endpoints of smaller ellipse
    %a,b,h,k,angle - parameters for larger ellipse
    inside = 0;
    pt_in_ellipse_tot = 0;
    dist_outside_pct = 0;
    %calculate ep
    ep = ellipse_abhkt_to_endpts( a1,b1,h1,k1,t1 );
    for idx=1:4
        [ pt_in_ellipse, dist_to_foci ] = is_pt_inside_ellipse( a2,b2,h2,k2,t2,ep(idx,2),ep(idx,1) );
        pt_in_ellipse_tot = pt_in_ellipse_tot + pt_in_ellipse;
        if dist_to_foci>2*a2
            if idx<=2
                divisor=a1;
            else
                divisor=b1;
            end
            dist_outside_pct = dist_outside_pct + (dist_to_foci-2*a2)*100/divisor;
        end
    end
    dist_outside_pct = dist_outside_pct/4;
    if pt_in_ellipse_tot == 4
        inside = 1;
    end
end

