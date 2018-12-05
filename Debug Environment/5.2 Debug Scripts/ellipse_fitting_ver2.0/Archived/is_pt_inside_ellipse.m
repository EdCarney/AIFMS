function [ pt_in_ellipse, dist_to_foci ] = is_pt_inside_ellipse( a,b,h,k,ang,x,y )
%is_pt_inside_ellipse: returns whether point is inside ellipse and sum of
%distances to foci
    focal_r = (a^2-b^2)^.5;
    focus(1,:) = [k+focal_r*sin(ang) h+focal_r*cos(ang)];
    focus(2,:) = [k-focal_r*sin(ang) h-focal_r*cos(ang)];
    dist_to_foci = 0;
    for idx2=1:2
        thisdistvec = [y x]-focus(idx2,:);
        thisdist = norm(thisdistvec,2);
        dist_to_foci = dist_to_foci + thisdist;
        a;
    end
    if dist_to_foci<=(2*a)
        pt_in_ellipse = 1;
    else
        pt_in_ellipse = 0;
    end
end

