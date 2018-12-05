function [ ra_safe,rb_safe,ang,centroid_safe ] = grow_ellipse_noElev( endpts, bin_image, oob_color, rb_init, dbl_ended_growth, concurrent_growth )
%grow_ellipse: form ellipses from endpoints of a branch and tune
%   Input - endpts of ellipse (2x2), binary image of current pixel-block
%   and surrounding out-of-bounds area (mxn)
%   Output - best fit ellipse for the endpoints that don't violate the OOB areas
%   Note that the first endpoint is the direction of ellipse "growth",
%   while the second endpoint is fixed for the most part.

% NOTE: endpts produce decimals in some runs, I think thats bad. still get
% negative endpts for some cases 
    debug = 0; %if debug==1, "safe" values are not used
    grow_len = 1;
    if concurrent_growth
        tune_itns = 1;
    else
        tune_itns = 2;
    end

    % get set of out-of-bounds points (high population)
    oob = find_oob_areas( bin_image,oob_color );
        
    oob_flg = 0;
    tune = 1; %1 for major axis, 2 for minor axis
    rb = rb_init; %user specified value
    % get initial oob state of ellipse
    [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpts(1,:), endpts(2,:), rb, 'r', bin_image );    
    [oob_flg, encr_pt] = detect_encroaching_points( elpse, oob );
    % decrease the semi-major & semi-minor axes until ellipse is in-bounds
    % counts are for avoiding infinite loops
    count1 = 0;
    flag = 0;
    while (oob_flg == 1)
        count1 = count1 + 1;
        if count1 < 4
            endpts(1,:) = endpts(1,:) - [grow_len*sin(ang) grow_len*cos(ang)]; % signs used to be switched, look out for this
            endpts(2,:) = endpts(2,:) + [grow_len*sin(ang) grow_len*cos(ang)];
            [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpts(1,:), endpts(2,:), rb, 'r', bin_image );    
            [oob_flg, encr_pt] = detect_encroaching_points( elpse, oob );
        else
            flag = 1;
            disp('broken loop 1')
            oob_flg = 0;
        end
    end
    if flag == 1
        ra_safe = 0;
        rb_safe = 0;
        ang = 0;
        centroid_safe = [0 0];
        return
    end
    ra_safe = ra;
    rb_safe = rb;
    centroid_safe = centroid; 
    
    % increase endpt2 until oob
    cnt3 = 0;
    [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpts(1,:), endpts(2,:), rb, 'r', bin_image );
    [oob_flg, encr_pt] = detect_encroaching_points( elpse, oob );
    flg = 0;
    while (oob_flg == 0) && endpts(1,1) < (size(bin_image,1)-1) && endpts(1,1) >= 1 && endpts(2,1) < (size(bin_image,1)-1) && endpts(2,1) >= 1 && endpts(1,2) < (size(bin_image,2)-1) && endpts(1,2) >= 1  && endpts(2,2) < (size(bin_image,2)-1) && endpts(2,2) >= 1; 
        cnt3 = cnt3 + 1;
        [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpts(1,:), endpts(2,:), rb, 'r', bin_image );
        [oob_flg, encr_pt] = detect_encroaching_points( elpse, oob );
        if cnt3 > 3
            disp('broken loop 2')
            centroid_safe = centroid;
            ra_safe = ra;
            rb_safe = rb;
            oob_flg = 1;
            flg = 1;
        end
        if oob_flg == 0 
            centroid_safe = centroid;
            ra_safe = ra;
            rb_safe = rb_init;
            endpts(2,:) = endpts(2,:) - [grow_len*sin(ang) grow_len*cos(ang)];
        elseif oob_flg == 1 && flg == 0
            endpts(2,:) = endpts(2,:) + 2*[grow_len*sin(ang) grow_len*cos(ang)];
            centroid_safe = centroid;
            ra_safe = ra;
            rb_safe = rb_init;
            oob_flg = 1;
        end
    end
    
    % increase endpt1 until oob
    cnt2 = 0;
    [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpts(1,:), endpts(2,:), rb, 'r', bin_image );    
    [oob_flg, encr_pt] = detect_encroaching_points( elpse, oob );
    flg = 0;
    while (oob_flg == 0) && endpts(1,1) < (size(bin_image,1)-1) && endpts(1,1) >= 1 && endpts(2,1) < (size(bin_image,1)-1) && endpts(2,1) >= 1 && endpts(1,2) < (size(bin_image,2)-1) && endpts(1,2) >= 1  && endpts(2,2) < (size(bin_image,2)-1) && endpts(2,2) >= 1
        cnt2 = cnt2 + 1;
        [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpts(1,:), endpts(2,:), rb, 'r', bin_image );
        [oob_flg, encr_pt] = detect_encroaching_points( elpse, oob );
        if cnt2 > 3
            disp('broken loop 3')
            centroid_safe = centroid;
            ra_safe = ra;
            rb_safe = rb;
            oob_flg = 1;
            flg = 1;
        end
        if oob_flg == 0 
            centroid_safe = centroid;
            ra_safe = ra;
            rb_safe = rb_init;
            endpts(1,:) = endpts(1,:) + [grow_len*sin(ang) grow_len*cos(ang)];
        elseif oob_flg == 1 & flg == 0
            endpts(1,:) = endpts(1,:) - 2*[grow_len*sin(ang) grow_len*cos(ang)];
            centroid_safe = centroid;
            ra_safe = ra;
            rb_safe = rb_init; 
            oob_flg = 1;
        end
    end
    
    % increase minor axis until oob
    cnt4 = 0;
    flg = 0;
    [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpts(1,:), endpts(2,:), rb, 'r', bin_image );
    [oob_flg, encr_pt] = detect_encroaching_points( elpse, oob );
    while (oob_flg == 0)
        cnt4 = cnt4 + 1;
        [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpts(1,:), endpts(2,:), rb, 'r', bin_image );
        [oob_flg, encr_pt] = detect_encroaching_points( elpse, oob );
        if cnt4 > 3
            disp('broken loop 4')
            endpts_safe = endpts;
            centroid_safe = centroid;
            ra_safe = ra;
            rb_safe = rb;
            oob_flg = 1;
            flg = 1;
        end
        if oob_flg == 0 & rb < 1
            endpts_safe = endpts;
            centroid_safe = centroid;
            ra_safe = ra;
            rb_safe = rb;
            rb = rb + 1;
        elseif oob_flg == 0 & rb > 1
            endpts_safe = endpts;
            centroid_safe = centroid;
            ra_safe = ra;
            rb_safe = rb; 
            oob_flg = 1;
        elseif oob_flg == 1 & flg == 1
            endpts_safe = endpts;
            centroid_safe = centroid;
            ra_safe = ra;
            rb_safe = rb_init;
            rb = rb_init;
            oob_flg = 1;
        end
    end
end
