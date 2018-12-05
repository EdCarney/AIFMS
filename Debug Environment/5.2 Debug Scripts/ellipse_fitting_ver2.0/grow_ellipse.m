function [ ra_safe,rb_safe,ang,centroid_safe ] = grow_ellipse( endpts, bin_image, oob_color, rb_init, dbl_ended_growth, concurrent_growth )
%grow_ellipse: form ellipses from endpoints of a branch and tune
%   Input - endpts of ellipse (2x2), binary image of current pixel-block
%   and surrounding out-of-bounds area (mxn)
%   Output - best fit ellipse for the endpoints that don't violate the OOB areas
%   Note that the first endpoint is the direction of ellipse "growth",
%   while the second endpoint is fixed for the most part.
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
    if ~exist('rb_init','var')
        rb = 2; %arbitrary
    else
        rb = rb_init; %user specified value
    end
    
    % get initial oob state of ellipse
    [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpts(1,:), endpts(2,:), rb, 'r', bin_image );    
    [oob_flg, encr_pt] = detect_encroaching_points( elpse, oob );
    % decrease the semi-major & semi-minor axes until ellipse is in-bounds
    % counts are for avoiding infinite loops
    count1 = 0;
    count2 = 0;
    count3 = 0;
    flag = 0;
    while (oob_flg == 1)
        if count1 < 25
            count1 = count1 + 1;
            endpts(1,:) = endpts(1,:) - [grow_len*sin(ang) grow_len*cos(ang)];
            endpts(2,:) = endpts(2,:) + [grow_len*sin(ang) grow_len*cos(ang)];
            rb = max(0, rb - 1);
            [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpts(1,:), endpts(2,:), rb, 'r', bin_image );    
            [oob_flg, encr_pt] = detect_encroaching_points( elpse, oob );
        else
            flag = 1;
            disp('broken loop 1')
            break
        end
    end
    if flag == 1
        ra_safe = 0;
        rb_safe = 0;
        ang = 0;
        centroid_safe = [0 0];
        return
    end

    % increase the semi-major & semi-minor axes until end of ellipse goes out-of-bounds
    while (oob_flg == 0) || (tune < tune_itns)
        if count2 < 20
            count2 = count2 + 1;
            % get points of ellipse from endpts
            [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpts(1,:), endpts(2,:), rb, 'r', bin_image );
            %check if it went out-of-bounds and extend the major axis    
            [oob_flg, encr_pt] = detect_encroaching_points( elpse, oob );
            if oob_flg == 0
                endpts_safe = endpts; %record ellipse if no encroachment
                centroid_safe = centroid;
                ra_safe = ra;
                rb_safe = rb;
            else
                if tune == 1
                    %done with tuning major axis, now tune minor axis
                    oob_flg = 0;
                    tune = tune+1;
                    endpts(1,:) = endpts(1,:) - 2*[grow_len*sin(ang) grow_len*cos(ang)]; %add a small space
                end
            end
            % actually grow the semimajor or semiminor axis
            if tune == 1 %one-sided or double-sided growth
                endpts(1,:) = endpts(1,:) + [grow_len*sin(ang) grow_len*cos(ang)];
                if dbl_ended_growth
                    endpts(2,:) = endpts(2,:) - [grow_len*sin(ang) grow_len*cos(ang)];
                end
                if concurrent_growth
                    rb = rb + 1;
                end
            else
                rb = rb + 1;
            end
        else 
            disp('broken loop 2')
            flag = 1;
            break
        end
    end
    
    if flag == 1
        ra_safe = 0;
        rb_safe = 0;
        ang = 0;
        centroid_safe = [0 0];
        return
    end
    %shift centroid in the direction away from semi-minor axis encroachment
    oob_flg = 0;
    if debug == 0
        rb = rb_safe;
        ra = ra_safe;
        centroid = centroid_safe;
    else %turn safety off so program doesn't crash
        rb_safe = rb-1;
        ra_safe = ra-1;
        centroid_safe = centroid;
    end
    shft_uvec = [sin(ang+pi/2) cos(ang+pi/2)]; %[i j] or [y x] unitvec
    
    %which side of ellipse centerline is the encroaching point and unitvec
    %endpoint on?
    side_encr = (encr_pt(2)-endpts(1,2)) * (endpts(2,1)-endpts(1,1)) - (encr_pt(1)-endpts(1,1)) * (endpts(2,2)-endpts(1,2));
    side_uvec = (shft_uvec(2)-endpts(1,2)) * (endpts(2,1)-endpts(1,1)) - (shft_uvec(1)-endpts(1,1)) * (endpts(2,2)-endpts(1,2));
    if side_encr*side_uvec > 0
        shft_uvec = -shft_uvec;
    end

    % shift minor axis until a collision, increasing semi-minor axis by the amount shifted
    while (oob_flg == 0)
        if count3 < 20
            count3 = count3 + 1;
            centroid = centroid + shft_uvec;  %performing the centroid shift...
            rb = rb + 1;
            % note: cleanup below repeated code at some point
            [x,y]=ellipse(ra,rb,ang,centroid(2),centroid(1));
            numpts = size(x,2);
            xprev = -99999; yprev = -99999;
            idx_e = 1;
            %pixelize the ellipse (int)
            clear elpse;
            for idx=1:numpts
                xcur = round( x(idx) ); ycur = round( y(idx) );
                if (xcur ~= xprev) || (ycur ~= yprev)
                    elpse(idx_e,:) = [ycur xcur];
                    idx_e = idx_e + 1;
                end
                xprev = xcur; yprev = ycur;
            end
            %check if it went out-of-bounds    
            [oob_flg, encr_pt] = detect_encroaching_points( elpse, oob );
            if oob_flg == 0
                endpts_safe = endpts; %record ellipse if no encroachment
                centroid_safe = centroid;
                ra_safe = ra;
                rb_safe = rb;        
            end
        else
            disp('broken loop 3')
            flag = 1;
        end
    end
    if flag == 1
        ra_safe = 0;
        rb_safe = 0;
        ang = 0;
        centroid_safe = [0 0];
        return 
    end
end
