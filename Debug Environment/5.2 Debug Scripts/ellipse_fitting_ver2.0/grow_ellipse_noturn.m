function [ ra_safe,rb_safe,ang,centroid_safe ] = grow_ellipse_noturn( endpts, endpts_tot, bin_image, oob_color, rb_init, dbl_ended_growth, concurrent_growth )
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
    oob = find_oob_areas_noturn( bin_image,oob_color );
    
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
    flag = 0;
    while (oob_flg == 1) % while ellipse is over high population 
        count1 = count1 + 1;
        if count1 < (size(endpts_tot,1) - 1)
            dist1 = sqrt((encr_pt(1,1) - endpts(1,1))^2 + (encr_pt(1,2) - endpts(1,2))^2);
            dist2 = sqrt((encr_pt(1,1) - endpts(2,1))^2 + (encr_pt(1,2) - endpts(2,2))^2);
            for i = 1:size(endpts_tot,1)
                dist2endpt(i) = sqrt((encr_pt(1,1) - endpts_tot(i,1))^2 + (encr_pt(1,2) - endpts_tot(i,2))^2);
            end
            [~,idx_endpts] = min(dist2endpt(:));
            if dist1 >= dist2                
                endpts(1,:) = endpts(1,:);
                endpts(2,:) = endpts_tot(idx_endpts(1)-1,:);
            elseif dist2 > dist1
                endpts(1,:) = endpts_tot(idx_endpts(1)+1,:);
                endpts(2,:) = endpts(2,:);
            end
            rb = max(0, rb - 1);
            [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpts(1,:), endpts(2,:), rb, 'r', bin_image );    
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
                if rb < 1
                    rb = rb + 1;
                else 
                    break
                end
            else 
                if rb_safe == 0
                    dist12 = sqrt((encr_pt(1,1) - endpts(1,1))^2 + (encr_pt(1,2) - endpts(1,2))^2);
                    dist22 = sqrt((encr_pt(1,1) - endpts(2,1))^2 + (encr_pt(1,2) - endpts(2,2))^2);
                     for i = 1:size(endpts_tot,1)
                        dist2endpt(i) = sqrt((encr_pt(1,1) - endpts_tot(i,1))^2 + (encr_pt(1,2) - endpts_tot(i,2))^2);
                    end
                    [~,idx_endpts] = min(dist2endpt(:));
                    if dist12 >= dist22
                        endpts(1,:) = endpts(1,:);
                        endpts(2,:) = endpts_tot(idx_endpts(1)-2,:);
                    elseif dist12 < dist22
                        endpts(1,:) = endpts_tot(idx_endpts(1)+2,:);
                        endpts(2,:) = endpts(2,:);
                    end
                else
                    break
                end
            end 
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
