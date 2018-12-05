function [ closest_pt, distance, tilt_angle ] = get_closest_pts( curr_locn, list_of_pts, numpts )
%get_closest_pts: finds the <numpts> closest points within <list_of_pts> to <curr_locn>.
%   Detailed explanation goes here
    sz = size(list_of_pts,1);
    if numpts > sz
        numpts = sz;
        disp('warning: number of points to find is greater than size of list_of_pts...setting to size of list_of_pts');
    end
    skp_list = zeros(numpts,1);
    closest_pt = zeros(numpts,2);
    distance = zeros(numpts,1);
    tilt_angle = zeros(numpts,1);
    
    for itn=1:numpts
        min_dist = 100000000000000000000000;
        for idx=1:sz
            skipthispt = 0;
            for kk=1:(itn-1)
                if idx == skp_list(kk)
                    skipthispt = 1;
                end
            end
            if skipthispt == 0
                dist_vec = list_of_pts(idx,:) - curr_locn(:)'; %fix indexing bug causing wrong distances
                dist = norm(dist_vec, 2);
                if dist < min_dist
                    min_dist = dist;
                    closest_pt(itn,:) = list_of_pts(idx,:);
                    distance(itn) = dist;
                    closest_pt_idx = idx;
                end
            end
        end
        skp_list(itn) = closest_pt_idx;
        dist_vec = closest_pt(itn,:) - curr_locn;
        slope = dist_vec(1) / dist_vec(2);
        tilt_angle(itn) = atand(slope);
    end
end

