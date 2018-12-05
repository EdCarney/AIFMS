function [ status_0_or_1, first_encroachment] = detect_encroaching_points( encroaching_pts_candidates, encroached_pts_candidates )
%UNTITLED Summary of this function goes here
%   Inputs: mx2 array of possible encroaching pts, nx2 array of possible pts encroached upon
%   Outputs: array of actual encroaching pts
sz = size(encroaching_pts_candidates);
m = sz(1);
sz = size(encroached_pts_candidates);
n = sz(1);

status_0_or_1 = 0;
first_encroachment = [0 0];
for idx1=1:m
    for idx2=1:n
        if encroaching_pts_candidates(idx1,1) == encroached_pts_candidates(idx2,1) && encroaching_pts_candidates(idx1,2) == encroached_pts_candidates(idx2,2)
            status_0_or_1 = 1;
            first_encroachment = encroaching_pts_candidates(idx1,:);
            break;
        end
    end
    if status_0_or_1==1
        break;
    end
end

end

