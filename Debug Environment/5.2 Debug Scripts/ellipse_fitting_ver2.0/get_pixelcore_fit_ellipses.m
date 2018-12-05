function [ra_, rb_, ang_, centroid_, nodes, endpts, area_, high_low_rat_, dist_, A_thin, xnew_, ynew_] = get_pixelcore_fit_ellipses ( blob,rb_initval, pop_map, x_earth, y_earth, imin, jmin, LS_USA, deg2met)
dist_crit_sm = 150; %150 tc3, 50 tc6 %75 tcs 1,6 B-matrix  %a larger number means fewer overlapping ellipses
dist_crit_lg = 35;
%% Initialize variables.
hi_pop = 100; %high population color
c_start = 40; %initial blob color
% c_incr = 1; %color increment for heat map
margin = 3; %enclose the pixels in a box
szi = max(blob, [], 1) + margin; %max values of y and x for window surrounding blob
A = hi_pop*ones(szi(1),szi(2)); %initialize image to all yellow.
A2 = A; %start with a blank skeleton

%% Generate baseline pixel map, and calculate minimum distance of each blob pixel to the edge.
szb = size(blob); [perimeter, numpixels] = find_perimeter( blob );
for ii=1:szb(1)
    %populate image array from blob pixel list
    A(blob(ii,1), blob(ii,2)) = c_start;
    %calculate distance from blob pixel to nearest border pixel.
    for jj=1:numpixels(2)
        dist(jj) = norm(blob(ii,:)-perimeter(jj,:));
    end
end
A(1,:) = hi_pop; A(:,1) = hi_pop; %2/15 bugfix for blobs that go off edge
A_orig = A; %original unmodified blob array for later use

%% Thin the core and obtain its nodes and endpoints.
A_thin = thin_skeleton(A_orig); % need skel for each blob, not entire image
[nodes, endpts, A_node_thin] = find_nodes( A_thin, -9999, 100, 35, 40); % good?
if size(nodes,1)==0 && size(endpts,1) == 0
    ra_=[]; rb_=[]; ang_=[]; centroid_=[]; nodes=[]; area_ = []; high_low_rat_ = []; dist_ = []; xnew_ = []; ynew_ = []; endpts = [];
    status = -1;
    return;
end

%% Fit ellipses to the region.
% Fit ellipses to endpoint-to-node line segments
num_endpts = size(endpts,1);
tot_pts = [endpts(:,:); nodes(:,:)];
idx_e = 1;

for idx=1:num_endpts
    for idx2 = 1:size(tot_pts,1)
        endpts2 = [endpts(idx,:); tot_pts(idx2,:)]; %2x2
        %check connectedness of nodes
        [lineseg, szl] = generate_line_pixels( endpts(idx,:), tot_pts(idx2,:) );
        neighbors(idx2) = 1;
        for idx_px=1:szl
            if A_orig( lineseg(idx_px,1), lineseg(idx_px,2) ) == hi_pop %if any pixels = 100 on linesegment, neighbor = 0
                neighbors(idx2) = 0;
            end
        end
    %% what if neighbors = 0? could maybe have a workaround for that 
        if neighbors(idx2) == 1 %if no high population pixels on line segment, grow ellipses
            [ ra(idx_e,1),rb(idx_e,1),ang(idx_e,1),centroid(idx_e,:) ] = grow_ellipse( endpts2, A_orig, hi_pop, rb_initval, 1, 0 ); %look into changing 0s to 1s to enable internal ellipses and concurrent growth
            if ra(idx_e,1)<rb(idx_e,1) %debug for when major axis is shorter than minor axis
                    [ ra(idx_e,1),rb(idx_e,1),ang(idx_e,1) ] = ellipse_flip_a_b_theta( ra(idx_e,1),rb(idx_e,1),ang(idx_e,1) );
            end
            idx_e = idx_e+1;
        end
    end
end
endpts_ = endpts; 
ellipse_leaves = [ra rb centroid(:,2) centroid(:,1) ang zeros(length(ra),4)]; %0's are placeholders
num_ellipses_endpt = length(ra);

% Fit ellipses to node-to-node line segments
szn = size(nodes,1);
idx_n = 1; %node ellipse number used in method 2 only
nodes_2 = zeros(szn,4); %list of connected internal node-pairs

for idx2=1:szn
    %draw ellipses connecting the current node to all its reachable neighbors
    for idx2_1=1:szn
        %determine all neighbors (nodes joinable by a lowpop only line)
        if idx2 ~= idx2_1
            [lineseg, szl] = generate_line_pixels( nodes(idx2,:), nodes(idx2_1,:) );
            neighbors = 1;
            for idx_px=1:szl
                if A_orig( lineseg(idx_px,1), lineseg(idx_px,2) ) == hi_pop
                    neighbors = 0;
                end
            end
            if neighbors == 1
                endpts2 = [nodes(idx2,:); nodes(idx2_1,:)]; %mx2x2 matrix
                endpts2r = [nodes(idx2,:) nodes(idx2_1,:)]; %mx4x1 matrix
                % grow an ellipse and add the endpoint pair to the "do not duplicate" register
                no_dup = 1;
                for idx_dup=1:size(nodes_2,1)
                    if sum(abs(endpts2r(1,:) - nodes_2(idx_dup,:)),2)<0.1
                        no_dup = 0;
                    end
                    if sum([abs(endpts2r(1,1:2) - nodes_2(idx_dup,3:4)) abs(endpts2r(1,3:4) - nodes_2(idx_dup,1:2))],2)<0.1
                        no_dup = 0;
                    end
                end
                if no_dup % if no_dup = 1
                    [ ra(num_ellipses_endpt+idx_n,1),rb(num_ellipses_endpt+idx_n,1),ang(num_ellipses_endpt+idx_n,1),centroid(num_ellipses_endpt+idx_n,:) ] = grow_ellipse( endpts2, A_orig, hi_pop, 0, 1, 1 );
                    if ra(num_ellipses_endpt+idx_n,1)<rb(num_ellipses_endpt+idx_n,1)
                        [ ra(num_ellipses_endpt+idx_n,1),rb(num_ellipses_endpt+idx_n,1),ang(num_ellipses_endpt+idx_n,1) ] = ellipse_flip_a_b_theta( ra(num_ellipses_endpt+idx_n,1),rb(num_ellipses_endpt+idx_n,1),ang(num_ellipses_endpt+idx_n,1) );
                    end
                    nodes_2(num_ellipses_endpt+idx_n,:) = [nodes(idx2,:) nodes(idx2_1,:)];
                    idx_n = idx_n + 1;
                    idx = idx + 1;
                end
            end
        end
    end
end
%%
num_ellipses_node = idx_n-1;
beg_idx=num_ellipses_endpt+1;
end_idx=num_ellipses_endpt+num_ellipses_node;
ellipse_branches = [ra(beg_idx:end_idx,1) rb(beg_idx:end_idx,1) centroid(beg_idx:end_idx,2) centroid(beg_idx:end_idx,1) ang(beg_idx:end_idx,1) nodes_2(beg_idx:end_idx,:)];
ellipse_candidates = [ellipse_leaves; ellipse_branches];
idx_=1; szec = size(ellipse_candidates,1);

ra_ = [];
rb_ = [];
centroid_ = [];
ang_ = [];
nodes_ = [];
for iii=1:szec
    if  ellipse_candidates(iii,1) && ellipse_candidates(iii,2)
        ra_(idx_,1) = ellipse_candidates(iii,1);
        rb_(idx_,1) = ellipse_candidates(iii,2);
        centroid_(idx_,2) = ellipse_candidates(iii,3);
        centroid_(idx_,1) = ellipse_candidates(iii,4);
        ang_(idx_,1) = ellipse_candidates(iii,5);
        nodes_(idx_,:) = ellipse_candidates(iii,6:9);
        idx_ = idx_ + 1;
    end
end

%% Test if an ellipse is inside another ellipse and remove.
[ra_, rb_, centroid_, ang_, area_, high_low_rat_, dist_, xnew_, ynew_] = final_ellipse_params(ra_, rb_, centroid_, ang_, pop_map, x_earth, y_earth, imin, jmin, LS_USA, deg2met);
area_ = area_(:);
high_low_rat_ = high_low_rat_(:);
dist_ = dist_(:);
xnew_ = xnew_(:);
ynew_ = ynew_(:);
status = 0;
end
