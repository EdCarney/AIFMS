function [ra_, rb_, ang_, centroid_, nodes, endpts, area_, high_low_rat_, dist_, A_thin, xnew_, ynew_] = get_pixelcore_fit_ellipses_noturn ( blob,rb_initval, pop_map, x_earth, y_earth, imin, jmin, LS_USA, deg2met, FGIF_p, d_deg)
dist_crit_sm = 150; %150 tc3, 50 tc6 %75 tcs 1,6 B-matrix  %a larger number means fewer overlapping ellipses
dist_crit_lg = 35;

%% Initialize variables.
high_pop_val = 10000;
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
[nodes, endpts, ~] = find_nodes( A_thin, -9999, 100, 35, 40); % good?

%% Restrain FGIF pts to given blob and extract potential endpts
endpts = gen_endpts(FGIF_p, blob, pop_map, high_pop_val);
if  size(endpts,1) == 0
    ra_=[]; rb_=[]; ang_=[]; centroid_=[]; nodes=[]; area_ = []; high_low_rat_ = []; dist_ = []; xnew_ = []; ynew_ = []; endpts = []; 
    status = -1;
    return;
end

%% Fit ellipses to the region.
% Fit ellipses to endpoint-to-node line segments
endpts2 = [endpts(1,:); endpts(end,:)];
[ ra, rb, ang, centroid(1,:) ] = grow_ellipse_noturn( endpts2, endpts, pop_map, high_pop_val, rb_initval, 1, 0 ); %look into changing 0s to 1s to enable internal ellipses and concurrent growth

iter = 0;
if rb == 0
    flag1 = 1;
    flag2 = 1;
else
    flag1 = 0;
    flag2 = 0;
end
while flag1 == 1
    iter = iter + 1;
    if iter >= (size(endpts,1) - 1)
        ra1 = [];
        rb1 = [];
        ang1 = [];
        centroid1 = [];
        break
    else
        endpts2 = [endpts(1+iter,:); endpts(end,:)];
    end
    [ ra1, rb1, ang1, centroid1(1,:) ] = grow_ellipse_noturn( endpts2, endpts, pop_map, high_pop_val, rb_initval, 1, 0 ); %look into changing 0s to 1s to enable internal ellipses and concurrent growth

    if rb1 == 0
        flag1 = 1;
    else
        flag1 = 0;
    end
end

iter2 = 0;
while flag2 == 1
    iter2 = iter2 + 1;
    if iter2 >= (size(endpts,1) - 1)
        ra = [];
        rb = [];
        ang = [];
        centroid = [];
        break
    else
        endpts2 = [endpts(1,:); endpts(end-iter2,:)];
    end
    [ ra2, rb2, ang2, centroid2(1,:) ] = grow_ellipse_noturn( endpts2, endpts, pop_map, high_pop_val, rb_initval, 1, 0 ); %look into changing 0s to 1s to enable internal ellipses and concurrent growth
    if rb2 == 0
        flag2 = 1;
    else
        flag2 = 0;
    end
end

if rb == 0
    ra = [ra1; ra2];
    rb = [rb1; rb2];
    ang = [ang1; ang2];
    centroid = [centroid1(:,:); centroid2(:,:)];
end

ellipse_candidates = [ra rb centroid(:,2) centroid(:,1) ang zeros(length(ra),4)]; %0's are placeholders
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
[ra_, rb_, centroid_, ang_, area_, high_low_rat_, dist_, xnew_, ynew_] = final_ellipse_params_noturn(ra_, rb_, centroid_, ang_, pop_map, x_earth, y_earth, imin, jmin, LS_USA, deg2met, FGIF_p);
area_ = area_(:);
high_low_rat_ = high_low_rat_(:);
dist_ = dist_(:);
xnew_ = xnew_(:);
ynew_ = ynew_(:);
status = 0;
end
