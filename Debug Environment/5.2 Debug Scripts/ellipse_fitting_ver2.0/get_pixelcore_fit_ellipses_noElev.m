function [ra_, rb_, ang_, centroid_, nodes, endpts, area_, high_low_rat_, dist_, A_thin, xnew_, ynew_] = get_pixelcore_fit_ellipses_noElev ( blob,rb_initval, pop_map, x_earth, y_earth, imin, jmin, LS_USA, deg2met, FGIF_p, FGIF_start, d_deg, psi0)
%% Initialize variables.
high_pop_val = 10000;
hi_pop = 10000; %high population color
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
% Get FGIF integer pts and see if oob
FGIF_p45i = [round(FGIF_p(:,1),0) round(FGIF_p(:,2),0)];
FGIF_p45i = unique(FGIF_p45i, 'rows');
FGIF_p45i(1:3:end,:) = [];
% get set of out-of-bounds points (high population)
centroid_opts1 = [];
angle = [];
for i = 1:size(FGIF_p45i,1)
    if pop_map(FGIF_p45i(i,2), FGIF_p45i(i,1)) == high_pop_val 
        continue
    else
        dy = FGIF_p45i(i,2) - FGIF_start(1,2);
        dx = FGIF_p45i(i,1) - FGIF_start(1,1);
        angle = [angle; wrapTo2Pi(atan2(dy,dx))];
        centroid_opts1 = [centroid_opts1; FGIF_p45i(i,2) FGIF_p45i(i,1)];
    end
end
centroid_opts = [centroid_opts1(:,:) angle(:)];
if  size(centroid_opts,1) == 0
    ra_=[]; rb_=[]; ang_=[]; centroid_=[]; nodes=[]; area_ = []; high_low_rat_ = []; dist_ = []; xnew_ = []; ynew_ = []; endpts = []; 
    status = -1;
    return;
end

%% Fit ellipses to the region.
% Fit ellipses to endpoint-to-node line segments
ra = [];
rb = [];
ang = [];
centroid = [];
for i = 1:size(centroid_opts,1)
    endpts2 = gen_endpts_noElev(centroid_opts(i,:), blob, pop_map, high_pop_val, psi0);
    % obtain endpts not oob within the given blob, then use grow ellipse to
    % expand the minor AND major axes
    if size(endpts2,1) == 2
        [ ra_pot, rb_pot, ang_pot, centroid_pot ] = grow_ellipse_noElev( endpts2, pop_map, high_pop_val, rb_initval, 1, 0); 
    else
        ra_pot = 0; rb_pot = 0; ang_pot = 0; centroid_pot = [0 0];
    end
    ra = [ra; ra_pot];
    rb = [rb; rb_pot];
    ang = [ang; ang_pot];
    centroid = [centroid; centroid_pot];
end
if  nnz(ra) == 0
    ra_=[]; rb_=[]; ang_=[]; centroid_=[]; nodes=[]; area_ = []; high_low_rat_ = []; dist_ = []; xnew_ = []; ynew_ = []; endpts = []; 
    status = -1;
    return;
end
ellipse_candidates = [ra(:) rb(:) centroid(:,2) centroid(:,1) angle(:) zeros(length(ra),4)]; %0's are placeholders
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
[ra_, rb_, centroid_, ang_, area_, high_low_rat_, dist_, xnew_, ynew_] = final_ellipse_params_noElev(ra_, rb_, centroid_, ang_, pop_map, x_earth, y_earth, imin, jmin, LS_USA, deg2met, FGIF_p45i);
area_ = area_(:);
high_low_rat_ = high_low_rat_(:);
dist_ = dist_(:);
xnew_ = xnew_(:);
ynew_ = ynew_(:);
status = 0;
end
