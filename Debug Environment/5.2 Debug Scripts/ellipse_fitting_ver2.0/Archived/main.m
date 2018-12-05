%Load generated mat-file.
clear all
%% Andrew's Notes %%
% decode_data: sets up compressed LS_USA data inbto useable matrix of 0's 
% and 1's

% minlat, maxlat, minlong, maxlong: parse large LS_USA matrix by specifying
% index ranges

% extract_local_data: parse original LS_USA matrix into smaller map based
% on specified lat and long ranges

% img_cleanup: fill singleton pixels or spursize blocks to smooth out
% image. Replaces those pixels with high population value of 1? Sometimes
% outputs plot that smooths out the image. 

% find_contiguous_areas: find border of pixels of same type next to each
% other. Provides bordergroup (borders of each blob), group (combination of
% pixels classified as a group), centoid (center of each blob), and
% numpixels (total number of pixels in each blob?)

% reshape: fits the actual ellipses to each blob. Provides ra (semimajor 
% axis radius), rb(semiminor axis radius ), ang (orientation of ellipse, 
% deg), centroid (centroid points of final ellipses), nodes ( ), and 
% status ( )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pop_thresh = 7;
tile_i=13; tile_j=8;
% minlat=1040; maxlat=minlat+50; minlong=300; maxlong=minlong+50; spursize=1; 
minlat=1110; maxlat=1160; minlong=100; maxlong=150; spursize=3; %finds ellipses ok

matfilename = strcat('bin_rle_thrsh_',int2str(pop_thresh),'_tile_',int2str(tile_i),'_',int2str(tile_j),'.mat');
map = decode_data(matfilename); map = map*100;
% map = flip(map,1); %pixel and latitude i-axes should be positive up

%Limit the area to 50x50
map_thresh_local = extract_local_data( map, minlat, maxlat, minlong, maxlong );

figure(1);
display_pixel_array( map_thresh_local );


%Filter the data to prepare it make ellipse fitting possible.
map_thresh_local_filtered = img_cleanup(map_thresh_local, spursize, 100);
figure(2);
display_pixel_array( map_thresh_local_filtered );

% % Test regionprops 
% figure(3)
% BW = map_thresh_local_filtered < 1;
% BW_orig = map_thresh_local < 1;
% 
% imagesc(BW_orig)
% hold on
% set(gca,'YDir','normal')
% cmap = [0.5843 0.8157 0.9882; 1 1 1];
% colormap(cmap);
% s = regionprops(BW, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
% phi = linspace(0, 2*pi, 50);
% cosphi = cos(phi);
% sinphi = sin(phi);
% if size(s) > 0
%     for k = 1:size(s)
%         xbar = s(k).Centroid(1);
%         ybar = s(k).Centroid(2);
%         a = s(k).MajorAxisLength/2;
%         b = s(k).MinorAxisLength/2;
%         ellipse_area(k) = pi*s(k).MajorAxisLength/2*s(k).MinorAxisLength/2;
%         theta = pi*s(k).Orientation/180;
%         R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
% 
%         xy = [a*cosphi; b*sinphi];
%         xy = R*xy;
% 
%         x = xy(1,:) + xbar; % Ellipse coordinates
%         y = xy(2,:) + ybar;
%         plot(x, y, 'y', 'LineWidth', 2)
%         plot(xbar, ybar, '*b', 'Markersize', 10)
%     end
% end
% hold off

%Get the border pixels, pixel-blocks, centroids, and sizes!
%Needs replacing with a much faster MATLAB function.
[blob_outlines, blobs, blob_centroids, blob_numpix] = find_contiguous_areas( map_thresh_local_filtered ); %just need blobs and blobs_numpix


%Obtain and display the ellipses...
figure(4);
tic
display_pixel_array(map_thresh_local); hold on;
for idx=1:size(blobs,1)
    tmpblob = blobs(idx,1:blob_numpix(idx,1),:);
    thisblob = reshape(tmpblob,size(tmpblob,2),size(tmpblob,3));
    [ra, rb, ang, centroid, nodes, status] = get_pixelcore_fit_ellipses ( thisblob, 2 );
end
aa = toc
