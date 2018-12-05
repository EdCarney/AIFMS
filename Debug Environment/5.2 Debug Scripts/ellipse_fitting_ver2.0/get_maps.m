%Load generated mat-file.
clear all
tic
pop_thresh = 7;
tile_i=13; tile_j=8;
min_border_pix = 4;

matfilename = strcat('bin_rle_thrsh_',int2str(pop_thresh),'_tile_',int2str(tile_i),'_',int2str(tile_j),'.mat');
map = decode_data(matfilename); map = map*100;

% %Limit the area to 50x50
% map_thresh_local = extract_local_data( map, minlat, maxlat, minlong, maxlong );
% 
% %Filter the data to prepare it make ellipse fitting possible.
% map_thresh_local_filtered = img_cleanup(map_thresh_local, spursize, 100);
% figure(1);
% display_pixel_array( map_thresh_local_filtered );
% 
% %Get the border pixels, pixel-blocks, centroids, and sizes!
% [blobs, blob_numpix, crash_options, flag] = find_contiguous_areas_bwmorph(map_thresh_local_filtered, min_border_pix);
% crash_pts = sortrows(crash_options, 3, 'descend');
% figure(3)
% display_pixel_array(map_thresh_local_filtered);
% hold on 
% for j = 1:size(blobs, 1)
%     x_out = blobs(j,:,2);
%     y_out = blobs(j,:,1);
%     plot(x_out, y_out, '.c', 'Markersize', 10);
% end 
% hold off
% 
% if flag == 1
%     disp('No low population zones!')
%     land_pts = [];
% elseif flag == 2
%     disp('All low population zones!')
%     land_pts = [];
% else
%     %Obtain the ellipses...
%     figure(2)
%     display_pixel_array(map_thresh_local); hold on;
%     ellipses = []; delimiter = -1;
%     for idx=1:size(blobs,1)
%         tmpblob = blobs(idx,1:blob_numpix(idx,1),:);
%         thisblob = reshape(tmpblob,size(tmpblob,2),size(tmpblob,3));
%         [ra, rb, ang, centroid, nodes, area, high_low_rat] = get_pixelcore_fit_ellipses ( thisblob, 0.5, map_thresh_local);
%         this_ellipse_set = [ra rb ang centroid area high_low_rat];
%         if size(this_ellipse_set,1)
%             tmp = ellipses;
%             ellipses = [tmp; this_ellipse_set];
%         end
%     end
%     
%     if size(ellipses, 1) > 0
%         land_pts = [ellipses(:,5) ellipses(:,4) ellipses(:,3).*180/pi ellipses(:,6) ellipses(:,7)]; %[x, y, orientation, area, high:low pop]
%     else 
%         land_pts = [];
%     end
%     land_pts = sortrows(land_pts, [5 4], {'ascend' 'descend'});
%     toc
% 
% %     ...and display them
%     x=[]; y=[];
%     for idx=1:size(ellipses,1)
%         if ellipses(idx,1)>0
%             [x_cur,y_cur]=ellipse(ellipses(idx,1),ellipses(idx,2),ellipses(idx,3),ellipses(idx,5),ellipses(idx,4),'c');
%             tmp_x = x; tmp_y = y;
%             x = [tmp_x; x_cur']; y = [tmp_y; y_cur'];
%         end
%     end
%     plot(x,y,'c.','LineWidth',1);
%     hold off
% end
