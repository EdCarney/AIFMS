%Load generated mat-file.
clear all
tic
cellsize = 0.00083333333333334; %pixels to latitude or longitude
xllcorner = -124.99999999999;
yllcorner = 24.25;
pop_thresh = 7;
tile_i=8; tile_j=50;

llcorner_this_tile = convert_pixels2latlon( [30300 69900],[20 60],[tile_i tile_j],[1 1], cellsize, [yllcorner xllcorner], [0 0] );

matfilename = strcat('bin_rle_thrsh_',int2str(pop_thresh),'_tile_',int2str(tile_i),'_',int2str(tile_j),'.mat');
map = decode_data(matfilename); map = map*100;

%Limit the area to 50x50
minlat = 1; maxlat = size(map,1); minlong = 1; maxlong = size(map,2); spursize = 3;
map_processed = extract_local_data( map, minlat, maxlat, minlong, maxlong );

save('map_processed_7_8_50.mat', 'map_processed', 'cellsize', 'llcorner_this_tile')

%Filter the data to prepare it make ellipse fitting possible.
map_processed_filtered = img_cleanup(map_processed, spursize, 100);

save('map_processed_filtered_7_8_50.mat', 'map_processed_filtered')