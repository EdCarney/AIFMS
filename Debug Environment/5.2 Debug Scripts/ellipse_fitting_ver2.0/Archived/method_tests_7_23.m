%% Test Other Polygon Fitting Methods %%
clear all

filename = 'bin_rle_thrsh_7_tile_13_8';
minlat=1110; maxlat=1160; minlong=100; maxlong=150; spursize=3; %finds ellipses ok
map = decode_data(filename); map = map*100;
parsed_map = extract_local_data( map, minlat, maxlat, minlong, maxlong );
parsed_map_filtered = img_cleanup(parsed_map, spursize, 100);

figure(1)
BW_parsed = parsed_map < 1;
imagesc(BW_parsed)
hold on 
set(gca,'YDir','normal')
cmap = [0.5843 0.8157 0.9882; 1 1 1];
colormap(cmap);
hold off

figure(2);
BW_parsed_filtered = parsed_map_filtered < 1;
imagesc(BW_parsed_filtered)
hold on 
set(gca,'YDir','normal')
cmap = [0.5843 0.8157 0.9882; 1 1 1];
colormap(cmap);
hold off