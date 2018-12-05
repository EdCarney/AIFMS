function [ coords_latlon,coords_pixels_global_i_j ] = convert_pixels2latlon( total_pixels_i_j,num_tiles_i_j,tile_i_j,local_pixel_list_i_j, scalefactor_pixel2latlon, ll_corner_lat_lon, extr_offset )
%convert_pixels2latlon: converts pixel values to latitude/longitude given:
    % total_pixels_i_j - 1x2 vector of image dimensions in pixels
    % num_tiles_i_j - 1x2 vector of #subdivisions in i and j directions
    % tile_i_j - which rectangular subdivision we are looking at
    % local_pixel_list_i_j - which pixels within the current subdivision we want to convert to lat/lon.
    % scalefactor_pixel2latlon - ratio of length of a degree lat/lon to pixels
    % ll_corner_lat_lon - 1x2 vector of latitude and longitude of lower left corner tile's lower left corner.
    % extr_offset - 1x2 vector of additional pixel offsets when extracting a block from a tile
    tile_i_j(1) = num_tiles_i_j(1)-tile_i_j(1)+1; %tiles numbers are flipped
    offset = (total_pixels_i_j./num_tiles_i_j).*(tile_i_j-[1 1]) + extr_offset;
    coords_pixels_global_i_j = [local_pixel_list_i_j(:,1)+offset(1,1) local_pixel_list_i_j(:,2)+offset(1,2)];
    coords_latlon = coords_pixels_global_i_j.*scalefactor_pixel2latlon;
    coords_latlon = [coords_latlon(:,1)+ll_corner_lat_lon(1,1) coords_latlon(:,2)+ll_corner_lat_lon(1,2)]; 
end
