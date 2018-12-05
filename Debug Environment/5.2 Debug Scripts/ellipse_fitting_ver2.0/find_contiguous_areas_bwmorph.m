function [blob_coord, blob_numpixels, crash_pts, flag] = find_contiguous_areas_bwmorph(img, LS_USA, min_border_pix, x_earth, y_earth, imin, jmin, deg2met)

IM = img < 150;
if ismember(1,IM) == 0
    flag = 1; %no low pop zones
elseif ismember(0,IM) == 0
    flag = 2; %all low pop zones
else
    flag = 0;
end

if flag == 0;
    s = regionprops(IM, 'Area', 'Perimeter', 'PixelList', 'Centroid');
    x_coord = [];
    y_coord = [];
    x_centroid = [];
    y_centroid = [];
    for i = 1:length(s)
        x_cent = s(i).Centroid(1);
        y_cent = s(i).Centroid(2);
        blob_area(i) = s(i).Area;
        blob_outline(i) = s(i).Perimeter;
        x_coord = [x_coord s(i).PixelList(:,1)'];
        y_coord = [y_coord s(i).PixelList(:,2)'];
        yrow = round(y_cent, 0) + min(imin) - 1;
        xcol = round(x_cent, 0) + min(jmin) - 1;
        ME2 = LS_USA(:,:,2); % check
        ME3 = LS_USA(:,:,3);
        xnew1 = ME2(sub2ind(size(ME2), yrow, xcol));
        xnew = xnew1(1);
        ynew1 = ME3(sub2ind(size(ME3), yrow, xcol));
        ynew = ynew1(1);
        dist(i) = sqrt((x_earth - xnew)^2 + (y_earth - ynew)^2)*deg2met;
        x_centroid = [x_centroid xnew];
        y_centroid = [y_centroid ynew];
    end

    blob_numpixels_temp = [blob_area(:) blob_outline(:)];
    row_remove = [];
    for k = 1:size(blob_numpixels_temp,1);
        if blob_numpixels_temp(k,1) < 4;
            row_remove(k) = k;
        end
    end
    row_remove(row_remove == 0) = [];

    blob_numpixels_temp(blob_numpixels_temp(:,1) < 4, :) = [];
    blob_numpixels = blob_numpixels_temp;

    blob_coord_temp = zeros(length(s), length(y_coord), 2);
    for j = 1:length(s)
        x_coord_temp = s(j).PixelList(:,1)';
        y_coord_temp = s(j).PixelList(:,2)';
        blob_coord_temp(j,1:length(y_coord_temp),1) = y_coord_temp;
        blob_coord_temp(j,1:length(x_coord_temp),2) = x_coord_temp;
    end

    blob_coord_temp2 = blob_coord_temp;
    blob_coord_temp2([row_remove], :, :) = [];
    blob_coord = blob_coord_temp2;
    
    crash_pts = [x_centroid(:) y_centroid(:) blob_area(:) dist(:)]; %high_low_rat not right
    
elseif flag == 2
    blob_coord = [];
    blob_numpixels = []; 
    s = regionprops(IM, 'Centroid', 'Area');
    x_centr = s.Centroid(1);
    y_centr = s.Centroid(2);
    
    yrow = round(y_centr, 0) + min(imin) - 1;
    xcol = round(x_centr, 0) + min(jmin) - 1;
    ME2 = LS_USA(:,:,2); % check
    ME3 = LS_USA(:,:,3);
    xnew1 = ME2(sub2ind(size(ME2), yrow, xcol));
    xnew = xnew1(1);
    ynew1 = ME3(sub2ind(size(ME3), yrow, xcol));
    ynew = ynew1(1);
    area = s.Area;
    
    crash_pts = [xnew ynew area 0];
    
else
    blob_coord = [];
    blob_numpixels = [];
    crash_pts = [];
end