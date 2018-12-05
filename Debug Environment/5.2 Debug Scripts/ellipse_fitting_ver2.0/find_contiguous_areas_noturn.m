function [blob_coord, blob_numpixels, crash_pts, flag] = find_contiguous_areas_noturn(img, LS_USA, min_border_pix, x_earth, y_earth, imin, jmin, deg2met, FGIF_p45, high_pop)
FGIF_p45i = [round(FGIF_p45(:,1),0) round(FGIF_p45(:,2),0)];
IM = img < 150;
if ismember(1,IM) == 0
    flag = 1; %no low pop zones
elseif ismember(0,IM) == 0
    flag = 2; %all low pop zones
else
    flag = 0;
end
% No using regionprops for finding crash pts
% Crossreference each FGIF pt with each blob
% Each FGIF pt within each blob has area associated with that blob 
% How to prioritize: - get number of high pop zones within blob
%                    - prioritize no high pop zones, then area, then dist

if flag == 0;
    s = regionprops(IM, 'Area', 'Perimeter', 'PixelList', 'Centroid');
    x_coord = [];
    y_coord = [];
    x_centroid = [];
    y_centroid = [];
    for i = 1:length(s)
        blob_area(i) = s(i).Area;
        blob_outline(i) = s(i).Perimeter;
        x_coord = [x_coord s(i).PixelList(:,1)'];
        y_coord = [y_coord s(i).PixelList(:,2)'];
    end
    
    % Find crash pts and distance 
    x_tmp = [];
    y_tmp = [];
    if size(FGIF_p45i,1) > 1
        for j = 1:size(FGIF_p45i,1)
            if img(FGIF_p45i(j,2), FGIF_p45i(j,1)) ~= high_pop
                x_tmp = [x_tmp; FGIF_p45i(j,1)];
                y_tmp = [y_tmp; FGIF_p45i(j,2)];
            end
        end
    end
    crash_loc = [x_tmp(:) y_tmp(:)];
    crash_loc = unique(crash_loc, 'rows');
    
    if isempty(crash_loc) == 0
        for k = 1:size(crash_loc,1)
            yrow2 = crash_loc(k,2) + min(imin) - 1;
            xcol2 = crash_loc(k,1) + min(jmin) - 1;
            ME22 = LS_USA(:,:,2); 
            ME32 = LS_USA(:,:,3);
            xnew12 = ME22(sub2ind(size(ME22), yrow2, xcol2));
            xnew2(k) = xnew12(1);
            ynew12 = ME32(sub2ind(size(ME32), yrow2, xcol2));
            ynew2(k) = ynew12(1);
            dist2(k) = sqrt((x_earth - xnew2(k))^2 + (y_earth - ynew2(k)).^2)*deg2met;
        end
        dist = dist2(:);
        x_centroid = xnew2(:);
        y_centroid = ynew2(:);       

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
        blob_area = zeros(size(x_centroid,1),1);
        crash_pts = [x_centroid(:) y_centroid(:) blob_area(:) dist(:)]; %high_low_rat not right
    else
        crash_pts = [];
        blob_coord = [];
        blob_numpixels = [];
    end
    
elseif flag == 2
    blob_coord = [];
    blob_numpixels = []; 
    s = regionprops(IM, 'Centroid', 'Area');
    x_centr = FGIF_p45i(end,1);
    y_centr = FGIF_p45i(end,2);
    
    yrow2 = round(y_centr,0) + min(imin) - 1   
    xcol2 = round(x_centr,0) + min(jmin) - 1
    ME22 = LS_USA(:,:,2); 
    ME32 = LS_USA(:,:,3);
    xnew12 = ME22(sub2ind(size(ME22), yrow2, xcol2));
    xnew = xnew12(1);
    ynew12 = ME32(sub2ind(size(ME32), yrow2, xcol2));
    ynew = ynew12(1);
    dist2 = sqrt((x_earth - xnew)^2 + (y_earth - ynew).^2)*deg2met;
    area = s.Area;
    crash_pts = [xnew ynew area dist2];
else
    blob_coord = [];
    blob_numpixels = [];
    crash_pts = [];
end