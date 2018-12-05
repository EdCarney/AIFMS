function [bordergroup, group, centroid, numpixels] = find_contiguous_areas( img )
% find_contiguous_areas - flood-fill algorithm to obtain contiguous blocks of pixels
maxnumblocks = 1000;
lo_pop = 0;
cyan = 25;
color_perim = 20;

A = img;
sz = size(A);
% figure;
% display_pixel_array(A);

%initialize
groupnum = 1;
found_pixel = 1;

%% find the first lo_pop pixel with a hi_pop pixel neighbor.
while (found_pixel && groupnum<maxnumblocks)
    found_pixel = 0;
    ii = 1;
    groupnum;
    while ii<sz(1)
        jj = 1;
        while jj<sz(2)
            if (A(ii,jj) <= lo_pop)
                if (A(ii,min(jj+1, sz(2))) > 0 || A(ii,max(jj-1, 1)) > 0 || A(max(ii-1, 1),jj) > 0 || A(min(ii+1, sz(1)),jj) > 0)
                    found_pixel = 1;
                    avoid = [ii jj];
                    ii = sz(1);
                    jj = sz(2);
                end
            end
            jj = jj + 1;
        end
        ii = ii + 1;
    end
    if (found_pixel == 0)
        break;
    end
    %initialize the first lo_pop pixel.
    i_init = avoid(1);
    j_init = avoid(2);
    pix_idx = 1; prev_row_idx = 1;
    group(groupnum, pix_idx,:) = [i_init j_init];
    A(i_init, j_init) = color_perim; %display only

    %% find all the lo_pop pixels connected to the above pixel. Save it to the 'group' array.
    for swp_i=1:2
        incr_y = [1 -1];
        i_lim = [1 sz(1)];
        i_lim2 = [sz(1) 1];
        for swp_j=1:2
            incr_x = [1 -1];
            j_lim = [1 sz(2)];
            j_lim2 = [sz(2) 1];
            for ii=i_lim(swp_i):incr_y(swp_i):i_lim2(swp_i)
                for jj=j_lim(swp_j):incr_x(swp_j):j_lim2(swp_j)
                    if A(ii, jj) <= lo_pop
                        if (A(min(sz(1), max(1, ii-incr_y(swp_i))), jj) == color_perim || A(min(sz(1), max(1, ii-incr_y(swp_i))),  max(1, (jj-1))) == color_perim || A(min(sz(1), max(1, ii-incr_y(swp_i))),  min(sz(2), (jj+1))) == color_perim || A(ii, min(sz(2), max(1, (jj-incr_x(swp_j))))) == color_perim)
%                         if (A(ii, min(sz(2), max(1, (jj-incr_x(swp_j))))) == color_perim || A(min(sz(1), max(1, ii-incr_y(swp_i))), jj) == color_perim || A(min(sz(1), max(1, ii-incr_y(swp_i))),  max(1, (jj-1))) == color_perim || A(min(sz(1), max(1, ii-incr_y(swp_i))),  min(sz(2), (jj+1))) == color_perim)
                            pix_idx = pix_idx + 1;
                            prev_row_idx = prev_row_idx + 1;
                            group(groupnum, pix_idx,:) = [ii jj];
                            A(ii, jj) = color_perim; %display only
                        end
                    end
                end
            end
        end
    end

    %% get the centroid of a "blob". Save it to the 'centroid' array.
    numpixels(groupnum,1) = pix_idx;
    centroid(groupnum,:) = [round(sum(group(groupnum, :, 1))/numpixels(groupnum,1)) round(sum(group(groupnum, :, 2))/numpixels(groupnum,1))];
    A(centroid(groupnum,1), centroid(groupnum,2)) = cyan;
    groupnum = groupnum + 1;
end

disp('done with finding blob...now get the perimeter...');
%% Get the borders of each contiguous pixel-block and save it to 'bordergroup' array.
szg = size(group);
bordergroup = zeros(szg(1), szg(2), 2);
for gp_idx=1:groupnum-1
    thisgroup = group(gp_idx,:,:);
    thisgroup = reshape(thisgroup,size(thisgroup,2),size(thisgroup,3));
    [bordergroup(gp_idx,:,:), numpixels_itn] = find_perimeter( thisgroup,1,1 );
    numpixels(gp_idx,2) = numpixels_itn(2);
end

% figure;
% image(A);
% set(gca,'YDir','normal');

end
