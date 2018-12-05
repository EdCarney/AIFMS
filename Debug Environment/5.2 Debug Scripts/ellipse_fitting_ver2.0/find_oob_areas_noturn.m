function [ oob ] = find_oob_areas_noturn( bin_image,oob_color )
%oob: get array of oob pixels from binary image
    sz_bin = size(bin_image);
    idx = 1;
    for ii=1:sz_bin(1)
        for jj=1:sz_bin(2)
            if bin_image(ii,jj) == oob_color
                oob(idx,:) = [ii jj];
                idx = idx + 1;
            end
        end
    end
end

