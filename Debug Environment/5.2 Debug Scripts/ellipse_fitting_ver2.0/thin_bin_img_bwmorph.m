function A_thin = thin_bin_img_bwmorph(A)

A_thin = bwmorph(A, 'thin');
% need to replace specific nodes with backgrnd_col
end