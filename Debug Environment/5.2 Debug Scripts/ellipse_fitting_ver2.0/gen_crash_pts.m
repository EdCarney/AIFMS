function crash_pts = gen_crash_pts(blob_coord, blob_numpix, blob_centroid, img)
IM = img < 100;
s = regionprops(IM, 'Centroid');
x_cent = [];
y_cent = [];
for i = 1:length(s)
    x_cent = [x_cent s(i).Centroid(1)];
    y_cent = [y_cent s(i).Centroid(2)];
end

for i = 1:size(blob_numpix, 1)
    area(i) = blob_numpix(i,1);
    row = blob_coord(i,:,1);
    col = blob_coord(i,:,2);
    row(row == 0) = [];
    col(col == 0) = [];
    if isempty(find(IM(row, col) == 0)) == 0
        idx_high_pop = find(IM(row, col) == 0);
        high_pop = length(idx_high_pop);
    end
    high_low_rat(i) = high_pop/(area(i) - high_pop);
end
% crash_pts = [x_cent(:) y_cent(:) area(:) high_low_rat(:)];

end