trimmed_MM = MM(26:50,26:35,:);
trimmed_MM(trimmed_MM==10000) = 1;
trimmed_MM(trimmed_MM==100) = 0;

% WGS84 ellipsoid constants:
a = 6378137;
e = 8.1819190842622e-2;

[rowNum, colNum, ~] = size(trimmed_MM);

ecefMM = zeros(rowNum, colNum, 3);

for i = 1:rowNum
    for j = 1:colNum
        alt = 0;
        lat = deg2rad(trimmed_MM(i, j, 6));
        lon = deg2rad(trimmed_MM(i, j, 5));
        N = a ./ sqrt(1 - e^2 .* sin(lat).^2);

        % convert to ECEF
        x = (N+alt) .* cos(lat) .* cos(lon);
        y = (N+alt) .* cos(lat) .* sin(lon);
        z = ((1-e^2) .* N + alt) .* sin(lat);
   
        ecefMM(i,j,1) = x;
        ecefMM(i,j,2) = y;
        ecefMM(i,j,3) = z;
    end
end

[north,east,down] = lla2ned(10, 10, 100, 10.1, 9.9, 0)

%highFidelityMM = 