function [east,north,up] = lla2enu(lat, lon, alt, lat_ref, lon_ref, alt_ref)
%LLA2NED converts LLA to ENU; assumes WGS84 ellipsoid model

% convert to radians
lat = deg2rad(lat);
lon = deg2rad(lon);
lat_ref = deg2rad(lat_ref);
lon_ref = deg2rad(lon_ref);

% define wgs84 properties
wgs84_a = 6378137.0;
wgs84_b = 6356752.3142;
wgs84_ecc = 8.1819190842622e-2;
wgs84_ecc_sqrd = 6.69437999014e-3;

% determine radius of earth at desired LLA and convert to ECEF
Rew = wgs84_a / sqrt(1.0 - wgs84_ecc_sqrd * sin(lat)^2);
Rns = wgs84_a * (1.0 - wgs84_ecc_sqrd) / (1.0 - wgs84_ecc_sqrd * sin(lat)^2)^1.5;

x = (Rew + alt) * cos(lat) * cos(lon);
y = (Rew + alt) * cos(lat) * sin(lon);
z = ((1 - wgs84_ecc_sqrd) * Rew + alt) * sin(lat);

% determine radius of earth at reference LLA and convert to ECEF
Rew = wgs84_a / sqrt(1.0 - wgs84_ecc_sqrd * sin(lat_ref)^2);
Rns = wgs84_a * (1.0 - wgs84_ecc_sqrd) / (1.0 - wgs84_ecc_sqrd * sin(lat_ref)^2)^1.5;

x0 = (Rew + alt_ref) * cos(lat_ref) * cos(lon_ref);
y0 = (Rew + alt_ref) * cos(lat_ref) * sin(lon_ref);
z0 = ((1 - wgs84_ecc_sqrd) * Rew + alt_ref) * sin(lat_ref);

% define differences
xdiff = x - x0;
ydiff = y - y0;
zdiff = z - z0;

% convert ECEF difference to NED
C = zeros(3);

C(1,1) = -sin(lat_ref) * cos(lon_ref);
C(1,2) = -sin(lat_ref) * sin(lon_ref);
C(1,3) = cos(lat_ref);

C(2,1) = -sin(lon_ref);
C(2,2) = cos(lon_ref);
C(2,3) = 0;

C(3,1) = -cos(lat_ref) * cos(lon_ref);
C(3,2) = -cos(lat_ref) * sin(lon_ref);
C(3,3) = -sin(lat_ref);

ned = C * [xdiff,ydiff,zdiff]';
north = ned(1);
east = ned(2);
up = -ned(3);

end
