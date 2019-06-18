function [lat,lon,alt] = enu2lla(east, north, up, lat_ref, lon_ref, alt_ref)
%NED2LLA converts NED to LLA; assumes WGS84 ellipsoid model

    % convert to radians
    lat_ref = deg2rad(lat_ref);
    lon_ref = deg2rad(lon_ref);

    % define wgs84 properties
    wgs84_a = 6378137.0;
    wgs84_b = 6356752.3142;
    wgs84_ecc = 8.1819190842622e-2;
    wgs84_ecc_sqrd = 6.69437999014e-3;

    % convert NED to ECEF
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

    ecef = C' * [north,east,-up]';

    x = ecef(1);
    y = ecef(2);
    z = ecef(3);

    % determine radius of earth at reference LLA and convert to ECEF
    Rew = wgs84_a / sqrt(1.0 - wgs84_ecc_sqrd * sin(lat_ref)^2);

    % convert LLA reference to ECEF
    x0 = (Rew + alt_ref) * cos(lat_ref) * cos(lon_ref);
    y0 = (Rew + alt_ref) * cos(lat_ref) * sin(lon_ref);
    z0 = ((1 - wgs84_ecc_sqrd) * Rew + alt_ref) * sin(lat_ref);

    % define differences
    xdiff = x + x0;
    ydiff = y + y0;
    zdiff = z + z0;

    % convert differences into LLA
    lon = atan2(ydiff,xdiff);

    p = hypot(xdiff, ydiff);
    lat = atan2(zdiff, p * (1 - wgs84_ecc_sqrd));
    alt = (p / cos(lat)) - Rew;

    % iterate to get latitude and altitude
    err = 1;

    while (abs(err) > 1e-10)
        Rew = wgs84_a / sqrt(1.0 - wgs84_ecc_sqrd * sin(lat)^2);
        alt = (p / cos(lat)) - Rew;

        err = atan2(zdiff * (1 + wgs84_ecc_sqrd * Rew * sin(lat) / zdiff), p) - lat;

        lat = lat + err;
    end

    % convert back to degrees
    lat = rad2deg(lat);
    lon = rad2deg(lon);

end