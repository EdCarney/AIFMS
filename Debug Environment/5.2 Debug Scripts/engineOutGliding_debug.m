function [flag, xnew, ynew, xnew_old, ynew_old, FM_x, FM_xmin, FM_y, FM_ymin, vel] = S1G_NoAcc(LandScan_data, LandScan_data_filt, d_deg, num_wps, home_llh, waypoints_w, x_earth, y_earth, Alt, v, phi0, psi0, alpha, alpha_min, m, Wsurf, Ar, cx0, Cl0, Wspan,  cx_alpha, xnew_old1, ynew_old1, DEM_filt, TRI_filt, LC_filt, Xq, Yq)
coder.extrinsic('cd')
coder.extrinsic('wind_effects')
coder.extrinsic('ellipse_fitting_main')
coder.extrinsic('pwd')
coder.extrinsic('strcat')

% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if psi0 == 0 % debugging
    psi0 = single(0.0001);
end
vel = double(v);
n = 30;
g = 9.8;
phi_max = 45*pi/180; % max roll
hi = Alt;
rho = 1.225;
CL = Cl0 + 2*pi*alpha_min;
W = m*g;
q_bar = 1/2*rho*v^2;
D = (cx0 + cx_alpha*alpha_min)*Wsurf*q_bar;
vs = D*v/W;
vturn = v*sqrt(sec(phi_max));
vs_turn = vs*(sec(phi_max))^(3/2); 
R = v^2/(g*tan(abs(phi_max)));

psi = zeros(1,361);
dpsi = zeros(1,361);
Larc = zeros(1,361);
dh_turn = zeros(1,361);
x0 = zeros(1,361);
y0 = zeros(1,361);
Dglide = zeros(1,361);
xd = zeros(1,361);
yd = zeros(1,361);
dx_rel = zeros(1,361);
dy_rel = zeros(1,361);
dx = zeros(1,361);
dy = zeros(1,361);
dxmin = zeros(1,361);
dymin = zeros(1,361);

%% Determination of FGIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rad = 6378160; %Earth radius in N 0°
ye = home_llh(3); %earth latitude, y
oa = 2*pi*Rad*cos(ye*pi/180); %circumference of the earth at the current location (latitude N 49)
deg2met = oa/360; %Full-time equivalent circuit for the length of one degree GPS
xnew_old = single(xnew_old1);
ynew_old = single(ynew_old1);

% Calculate FGIF with GIHM
psi(1) = psi0;
dpsi(1) = psi(1) - psi0;
Larc(1) = R*abs(dpsi(1));
dh_turn(1) = Larc(1)*vs_turn/vturn;
x0(1) = R*cos(pi/2 - dpsi(1));
y0(1) = R*sin(pi/2 - dpsi(1));
Dglide(1) = (hi - dh_turn(1))*v/vs;
xd(1) = Dglide(1)*sin(dpsi(1));
yd(1) = Dglide(1)*cos(dpsi(1));
dx_rel(1) = x0(1) + xd(1);
dy_rel(1) = y0(1) + yd(1);
if dh_turn(1) > hi
    dx(1) = 0;
    dy(1) = 0;
    dxmin(1) = 0;
    dymin(1) = 0;
else 
    dx(1) = (dx_rel(1)*cos(pi/2 - psi0) + dy_rel(1)*sin(pi/2 - psi0))/(cos(pi/2 - psi0)^2 + sin(pi/2 - psi0)^2);
    dy(1) = (dy_rel(1)*cos(pi/2 - psi0) - dx_rel(1)*sin(pi/2 - psi0))/(cos(pi/2 - psi0)^2 + sin(pi/2 - psi0)^2);
    dxmin(1) = (x0(1)*cos(pi/2 - psi0) + y0(1)*sin(pi/2 - psi0))/(cos(pi/2 - psi0)^2 + sin(pi/2 - psi0)^2);
    dymin(1) = (y0(1)*cos(pi/2 - psi0) - x0(1)*sin(pi/2 - psi0))/(cos(pi/2 - psi0)^2 + sin(pi/2 - psi0)^2);
end 

dpi = 0.01745329251; 
theta_range = 0:dpi:2*pi;

for i = 2:length(theta_range);
    psi(i) = psi(i-1) + dpi;
    dpsi(i) = psi(i) - psi0;
    while dpsi(i) > pi
       dpsi(i) = dpsi(i) - 2*pi;
    end
    Larc(i) = R*abs(dpsi(i));
    dh_turn(i) = Larc(i)*vs_turn/vturn;    
    x0(i) = R*cos(pi/2 - dpsi(i));
    y0(i) = R*sin(pi/2 - dpsi(i));
    Dglide(i) = (hi - dh_turn(i))*v/vs;
    xd(i) = Dglide(i)*sin(dpsi(i));
    yd(i) = Dglide(i)*cos(dpsi(i));
    dx_rel(i) = x0(i) + xd(i);
    dy_rel(i) = y0(i) + yd(i);
    if dh_turn(i) > hi
        dx(i) = 0;
        dy(i) = 0;
        dxmin(i) = 0;
        dymin(i) = 0;
    else 
        dx(i) = (dx_rel(i)*cos(pi/2 - psi0) + dy_rel(i)*sin(pi/2 - psi0))/(cos(pi/2 - psi0)^2 + sin(pi/2 - psi0)^2);
        dy(i) = (dy_rel(i)*cos(pi/2 - psi0) - dx_rel(i)*sin(pi/2 - psi0))/(cos(pi/2 - psi0)^2 + sin(pi/2 - psi0)^2);
        dxmin(i) = (x0(i)*cos(pi/2 - psi0) + y0(i)*sin(pi/2 - psi0))/(cos(pi/2 - psi0)^2 + sin(pi/2 - psi0)^2);
        dymin(i) = (y0(i)*cos(pi/2 - psi0) - x0(i)*sin(pi/2 - psi0))/(cos(pi/2 - psi0)^2 + sin(pi/2 - psi0)^2);
    end  
end

FM_x = dx + x_earth*deg2met;
FM_xmin = dxmin + x_earth*deg2met;
FM_y = dy + y_earth*deg2met;
FM_ymin = dymin + y_earth*deg2met;

% Populate FGIF
dn_x = zeros(361, n); 
dx_n = zeros(360, n);
dn_y = zeros(361, n);
dy_n = zeros(360, n);
for row1 = 1:numel(dx)
    dn_x(row1,:) = linspace(dxmin(row1), dx(row1), n);
end
dx_n = dn_x;

for row2 = 1:numel(dy)
    dn_y(row2,:) = linspace(dymin(row1), dy(row2), n);
end
dy_n = dn_y;

figure(2)
plot(dx,dy)
% Create and refine FGIF
FGIF = [0 0];
FGIF = [dx_n(:) dy_n(:)]; 
FGIF(any(FGIF==0,2),:) = [];

dx_deg = FGIF(:,1)/deg2met;
dy_deg = FGIF(:,2)/deg2met;
x_ned1 = x_earth + dx_deg; %deg
x_ned = round(x_ned1*10^8)/10^8;
y_ned1 = y_earth + dy_deg; %deg
y_ned = round(y_ned1*10^8)/10^8;

FGIF = [x_ned(:) y_ned(:)];
FGIF = unique(FGIF(:,1:2), 'rows'); % [x, y]

%% Safest Response %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0_deg = x0(:)/deg2met;
y0_deg = y0(:)/deg2met;
x_ned01 = x_earth + x0_deg;
x_ned0 = round(x_ned01*10^8)/10^8;
y_ned01 = y_earth + y0_deg;
y_ned0 = round(y_ned01*10^8)/10^8;
xhome = home_llh(2);
yhome = home_llh(3);
xfinwpt = waypoints_w(num_wps, 2);
yfinwpt = waypoints_w(num_wps, 3);
[in1, on1] = inpolygon(xhome, yhome, x_ned, y_ned);
[in2, on2] = inpolygon(xhome, yhome, x_ned0, y_ned0);
[in3, on3] = inpolygon(xfinwpt, yfinwpt, x_ned, y_ned);
[in4, on4] = inpolygon(xfinwpt, yfinwpt, x_ned0, y_ned0);
test1 = numel(xhome(in1));
test2 = numel(xhome(in2));
test3 = numel(xfinwpt(in3));
test4 = numel(xfinwpt(in4));

if  test3 == 1 & test4 == 0 % if end wp is within FGIF
    xnew = waypoints_w(num_wps,2);
    ynew = waypoints_w(num_wps,3);
    
elseif test1 == 1 & test2 == 0 % if home wp is within FGIF
    xnew = home_llh(2);
    ynew = home_llh(3);

else % home or end wp not in FIGF
    [FGIF_r, FGIF_c] = size(FGIF);
    row = zeros(1, 150000);
    col = zeros(1, 150000);
    tol = 0.00024;
    tol1 = 0.0005;
    
    maxlat = max(FGIF(:,2));
    [imax, ~] = find(abs(LandScan_data(:,:,3) - (maxlat + 30/deg2met)) <= tol1);
    minlat = min(FGIF(:,2));
    [imin, ~] = find(abs(LandScan_data(:,:,3) - (minlat - 30/deg2met)) <= tol1);
    maxlon = max(FGIF(:,1));
    [~, jmin] = find(abs(LandScan_data(:,:,2) - (maxlon + 30/deg2met)) <= tol1);
    minlon = min(FGIF(:,1));
    [~, jmax] = find(abs(LandScan_data(:,:,2) - (minlon - 30/deg2met)) <= tol1);
 
    N_FGIF = LandScan_data(min(imin):max(imax), min(jmin):max(jmax), :);
    N_FGIF_filt = LandScan_data_filt(min(imin):max(imax), min(jmin):max(jmax), :);
    
    
    %% LZ Determination

    % Set minimum and maximum latitude and longitude values based on N_FGIF
    
    minLon = min(min(N_FGIF_filt(:,:,2)));
    maxLon = max(max(N_FGIF_filt(:,:,2)));
    minLat = min(min(N_FGIF_filt(:,:,3)));
    maxLat = max(max(N_FGIF_filt(:,:,3)));
    
    % Get corresponding indices for the DEM and TRI matrices
    
    c = 0;
    minLonInd = 0;
    maxLonInd = 0;
    minLatInd = 0;
    maxLatInd = 0;
    
    [c minLonInd] = min(abs(Xq(1,:)' - minLon));
    [c maxLonInd] = min(abs(Xq(1,:)' - maxLon));
    [c minLatInd] = min(abs(Yq(:,1)  - minLat));
    [c maxLatInd] = min(abs(Yq(:,1)  - maxLat));
    
    % Set subsets of the DEM, TRI, LC, Xq, and Yq data to align with the
    % FGIF
    
    DEM_FGIF_filt = DEM_filt(minLatInd:maxLatInd,maxLonInd:minLonInd);
    TRI_FGIF_filt = TRI_filt(minLatInd:maxLatInd,maxLonInd:minLonInd);
    LC_FGIF_filt = LC_filt(minLatInd:maxLatInd,maxLonInd:minLonInd);
    Xq_FGIF = Xq(minLatInd:maxLatInd,maxLonInd:minLonInd);
    Yq_FGIF = Yq(minLatInd:maxLatInd,maxLonInd:minLonInd);
    
    % Debug plots to check that chosen LZ makes sense. Can comment out if
    % necessary
    
    figure(5)
    hold on
    surf(-Xq_FGIF,Yq_FGIF,DEM_FGIF_filt);
    title('DEM Data - FGIF')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    h = colorbar;
    ylabel(h, 'Elevation (m)')
    hold off
    
    figure(6)
    hold on
    surf(-Xq_FGIF,Yq_FGIF,TRI_FGIF_filt);
    title('TRI Data - FGIF')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    h = colorbar;
    ylabel(h, 'Rugedness (m)')
    hold off
    
    figure(7)
    hold on
    surf(-Xq_FGIF,Yq_FGIF,LC_FGIF_filt);
    title('LC Data - FGIF')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    h = colorbar;
    ylabel(h, 'Land Cover Type')
    hold off
    
    % Check to make sure that FGIF dimensions do not exceed the domain of
    % the available TRI, DEM, or LC data; currently if it does the LZ
    % determination will be skipped. Could theoretically cut down the FGIF
    % to align with the above datasets in the future.
    
    [m_N,n_N,p_N] = size(N_FGIF_filt);
    [m_check,n_check] = size(DEM_FGIF_filt);
        
    if m_check ~= m_N || n_check ~= n_N
        disp('WARN: FGIF exceeds domain of of LC, TRI, or DEM data.');
        disp('INFO: Trimming FGIF to fit available data sets.');
        
        [c minLonIndNew] = min(abs(N_FGIF_filt(1,:,2)' - min(min(Xq_FGIF))));
        [c maxLonIndNew] = min(abs(N_FGIF_filt(1,:,2)' - max(max(Xq_FGIF))));
        [c minLatIndNew] = min(abs(N_FGIF_filt(:,1,3)  - min(min(Yq_FGIF))));
        [c maxLatIndNew] = min(abs(N_FGIF_filt(:,1,3)  - max(max(Yq_FGIF))));
        N_FGIF_filt = N_FGIF_filt(minLatIndNew:maxLatIndNew,maxLonIndNew:minLonIndNew,:);
    end
    
    figure(8)
    hold on
    surf(-Xq_FGIF,Yq_FGIF,N_FGIF_filt(:,:,1));
    title('LS Data - FGIF')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    hold off
    
    master_matrix = cat(6,N_FGIF_filt(:,:,1),DEM_FGIF_filt,TRI_FGIF_filt,LC_FGIF_filt,Xq_FGIF,Yq_FGIF);
    bestStrip = LZ_Calculation_FGIF(master_matrix, x_earth, y_earth, deg2met, Alt);
    
    currentPosition = [x_earth, y_earth];
    alt = 250;
    delta = 100;
    numIterations = 10;
    minAngle = 145;
    heading = 180;
    
    % badAreas = getBadAreas(master_matrix);
    
    [Waypoints, totalDistance] = BiRRT(bestStrip, currentPosition, alt, delta, numIterations, {}, minAngle, heading)

    %% New Ellipse Fitting Algorithm %% 
    oldFolder = pwd;
    newFolder = strcat(pwd, '\ellipse_fitting_ver2.0');
    cd(newFolder)
    LandScan_USA = LandScan_data(:,:,:);
    img = N_FGIF(:,:,1);
    img_filt = N_FGIF_filt(:,:,1);
    img = N_FGIF(:,:,1);
    img_filt = N_FGIF_filt(:,:,1);
    img_45m = [];
    img_filt_45m = [];
    LandScan_USA_45m = [];
    imin_45m = [];
    jmin_45m = [];
    FM = 0;
    [xnewd ynewd] = ellipse_fitting_main(img, img_filt, img_45m, img_filt_45m, LandScan_USA, LandScan_USA_45m, x_earth, y_earth, imin, jmin, imin_45m, jmin_45m, deg2met, FM, FGIF, d_deg, minlat, minlon, psi0); %need to check crash point algorithm
    xnew = single(xnewd)
    ynew = single(ynewd)
    cd(oldFolder)
end

flag = 1;   