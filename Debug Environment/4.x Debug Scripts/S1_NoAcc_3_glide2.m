% State 1: No acceleration, phi and theta changes 
clear all
tic
%% Initialize
psi0 = 0*pi/180; % initial yaw
theta0 = 5*pi/180; % initial pitch 
phi0 = 3*pi/180; % initial bank
phi_max = 45*pi/180; % max bank
phi_min = -phi_max;
alpha = -2.3*pi/180; % angle of attack
m = 1.128;
g = 9.8;
v = 20;
hi = 150;
rho = 1.225;
s = 0.310678;
Ar = 6.41;
Cd0 = 0.0315;
Cl0 = 0.26;
CL = Cl0 + 2*pi*alpha;
e = 0.8;
Cdi = CL^2/(pi*Ar*e);
Cd = Cd0 + Cdi;
k = 1/e;
W = m*g;
vs = 0.5*rho*v^3*s*Cd/W;
vturn = v*sqrt(sec(phi_max));
vs_turn = vs*(sec(phi_max))^(3/2); 
R = vturn^2/(g*tan(abs(phi_max)));
n = 50;
flag = 1;

x_earth = 76.94695;
y_earth = 38.99479;

d_deg = 0.008333333;
long_min = 76.8426; % x
long_max = 77.0426;
lat_min = 38.8869; % y
lat_max = 39.0869;
xs = floor((80 - long_max)/d_deg);
ys = floor((lat_min - 35)/d_deg);
xe = ceil((80 - long_min)/d_deg);
ye = ceil((lat_max - 35)/d_deg);
load('Fake_LandScan.mat');
[x,y] = size(N_new);
d_deg_new = (long_max - long_min)/x;
H1 = long_max:-d_deg_new:(long_min + d_deg_new);
V1 = lat_min:d_deg_new:(lat_max - d_deg_new);
V2 = V1';
N_new(:,:,2) = H1(ones(x,1),:);
N_new(:,:,3) = V2(:,ones(y,1));

%% Determination of FGIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag == 0 % for when FGIF is not used
    x0 = R*cos(pi/2);
    y0 = R*sin(pi/2);
    Dglide0 = hi*v/vs
    xd0 = Dglide0*sin(0);
    yd0 = Dglide0*cos(0);
    dx_rel = x0 + xd0;
    dy_rel = y0 + yd0;
    dx = (dx_rel*cos(psi0) + dy_rel*sin(psi0))/(cos(psi0)^2 + sin(psi0)^2);
    dy = (dy_rel*cos(psi0) - dx_rel*sin(psi0))/(cos(psi0)^2 + sin(psi0)^2);
    FGIF = [dx dy];
    
else % generate FGIF 
    psi(1) = psi0;
    dpsi(1) = psi(1) - psi0;
    Larc(1) = R*abs(dpsi(1));
    dh_turn(1) = Larc*vs_turn/vturn;
    x0(1) = R*cos(pi/2 - dpsi(1));
    y0(1) = R*sin(pi/2 - dpsi(1));
    Dglide(1) = (hi - dh_turn(1))*v/vs;
    xd(1) = Dglide(1)*sin(dpsi(1));
    yd(1) = Dglide(1)*cos(dpsi(1));
    dx_rel(1) = x0(1) + xd(1);
    dy_rel(1) = y0(1) + yd(1);
    
    if dh_turn(1) > hi
        dx(1) = nan;
        dy(1) = nan;
        dxmin(1) = nan;
        dymin(1) = nan;
    else 
        dx(1) = (dx_rel(1)*cos(psi0) + dy_rel(1)*sin(psi0))/(cos(psi0)^2 + sin(psi0)^2);
        dy(1) = (dy_rel(1)*cos(psi0) - dx_rel(1)*sin(psi0))/(cos(psi0)^2 + sin(psi0)^2);
        dxmin(1) = (x0(1)*cos(psi0) + y0(1)*sin(psi0))/(cos(psi0)^2 + sin(psi0)^2);
        dymin(1) = (y0(1)*cos(psi0) - x0(1)*sin(psi0))/(cos(psi0)^2 + sin(psi0)^2);
        dx_n(1,:) = linspace(dxmin(1), dx(1), n);
        dy_n(1,:) = linspace(dymin(1), dy(1), n);
    end 
    
    dpi = 0.01745329251; 
    theta_range = 0:dpi:2*pi;
    
    for i = 2:length(theta_range)
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
            dx(i) = nan;
            dy(i) = nan;
            dxmin(i) = nan;
            dymin(i) = nan;
        else 
            dx(i) = (dx_rel(i)*cos(psi0) + dy_rel(i)*sin(psi0))/(cos(psi0)^2 + sin(psi0)^2);
            dy(i) = (dy_rel(i)*cos(psi0) - dx_rel(i)*sin(psi0))/(cos(psi0)^2 + sin(psi0)^2);
            dxmin(i) = (x0(i)*cos(psi0) + y0(i)*sin(psi0))/(cos(psi0)^2 + sin(psi0)^2);
            dymin(i) = (y0(i)*cos(psi0) - x0(i)*sin(psi0))/(cos(psi0)^2 + sin(psi0)^2);
            dx_n(i,:) = linspace(dxmin(i), dx(i), n);
            dy_n(i,:) = linspace(dymin(i), dy(i), n);
        end    
    end
    dx_n( dx_n == 0) = [];
    dy_n( dy_n == 0) = [];
    
    FGIF = [dx_n(:) dy_n(:)];
    FGIF(any(FGIF==0,2),:) = [];
end 
% 
% 
Rad = 6378160; %Earth radius in N 0°
ye = 38.9897; %earth latitude, y
oa = 2*pi*Rad*cos(ye*pi/180); %circumference of the earth at the current location (latitude N 49)
deg2met = oa/360; %Full-time equivalent circuit for the length of one degree GPS
dx_deg = FGIF(:,1)/deg2met;
dy_deg = FGIF(:,2)/deg2met;
x_ned1 = x_earth + dx_deg; %deg
x_ned = round(x_ned1*10^8)/10^8;
y_ned1 = y_earth + dy_deg; %deg
y_ned = round(y_ned1*10^8)/10^8;
FGIF = [x_ned(:) y_ned(:)];
FGIF = unique(FGIF(:,1:2), 'rows'); % [x, y]

%% Safest Response %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('N035W080.mat');
% N = map;
% 
% Map = map(ys:ye, xs:xe);
% Z = 34;
% N_new = cell2mat(arrayfun(@(Map) repmat(Map,Z,Z),Map,'uniform',0))./Z;
% save('/Users/andrewpoissant/Documents/Research/UAS/LandScan/SubLandScan.mat','N_new');

if flag == 0
    xnew = x_ned;
    ynew = y_ned;
    
else
    [FGIF_r, FGIF_c] = size(FGIF);
    tol = 0.0001;
    tol1 = 0.0002;
    
    maxlat = max(FGIF(:,2));
    [imax, ~] = find(abs(N_new(:,:,3) - (maxlat + 30/deg2met)) <= tol1);
    minlat = min(FGIF(:,2));
    [imin, ~] = find(abs(N_new(:,:,3) - (minlat - 30/deg2met)) <= tol1);
    maxlon = max(FGIF(:,1));
    [~, jmin] = find(abs(N_new(:,:,2) - (maxlon + 30/deg2met)) <= tol1);
    minlon = min(FGIF(:,1));
    [~, jmax] = find(abs(N_new(:,:,2) - (minlon - 30/deg2met)) <= tol1);
    
    N_FGIF = N_new(min(imin):max(imax), min(jmin):max(jmax), :);
    
    for i = 1:FGIF_r
        [row_temp, col_temp] = find(abs(N_FGIF(:,:,2) - FGIF(i,1)) <= tol & abs(N_FGIF(:,:,3) - FGIF(i,2)) <= tol);
        if isempty(row_temp) | isempty(col_temp)
            continue
        else 
            row(i) = row_temp + min(imin) - 1;
            col(i) = col_temp + min(jmin) - 1;
        end
    end

    SRM = [row(:) col(:)];
    SRM = SRM(any(SRM,2),:);
    nRows = size(SRM, 1);
    keep = true(nRows, 1);
    for rowId = 2:nRows
        if any(all(repelem(SRM(rowId,:),rowId-1,1) == SRM(1:rowId-1,:),2));
            keep(rowId) = false;
        end
    end
    SRM = SRM(keep,:);
  
    % Add Safety Factor
    sf = 0.003;
    
%     for i = 1:length(SRM(:,1))
%         SRMrr(i) = SRM(i,1) - 1;
%         SRMra(i) = SRM(i,1) + 1;
%         SRMcr(i) = SRM(i,2) - 1;
%         SRMca(i) = SRM(i,2) + 1;
%         SRMsf_tempr = SRM(i,1);
%         SRMsf_tempc = SRM(i,2);
%         if N_new(SRMrr(i), SRM(i,2)) >= sf | N_new(SRMra(i), SRM(i,2)) >= sf | N_new(SRM(i,1), SRMcr(i)) >= sf | N_new(SRM(i,1), SRMca(i)) >= sf
%             continue
%         else 
%             SRMsfr(i) = SRMsf_tempr;
%             SRMsfc(i) = SRMsf_tempc;
%         end
%     end
   

    for i = 1:length(SRM(:,1))
        SRMsf_tempr = SRM(i,1);
        SRMsf_tempc = SRM(i,2);
        N_new(SRMsf_tempr,SRMsf_tempc,1);
        if N_new(SRMsf_tempr,SRMsf_tempc,1) >= sf
            continue
        else
            SRMsfr(i) = SRMsf_tempr;
            SRMsfc(i) = SRMsf_tempc;
        end
    end
   
    clear SRM
    
    SRM = [SRMsfr' SRMsfc'];
    SRM = SRM(any(SRM,2),:);
    
    idx = sub2ind(size(N_new), SRM(:,1), SRM(:,2));                          
    minval = min(N_new(idx));                                               
    minidx = find(N_new(:)==minval);                                                                         
    [rmin, cmin] = ind2sub(size(N_new),intersect(minidx,idx));    
    C = [rmin(:) cmin(:)];
    ME2 = N_new(:,:,2);
    ME3 = N_new(:,:,3);
    xnew_p = ME2(sub2ind(size(ME2), rmin, cmin));
    ynew_p = ME3(sub2ind(size(ME3), rmin, cmin));
    ME4 = [xnew_p(:) ynew_p(:)];
    xnew = xnew_p(end);
    ynew = ynew_p(end);
  
%     plotting
%     figure(3)
%     hold on
%     plot(ME4(:,1), ME4(:,2), 'r*')
%     plot(xnew, ynew, 'g*')
end

tic
%% Determine Casualty Expectation (CE)
if flag ~= 0
    tol = 0.0002;
    [row, col] = find(abs(N_new(:,:,2) - xnew) <= tol & abs(N_new(:,:,3) - ynew) <= tol)
    idx = sub2ind(size(N_new), row, col)       
    PD = min(N_new(idx)/(d_deg_new*deg2met));
    Hp = 1.6256; 
    PF = 0.0217; 
    L = 1.8288; 
    W = 3.048;
    B = 0.3048; 
    COF = 0.68;  
    m = (ynew - y_earth)/(xnew - x_earth);
    psi_CE = atan(m);
    dpsi_CE = pi/2 - (psi_CE - psi0);
    DS = v^2/(2*COF*9.8);
    Dglide_CE = Hp*v/vs;
    x_CE = Dglide_CE*sin(dpsi_CE);
    y_CE = Dglide_CE*cos(dpsi_CE);
    DG = sqrt(x_CE^2 + y_CE^2);
    A_L = (L + DG + DS + 2*B)*(W + 2*B);
    PK = 1; 
    S = 1; 
    CE = PF*A_L*PK*S*PD;
end
 
if flag == 0
    tol = 0.0002;
    [row, col] = find(abs(N_new(:,:,2) - xnew) <= tol & abs(N_new(:,:,3) - ynew) <= tol);
    idx = sub2ind(size(N_new), row, col);                          
    PD_old = max(N_new(idx)/(d_deg_new*deg2met));
    Hp = 1.6256; 
    dpsi_CE_old = 0;
    PF = 0.47; 
    L = 1.8288; 
    W = 3.048;
    B = 0.3048; 
    COF = 0.68;  
    DS = v^2/(2*COF*9.8);
    Dglide_CE = Hp*v/vs;
    x_CE_old = Dglide_CE*sin(dpsi_CE_old);
    y_CE_old = Dglide_CE*cos(dpsi_CE_old);
    DG = sqrt(x_CE_old^2 + y_CE_old^2);
    A_L = (L + DG + DS + 2*B)*(W + 2*B);
    PK = 1; 
    S = 1; 
    CE_old = PF*A_L*PK*S*PD_old;
end

toc