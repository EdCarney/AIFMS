clc
clear all
close all

%% Load processed LandScan USA Data
prep_LandScan;

%% Initialize variables for LZ determination
    % Load load full data set
    load('DEM_processed.mat')
    load('TRI_processed.mat')
    load('LC_processed.mat')
  
    % Check which dataset has the largest minimum and smallest maximum and
    % constrain the master matrix to those bounds.
    
    minLonDEM = min(min(DEM_filt(:,:,2)));
    minLonTRI = min(min(TRI_filt(:,:,2)));
    minLonLC  = min(min(LC_filt(:,:,2)));
    minLonLS  = long_min;
    
    maxLonDEM = max(max(DEM_filt(:,:,2)));
    maxLonTRI = max(max(TRI_filt(:,:,2)));
    maxLonLC  = max(max(LC_filt(:,:,2)));
    maxLonLS  = long_max;
    
    minLatDEM = min(min(DEM_filt(:,:,3)));
    minLatTRI = min(min(TRI_filt(:,:,3)));
    minLatLC  = min(min(LC_filt(:,:,3)));
    minLatLS  = lat_min;
    
    maxLatDEM = max(max(DEM_filt(:,:,3)));
    maxLatTRI = max(max(TRI_filt(:,:,3)));
    maxLatLC  = max(max(LC_filt(:,:,3)));
    maxLatLS  = lat_max;
    
    minLon = max([minLonDEM,minLonTRI,minLonLC,minLonLS]);
    maxLon = min([maxLonDEM,maxLonTRI,maxLonLC,maxLonLS]);
    minLat = max([minLatDEM,minLatTRI,minLatLC,minLatLS]);
    maxLat = min([maxLatDEM,maxLatTRI,maxLatLC,maxLatLS]);
    
    % Find closest data point in the LS dataset for the absolute max and
    % min lat/lon (we need to shift all datasets to align with exact data
    % points for the N_filt matrix.
    
    [c minLonIndN] = min(abs(N_filt(1,:,2)' - minLon));
    [c maxLonIndN] = min(abs(N_filt(1,:,2)' - maxLon));
    [c minLatIndN] = min(abs(N_filt(:,1,3)  - minLat));
    [c maxLatIndN] = min(abs(N_filt(:,1,3)  - maxLat));
    
    % Shift the absolute min and max values to the closest LS data point;
    % need to account for the case where the closest LS point exceeds the
    % domain of the absolute values, and shift the index accordingly.
    
    if N_filt(1,minLonIndN,2) < minLon
        minLon = N_filt(1,minLonIndN - 1,2);
    else
        minLon = N_filt(1,minLonIndN,2);
    end
    
    if N_filt(1,maxLonIndN,2) > maxLon
        maxLon = N_filt(1,maxLonIndN + 1,2);
    else
        maxLon = N_filt(1,maxLonIndN,2);
    end
    
    if N_filt(minLatIndN,1,3) < minLat
        minLat = N_filt(minLatIndN + 1,1,3);
    else
        minLat = N_filt(minLatIndN,1,3);
    end
    
    if N_filt(maxLatIndN,1,3) > maxLat
        maxLat = N_filt(maxLatIndN - 1,1,3);
    else
        maxLat = N_filt(maxLatIndN,1,3);
    end
    
    % Get corresponding indices for the DEM, TRI, and LC matrices
    
    [c minDemLonInd] = min(abs(DEM_filt(1,:,2)' - (minLon - d_deg)));
    [c maxDemLonInd] = min(abs(DEM_filt(1,:,2)' - (maxLon + d_deg)));
    [c minDemLatInd] = min(abs(DEM_filt(:,1,3) - (minLat - d_deg)));
    [c maxDemLatInd] = min(abs(DEM_filt(:,1,3) - (maxLat + d_deg)));

    [c minTriLonInd] = min(abs(TRI_filt(1,:,2)' - (minLon - d_deg)));
    [c maxTriLonInd] = min(abs(TRI_filt(1,:,2)' - (maxLon + d_deg)));
    [c minTriLatInd] = min(abs(TRI_filt(:,1,3) - (minLat - d_deg)));
    [c maxTriLatInd] = min(abs(TRI_filt(:,1,3) - (maxLat + d_deg)));

    [c minLcLonInd] = min(abs(LC_filt(1,:,2)' - (minLon - d_deg)));
    [c maxLcLonInd] = min(abs(LC_filt(1,:,2)' - (maxLon + d_deg)));
    [c minLcLatInd] = min(abs(LC_filt(:,1,3) - (minLat - d_deg)));
    [c maxLcLatInd] = min(abs(LC_filt(:,1,3) - (maxLat + d_deg)));

    % Set the size of the latitude and longitude span to start and stop at
    % the determined min and max values and to be separated by the delta
    % determined in the LS data.
    
    LonArraySize = int32(ceil(((maxLon - minLon)/d_deg)));
    LatArraySize = int32(ceil(((maxLat - minLat)/d_deg)));
    LonArray = zeros(1,LonArraySize);
    LatArray = zeros(1,LatArraySize);
    for i = 1:LonArraySize
        LonArray(i) = maxLon - d_deg * (double(i) - 1);
    end
    for i = 1:LatArraySize
        LatArray(i) = minLat + d_deg * (double(i) - 1);
    end
    
    % Set Xq and Yq for the linear interpolation to align with the above
    % arrays.
    
    [Xq,Yq] = meshgrid(LonArray,LatArray);

    % Use interp2 to align DEM data subset with LS data
        X_dem = DEM_filt(:,:,2);
        Y_dem = DEM_filt(:,:,3);
        V_dem = DEM_filt(:,:,1);

        % Reorder Y and V matrices to align with N_FGIF
        Y_dem = flipud(Y_dem);
        V_dem = flipud(V_dem);
        Vq_dem = interp2(X_dem,Y_dem,V_dem,Xq,Yq,'linear');

    % Use interp2 to align TRI data subset with LS data
        X_tri = TRI_filt(:,:,2);
        Y_tri = TRI_filt(:,:,3);
        V_tri = TRI_filt(:,:,1);

        % Reorder Y and V matrices to align with N_FGIF
        Y_tri = flipud(Y_tri);
        V_tri = flipud(V_tri);
        Vq_tri = interp2(X_tri,Y_tri,V_tri,Xq,Yq,'linear');

    % Use interp2 to align LC data subset with LS data
        X_lc = LC_filt(:,:,2);
        Y_lc = LC_filt(:,:,3);
        V_lc = LC_filt(:,:,1);

        % Reorder Y and V matrices to align with N_FGIF
        Y_lc = flipud(Y_lc);
        V_lc = flipud(V_lc);
        Vq_lc = round(interp2(X_lc,Y_lc,V_lc,Xq,Yq,'linear'));

%% Initialize debug variables
d_deg = 8.3333e-04;
num_wps = 5;
home_llh = [0; 76.8903; 38.9582; 20.0000];
waypoints_w = [
   50.0000   76.9063   38.9701   20.0000
   50.0000   76.9260   38.9809   20.0000
   50.0000   76.9571   38.9839   20.0000
   50.0000   76.9851   38.9866   20.0000
    0.1000   77.0155   38.9969   20.0000
         0         0         0         0
         0         0         0         0
         0         0         0         0
         0         0         0         0
         0         0         0         0
         0         0         0         0
         0         0         0         0
         0         0         0         0
         0         0         0         0
         0         0         0         0
		 ];
x_earth = 76.9448;
y_earth = 38.9826;
Alt = 175.0000;
v = 19.3289;
phi0 = 0.0339;
psi0 = 0.1057;
alpha = 0.0160;
alpha_min = -0.0528;
m = 1.1280;
Wsurf = 0.3107;
Ar = 6.4000;
cx0 = 0.0315;
Cl0 = 0.2600;
Wspan = 1.4100;
cx_alpha = 0.3053;
xnew_old1 = 76.9094;
ynew_old1 = 38.9669;
LandScan_data = N_processed;
LandScan_data_filt = N_filt;

[flag, xnew, ynew, xnew_old, ynew_old, FM_x, FM_xmin, FM_y, FM_ymin, vel] = engineOutGliding_debug(LandScan_data, LandScan_data_filt, d_deg, num_wps, home_llh, waypoints_w, x_earth, y_earth, Alt, v, phi0, psi0, alpha, alpha_min, m, Wsurf, Ar, cx0, Cl0, Wspan,  cx_alpha, xnew_old1, ynew_old1, Vq_dem, Vq_tri, Vq_lc, Xq, Yq)
