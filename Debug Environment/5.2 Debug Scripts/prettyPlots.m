clear all
close all
clc

%%
% Script to generate pretty plots for the SciTech paper.

%% Load processed LandScan USA Data
prep_LandScan;

%% Load load full datasets

load('DEM_processed.mat')
load('TRI_processed.mat')
load('LC_processed.mat')

%% Get data sets in same resolution
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
        
        
%% Get LandScan data to plot
N_plot = N_filt;
N_plot(N_plot(:,:,1) == 100) = 0;;
N_plot(N_plot(:,:,1) > 0) = 1;

%% Plot full and equivalent resolution datasets

    figure(7)
    hold on
    title('LS Data - Full Resolution')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    zlabel('Population')
    mesh(-N_plot(:,:,2),N_plot(:,:,3),N_plot(:,:,1))
    h = colorbar;
    ylabel(h, 'Population')
    hold off    


    figure(1)
    hold on
    title('DEM Data - Full Resolution')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    zlabel('Elevation (m)')
    mesh(-DEM_filt(:,:,2),DEM_filt(:,:,3),DEM_filt(:,:,1))
    h = colorbar;
    ylabel(h, 'Elevation (m)')
    hold off

    figure(2)
    hold on
    title('TRI Data - Full Resolution')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    zlabel('Ruggedness (m)')
    mesh(-TRI_filt(:,:,2),TRI_filt(:,:,3),TRI_filt(:,:,1))
    h = colorbar;
    ylabel(h, 'Ruggedness (m)')
    hold off

    figure(3)
    hold on
    title('LC Data - Full Resolution')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    zlabel('Land Cover Type')
    mesh(-LC_filt(:,:,2),LC_filt(:,:,3),LC_filt(:,:,1))
    h = colorbar;
    ylabel(h, 'Land Cover Type')
    hold off
        
    figure(4)
    hold on
    title('DEM Data - Equivalent Resolution')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    zlabel('Elevation (m)')
    mesh(-Xq,Yq,Vq_dem)
    h = colorbar;
    ylabel(h, 'Elevation (m)')
    hold off

    figure(5)
    hold on
    title('TRI Data - Equivalent Resolution')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    zlabel('Ruggedness (m)')
    mesh(-Xq,Yq,Vq_tri)
    h = colorbar;
    ylabel(h, 'Ruggedness (m)')
    hold off

    figure(6)
    hold on
    title('LC Data - Equivalent Resolution')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    zlabel('Land Cover Type')
    mesh(-Xq,Yq,Vq_lc)
    h = colorbar;
    ylabel(h, 'Land Cover Type')
    hold off
