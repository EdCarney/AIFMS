%% Prepare LandScan USA Data %%

load('C:\Users\ecarney\Documents\School\AIFMS Research\UAV Sim ver 5.2\UAV Sim ver 5.2\ellipse_fitting_ver2.0\map_processed_filtered_7_8_50.mat')
load('C:\Users\ecarney\Documents\School\AIFMS Research\UAV Sim ver 5.2\UAV Sim ver 5.2\ellipse_fitting_ver2.0\map_processed_7_8_50.mat')
load('C:\Users\ecarney\Documents\School\AIFMS Research\UAV Sim ver 5.2\UAV Sim ver 5.2\ellipse_fitting_ver2.0\map_processed_filtered_45m_7_8_50.mat')
load('C:\Users\ecarney\Documents\School\AIFMS Research\UAV Sim ver 5.2\UAV Sim ver 5.2\ellipse_fitting_ver2.0\map_processed_45m_7_8_50.mat')

%Define Mission Envelope
A_filt = map_processed_filtered;
A_filt_45m = map_processed_filtered_45m;
A_processed = map_processed;
A_processed_45m = map_processed_45m;
N_filt = zeros(size(A_filt,1), size(A_filt,2));
N_filt_45m = zeros(size(A_filt_45m,1), size(A_filt_45m,2));
N_processed = cat(3, N_filt, N_filt, N_filt);
N_processed_45m = cat(3, N_filt_45m, N_filt_45m, N_filt_45m);
N_filt = cat(3, N_filt, N_filt, N_filt);
N_filt_45m = cat(3, N_filt_45m, N_filt_45m, N_filt_45m);
N_filt(:,:,1) = A_filt*100; 
N_filt_45m(:,:,1) = A_filt_45m*100;
N_processed(:,:,1) = A_processed*100;
N_processed_45m(:,:,1) = A_processed_45m*100;
d_deg = cellsize;
d_deg_45m = cellsize_45m;
[x,y] = size(N_filt(:,:,1));
[x45,y45] = size(N_filt_45m(:,:,1));

long_max = abs(llcorner_this_tile(2));
long_min = long_max - size(A_filt,2)*d_deg;
lat_max = llcorner_this_tile(1);
lat_min = lat_max - size(A_filt,1)*d_deg;
H1 = long_max:-d_deg:(long_min + d_deg);
H1_45m = long_max:-d_deg_45m:(long_min + d_deg_45m);
V1 = lat_min:d_deg:(lat_max - d_deg);
V1_45m = lat_min:d_deg_45m:(lat_max - d_deg_45m);
V2 = V1';
V2_45m = V1_45m';
N_filt(:,:,2) = repmat(H1, x, 1);
N_filt_45m(:,:,2) = repmat(H1_45m, x45,1);
N_filt(:,:,3) = repmat(V2, 1, y); % 100 = abv thresh, need to swap for regionprops 
N_filt_45m(:,:,3) = repmat(V2_45m, 1, y45);
N_filt(N_filt(:,:,1) == 100) = 5;
N_filt_45m(N_filt_45m(:,:,1) == 100) = 5;
N_filt(N_filt(:,:,1) == 0) = 100;
N_filt_45m(N_filt_45m(:,:,1) == 0) = 100;
N_filt(N_filt(:,:,1) == 5) = 0; % now 100 = below thresh
N_filt_45m(N_filt_45m(:,:,1) == 5) = 0;

N_processed(:,:,2) = repmat(H1, x, 1);
N_processed_45m(:,:,2) = repmat(H1_45m, x45,1);
N_processed(:,:,3) = repmat(V2, 1, y); % 100 = abv thresh, need to swap for regionprops 
N_processed_45m(:,:,3) = repmat(V2_45m, 1, y45);
N_processed(N_processed(:,:,1) == 100) = 5;
N_processed_45m(N_processed_45m(:,:,1) == 100) = 5;
N_processed(N_processed(:,:,1) == 0) = 100;
N_processed_45m(N_processed_45m(:,:,1) == 0) = 100;
N_processed(N_processed(:,:,1) == 5) = 0; % now 100 = below thresh
N_processed_45m(N_processed_45m(:,:,1) == 5) = 0; % now 100 = below thresh