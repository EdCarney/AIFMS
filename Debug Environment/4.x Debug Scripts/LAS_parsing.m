
clear all
close all
clc

% Import all LAS data for the four applicable regions.

LAS1 = lasdata('USGS_LPC_MD_VA_Sandy_NCR_2014_18SUJ330317_LAS_2015.las', 'loadall');
LAS2 = lasdata('USGS_LPC_MD_VA_Sandy_NCR_2014_18SUJ330318_LAS_2015.las', 'loadall');
LAS3 = lasdata('USGS_LPC_MD_VA_Sandy_NCR_2014_18SUJ331317_LAS_2015.las', 'loadall');
LAS4 = lasdata('USGS_LPC_MD_VA_Sandy_NCR_2014_18SUJ331318_LAS_2015.las', 'loadall');

% dLong = (maxLong - minLong)/length(unX);
% dLat = (maxLat - minLat)/length(unX);
% 
% longs = [minLong:dLong:maxLong]';
% lats = [minLat:dLat:maxLat]';

% Combine all four LIDAR region arrays

x = unique([LAS1.x;LAS2.x;LAS3.x;LAS4.x]);
y = unique([LAS1.y;LAS2.y;LAS3.y;LAS4.y]);
z = unique([LAS1.z;LAS2.z;LAS3.z;LAS4.z]);
c = z;

% Compute average of each 30x30 square meter section of the high-fidelity
% LAS data

avg = 0;
xMin = min(x);
yMin = min(y);

for j = 1:100
    for i = 1:100
        clear IxL IxG Ix
        avg = 0;
        coord = [];
        coord_new = [];
        IxL = find(x<=(xMin+i*30));
        IxG = find(x>(xMin+30*(i-1)));
        Ix = intersect(IxG,IxL);
        for m = 1:length(Ix)
            if y(Ix(m)) > (30*j+yMin) || y(Ix(m)) < (30*(j-1)+yMin)
                continue
            else
                coord(m) = z(Ix(m));
            end
        end
        coord_new = nonzeros(coord);
        avg_all(i,j) = sum(coord_new)/numel(coord_new);
        std_all(i,j) = std(coord_new);
    end
end

avg_all = avg_all';

% Use a 3D scatter plot for the points and set the colormap to show
% altitude variations.

dx = [min(x):30:max(x)]';
dy = [min(y):30:max(y)]';

figure(1)
scatter3(x,y,z,0.1,c,'.')
hold on
colormap('jet')
colorbar
xlabel('X Coordinate (m)')
ylabel('Y Coordinate (m)')
zlabel('Altitude (m)')

% Use a mesh plot to show the result of the 100x100 matrix of resultant
% average altitudes

figure(2)
mesh(dx,dy,avg_all)
xlabel('X Coordinate (m)')
ylabel('Y Coordinate (m)')
zlabel('Altitude (m)')

% Use reduced matrix to compute 

Rad = 6378160; % earth radius in N 0°
ye = 38.9897; % earth latitude, y
oa = 2*pi*Rad*cos(ye*pi/180); % circumference of the earth at the current location (latitude N 49)
deg2met = 0.1*oa/360; %fFull-time equivalent circuit for the length of one degree GPS

minLong = -76.9635009140166;
maxLong = -76.9281357760136;
maxLat = 38.9995885774388;
minLat = 38.9855019063119;

dLat = (maxLat-minLat)/100;
dLong = (maxLong-minLong)/100;

dLats = [minLat:dLat:maxLat-dLat]';
dLongs = [minLong:dLong:maxLong-dLong]';

figure(3)
mesh(dLongs,dLats,avg_all)
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
zlabel('Altitude (m)')
