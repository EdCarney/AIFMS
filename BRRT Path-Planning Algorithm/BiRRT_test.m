clear all
close all
clc

Indices = [2,13;2,14;2,15;3,13;3,14;3,15;4,13;4,14;4,15;5,13;5,14;5,15;6,13;6,14;6,15];

LonLat = [
    -76.965,38.954;
    -76.964,38.954;
    -76.963,38.954;
    -76.965,38.955;
    -76.964,38.955;
    -76.963,38.955;
    -76.965,38.956;
    -76.964,38.956;
    -76.963,38.956;
    -76.965,38.957;
    -76.964,38.957;
    -76.963,38.957;
    -76.965,38.958;
    -76.964,38.958;
    -76.963,38.958
    ];
  
DEM = [10.293;7.563;4.5178;12.418;12.765;7.5987;13.852;14.323;15.007;15.545;14.917;15.481;14.685;14.623;15.642];
   
TRI = [1.43;0.54163;2.1454;0.39515;1.3927;2.1613;0.45747;0.39242;0.49842;0.47663;0.38665;0.16039;0.3187;0.22491;0.31633];
  
LC = [7;7;7;7;7;7;3;7;7;4;7;7;7;7;7];

deg2met = 86563;

x_UAV = -76.94;
y_UAV = 38.983;

Strip = LandingStrip(Indices, LonLat, DEM, TRI, LC, x_UAV, y_UAV, deg2met);

currentPosition = [x_UAV, y_UAV];

badArea1 = [
    -76.955, 38.9615;
    -76.960, 38.9615;
    -76.960, 38.9665;
    -76.955, 38.9665
    ];

badArea2 = [
    -76.951, 38.970;
    -76.946, 38.970;
    -76.946, 38.975;
    -76.951, 38.975
    ];

badArea3 = [
    -76.962, 38.972;
    -76.957, 38.972;
    -76.957, 38.977;
    -76.962, 38.977
    ];

alt = 250;
delta = 250;
numIterations = 5;
minAngle = 110;
heading = 200.0;

tic()
[Waypoints, totalDistance] = BiRRT(Strip, currentPosition, alt, delta, ...
    numIterations, {badArea1 badArea2 badArea3}, minAngle, heading, 0)
toc()