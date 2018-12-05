function [ ep ] = ellipse_abhkt_to_endpts( a,b,h,k,t )
%ellipse_abhk_to_endpts - calculates ellipse endpoints from a,b,h,k,t
    ep(1,:) = [k+a*sin(t) h+a*cos(t)];
    ep(2,:) = [k-a*sin(t) h-a*cos(t)];
    ep(3,:) = [k+b*cos(t+pi/2) h+b*sin(t+pi/2)];
    ep(4,:) = [k-b*cos(t+pi/2) h-b*sin(t+pi/2)];
end

