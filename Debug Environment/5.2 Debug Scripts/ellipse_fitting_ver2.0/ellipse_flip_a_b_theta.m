function [ a_new,b_new,ang_new ] = ellipse_flip_a_b_theta( a,b,ang )
%ellipse_flip_a_b_theta: adjust ellipse parameters if b>a
    a_new = b;
    b_new = a;
    ang_new = ang+pi/2;
end

