function [ elpse,centroid,ra,ang ] = generate_ellipse_pixels( endpt1, endpt2, rb, col, backgrnd_img )
%generate_line: get all integer points on an ellipse given 2 endpoints
%of form [y x], and the semi-minor axis
    lineseg = endpt1 - endpt2;
    ra = norm(lineseg, 2)/2;              %major axis
    rb;                                   %minor axis
    centroid = (endpt1 + endpt2)/2;       %centroid
    ang = atan2(lineseg(1), lineseg(2));  %tilt angle
    %display the ellipse
    [x,y]=ellipse(ra,rb,ang,centroid(2),centroid(1),col);
    %pixelize the ellipse    
    numpts = size(x,2);
    xprev = -99999; yprev = -99999;
    idx_e = 1;
    clear elpse;
    for idx=1:numpts
        xcur = round( x(idx) ); ycur = round( y(idx) );
        if (xcur ~= xprev) || (ycur ~= yprev)
            elpse(idx_e,:) = [ycur xcur];
            idx_e = idx_e + 1;
        end
        xprev = xcur; yprev = ycur;
    end
end

