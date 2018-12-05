function [ pt_list, numpixels ] = generate_line_pixels( pt1, pt2 )
%generate_line: get all integer points on a line given 2 endpoints of form [y x]
%   Detailed explanation goes here
    if pt1(1,1)==pt2(1,1) && pt1(1,2)==pt2(1,2)
        pt_list = pt1;
        numpixels = 1;
    else
        dist = norm(pt2-pt1, 2);
        resolution = .2; %pixels
        numpts = dist/resolution;
        yout = linspace( pt1(1),pt2(1),numpts );
        xout = linspace( pt1(2),pt2(2),numpts );

        %pixelize the line (int)
        clear pt_list;
        idx_l = 1;
        xprev = -99999; yprev = -99999;
        for idx=1:numpts
            xcur = round( xout(idx) ); ycur = round( yout(idx) );
            if (xcur ~= xprev) || (ycur ~= yprev)
                pt_list(idx_l,:) = [ycur xcur];
                idx_l = idx_l + 1;
            end
            xprev = xcur; yprev = ycur;
        end
        numpixels = idx_l - 1;
    end
end