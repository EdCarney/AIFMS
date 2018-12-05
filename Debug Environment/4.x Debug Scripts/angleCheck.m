function isValid = angleCheck(segment1, segment2, minAngle)
%ANGLECHECK Function to determine if the angle subtended by two line
%segments meets a minimum angle requirement.
%	Original Author:
%       Edward Carney
%       Systems Engineer, Orbit Logic Inc.
%       Graduate Student, Aerospace Engineering, University of Maryland
%   Creation Date:
%       20 Oct 2018
%   Last Modified By:
%       Edward Carney
%   Last Modification Date:
%       1 Nov 2018
%
%   Takes two LineSegment objects (segment1, segment2) and a minimum angle
%   in degrees (minAngle). Outputs a boolean 1 is the segments meet the
%   requirement and a 0 otherwise.

    isValid = true;
    u = [(segment1.stopPoint.east - segment1.startPoint.east),...
        (segment1.stopPoint.north - segment1.startPoint.north),...
        0.0];
    v = [(segment2.stopPoint.east - segment2.startPoint.east),...
        (segment2.stopPoint.north - segment2.startPoint.north),...
        0.0];
    ThetaInDegrees = atan2d(norm(cross(u,v)),dot(u,v));
    if ThetaInDegrees > (180 - minAngle) && ThetaInDegrees < (180 + minAngle)
        isValid = false;
    end
end

