function isValid = validPoint(newPoint, closestNode, badWaypoints)
%VALIDPOINT Determines whether a point is valid for the BiRRT algorithm.
%	Original Author:
%       Edward Carney
%       Systems Engineer, Orbit Logic Inc.
%       Graduate Student, Aerospace Engineering, University of Maryland
%   Creation Date:
%       20 Oct 2018
%   Last Modified By:
%       Edward Carney
%   Last Modification Date:
%       3 Nov 2018
%
%   Takes two Point objects (newPoint, closestNode) and a 1 x n
%   cell array (badWaypoints) where n is the number of bad areas and each
%   cell consists on an m x 1 vector of the bounding waypoints. The new
%   point and closest node, and the line segment formed between them, will
%   be checked to see if they intersect any of the bad areas.

    % Will set is valid to true and nly set to false if it violates a
    % condition
    isValid = true;
    for i = 1:length(badWaypoints)
        intersect = false;
        badWaypointsSet = badWaypoints{i};
        badArea = [];
        for j = 1:length(badWaypointsSet)
            badArea = [badArea; badWaypointsSet(j).east, badWaypointsSet(j).north];
        end
        % If bad area is more than two points, and is not the last area
        % (which will be the heading bad area and must NOT be a closed
        % region), will use inpoly and polyxpoly to check conditions.
        if length(badArea) > 2 && i ~= length(badWaypoints)
            badArea = [badArea; badWaypointsSet(1).east, badWaypointsSet(1).north];
            [in,on] = inpolygon(newPoint.east,newPoint.north,badArea(:,1),badArea(:,2));
            [eastIntercept, northIntercept] = polyxpoly(...
                [closestNode.east; newPoint.east],...
                [closestNode.north; newPoint.north],...
                badArea(:,1),...
                badArea(:,2)...
                );
        % Otherwise, use same method for line segment objects; will defualt
        % in, on, eastIntercept, northIntercept so that they don't
        % erroneously eliminate potential points.
        elseif i ~= length(badWaypoints)
            in = 0;
            on = 0;
            eastIntercept = 0;
            northIntercept = 0;
            badLineSegment = LineSegment(badWaypointsSet(1), badWaypointsSet(2));
            newLineSegment = LineSegment(newPoint, closestNode);
            if badLineSegment.doIntersect(newLineSegment)
                intersect = true;
            end
         % The only remaining condition is that the considered bad area is
         % the last one in the set and is therefore the one reltaing to the
         % heading open rectangle; will treat this as a set of three line
         % segments.
        else
            in = 0;
            on = 0;
            eastIntercept = 0;
            northIntercept = 0;
            newLineSegment = LineSegment(newPoint, closestNode);
            for k = 1:length(badArea) - 1
                badLineSegment = LineSegment(badWaypointsSet(k), badWaypointsSet(k+1));
                if badLineSegment.doIntersect(newLineSegment)
                    intersect = true;
                    break
                end
            end
        end
        if in || on || intersect || length(eastIntercept) > 1 || length(northIntercept) > 1
            isValid = false;
            break
        end
    end
end

