function [goodIntersect, checkPoint, interBranchInd, interWaypointInd] = checkBranchIntersect(branches, newSegment, referencePoint, minAngle)
%CHECKBRANCHINTERSECT Summary of this function goes here
%   Detailed explanation goes here

    goodIntersect = false;
    interBranchInd = 0;
    interWaypointInd = 0;
    checkPoint = false;
    
    for i = 1:length(branches)
        for j = 2:length(branches(i).waypoints)
            if newSegment.doIntersect(branches(i).segments(j-1))
                % If intersecting segments can be found, check
                % if the resulting trimmed segments meet the
                % minimum angle requirement
                [checkEast, checkNorth] = newSegment.getIntersection(branches(i).segments(j-1));
                checkPoint = Waypoint(checkNorth, checkEast, 0.0, 1, referencePoint);
                newSegment_Trimmed = LineSegment(checkPoint, newSegment.startPoint);
                existingSegment_Trimmed = LineSegment(branches(i).segments(j-1).startPoint, checkPoint);
                if angleCheck(newSegment_Trimmed, existingSegment_Trimmed, minAngle)
                    interBranchInd = i;
                    interWaypointInd = j;
                    goodIntersect = true;
                    
%                     plot([newSegment.startPoint.east, newSegment.stopPoint.east], [newSegment.startPoint.north, newSegment.stopPoint.north], 'r-', ...
%                         [branches(i).segments(j-1).startPoint.east, branches(i).segments(j-1).stopPoint.east], [branches(i).segments(j-1).startPoint.north, branches(i).segments(j-1).stopPoint.north], 'b-', ...
%                         [checkEast], [checkNorth], 'kx',...
%                         [newSegment.startPoint.east, branches(i).segments(j-1).startPoint.east], [newSegment.startPoint.north, branches(i).segments(j-1).startPoint.north], 'ko', ...
%                         [newSegment.stopPoint.east, branches(i).segments(j-1).stopPoint.east], [newSegment.stopPoint.north, branches(i).segments(j-1).stopPoint.north], 'k+')
                    
                    return
                end
            end
        end
    end
end
