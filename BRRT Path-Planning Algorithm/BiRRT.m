function [Waypoints, totalDistance] = BiRRT(LandingStrip, currentPosition, alt, delta, numIterations, badAreas, minAngle, heading)
%BIRRT Builds a path between the LzStrip and UAV using a Bi-directional
%Rapidly Exploring Random Tree algorithm.
%	Original Author:
%       Edward Carney
%       Systems Engineer, Orbit Logic Inc.
%       Graduate Student, Aerospace Engineering, University of Maryland
%   Creation Date:
%       20 Oct 2018
%   Last Modified By:
%       Edward Carney
%   Last Modification Date:
%       28 Nov 2018
%
%   Takes seven arguments:
%       LandingStrip: the desired landing strip input as a LandingStrip
%         object.
%       currentPosition: the current UAV position as a 1 x 2 vector of
%         [Lon Lat].
%       alt: the current UAV altitude as a double.
%       delta: the delta value (in meters) for the growth of the tree as a
%         double.
%       numIterations: the number of iterations as an integer (more
%         iterations yield a better solution at the cost of a longer
%         runtime).
%       badAreas: areas that the final path must avoid provided in a 1 x n
%         cell array where n is the number of bad areas and each cell
%         consists of an m x 2 matrix of the bounding longitude and latitude
%         coordinates for the bad area where m is the number of bounding
%         points.
%       minAngle: the minimum angle as a double that must be subtended by
%         two line segments for a new branch to be valid; must be between
%         0 and 180 (larger angles yield 'smoother' paths at the expense of
%         longer runtimes).
%       heading: the current heading of the UAV assumed to be measured in
%         degrees clockwise from north; must be between 0 and 360.
%   Outputs two parameters:
%       Waypoints: an n x 1 vector of Waypoint objects, where n is the
%         number of waypoints between the UAV position and LZ strip
%         centriod.
%       totalDistance: the total distance of the path in meters as a
%         double.



    % Will define the centriod of the LzStrip as the reference point for all
    % conversions to NED; defined as Lat, Lon, Alt.
    LandingStrip.LonLat = [LandingStrip.LonLat(:,1) * -1, LandingStrip.LonLat(:,2)];
    LandingStrip.centroid = [LandingStrip.centroid(1) * -1, LandingStrip.centroid(2)];
    currentPosition = [currentPosition(1) * -1, currentPosition(2)];
    referencePoint = [LandingStrip.centroid(2),LandingStrip.centroid(1), LandingStrip.avgDEM];

    % Define start and stop points for the algorithm.
    startPoint = Waypoint(currentPosition(2), currentPosition(1), alt, 0, referencePoint);
    stopPoint  = Waypoint(LandingStrip.centroid(2), LandingStrip.centroid(1), LandingStrip.avgDEM, 0, referencePoint);

    % Set domain of algorithm freespace as [minE maxE ; minN maxN].
    minE = Waypoint(referencePoint(1), min(min(LandingStrip.LonLat(:,1), startPoint.lon)), 0.0, 0, referencePoint).east;
    maxE = Waypoint(referencePoint(1), max(max(LandingStrip.LonLat(:,1), startPoint.lon)), 0.0, 0, referencePoint).east;
    minN = Waypoint(min(min(LandingStrip.LonLat(:,2), startPoint.lat)), referencePoint(2), 0.0, 0, referencePoint).north;
    maxN = Waypoint(max(max(LandingStrip.LonLat(:,2), startPoint.lat)), referencePoint(2), 0.0, 0, referencePoint).north;
    % Determine an approximate minimum possible number of line segments
    % required to execute a 90 degree turn with the provided minAngle;
    % will use this to expand the freespace limits to allow UAV to execute
    % a 180 degree turn from current heading if necessary. Note here that
    % the startPoint is the current UAV position.
    angleDiffFromStraight = 180 - minAngle;
    minNumDeltas = ceil(90 / angleDiffFromStraight) + 4;
    if startPoint.east - minNumDeltas * delta < minE
        minFreespaceE = startPoint.east - minNumDeltas * delta;
    else
        minFreespaceE = minE;
    end
    
    if startPoint.east + minNumDeltas * delta > maxE
        maxFreespaceE = startPoint.east + minNumDeltas * delta;
    else
        maxFreespaceE = maxE;
    end
    
    if startPoint.north - minNumDeltas * delta < minN
        minFreespaceN = startPoint.north - minNumDeltas * delta;
    else
        minFreespaceN = minN;
    end
    
    if startPoint.north + minNumDeltas * delta > maxN
        maxFreespaceN = startPoint.north + minNumDeltas * delta;
    else
        maxFreespaceN = maxN;
    end
    
    freespace = [minFreespaceE, maxFreespaceE; minFreespaceN, maxFreespaceN];
    
    % Add LZ sides as bad areas; include a buffer of half the LZ
    maxLat = max(LandingStrip.LonLat(:,2));
    minLat = min(LandingStrip.LonLat(:,2));
    maxLon = max(LandingStrip.LonLat(:,1));
    minLon = min(LandingStrip.LonLat(:,1));
    if length(unique(LandingStrip.LonLat(:,2))) > length(unique(LandingStrip.LonLat(:,1)))
        LzLimits1 = [minLon * -1, minLat; minLon * -1, maxLat];
        LzLimits2 = [maxLon * -1, minLat; maxLon * -1, maxLat];
    else
        LzLimits1 = [minLon * -1, minLat; maxLon * -1, minLat];
        LzLimits2 = [minLon * -1, maxLat; maxLon * -1, maxLat];
    end
    maxBadAreaInd = length(badAreas);
    badAreas{maxBadAreaInd + 1} = LzLimits1;
    maxBadAreaInd = maxBadAreaInd + 1;
    badAreas{maxBadAreaInd + 1} = LzLimits2;
    maxBadAreaInd = maxBadAreaInd + 1;
    
    % Add heading limitations as bad areas; assuming heading is given as
    % azimuth angle measured clockwise from north
    % Define nominal points; assuming the UAV position as origin, and the
    % point definition as [east; north], these nominal points correspond to
    % [-15; 50], [-15; 0], [15; 0], and [15; 50] (measured in meters).
    % These will form an open rectangle bad area the UAV must navigate out
    % of.
    rotMat = [cosd(heading) sind(heading);-sind(heading) cosd(heading)];
    % Heading array values are assumed as [east; north]
    headingVals1 = rotMat * [-50; 150];
    headingVals2 = rotMat * [-50;  -50];
    headingVals3 = rotMat * [ 50;  -50];
    headingVals4 = rotMat * [ 50; 150];
    % Create points using above heading value definitions
    headingPoint1 = Waypoint(startPoint.north + headingVals1(2), startPoint.east + headingVals1(1), startPoint.down, 1, referencePoint);
    headingPoint2 = Waypoint(startPoint.north + headingVals2(2), startPoint.east + headingVals2(1), startPoint.down, 1, referencePoint);
    headingPoint3 = Waypoint(startPoint.north + headingVals3(2), startPoint.east + headingVals3(1), startPoint.down, 1, referencePoint);
    headingPoint4 = Waypoint(startPoint.north + headingVals4(2), startPoint.east + headingVals4(1), startPoint.down, 1, referencePoint);
    % Create array of badArea lon, lat combinations and add cell to badArea
    % cell array.
    % NOTE: ASSUMPTIONS ARE MADE IN THE VALIDPOINT FUNCTION SUCH THAT IT IS
    % IMPORTANT THAT THE HEADING BAD AREA BE LAST IN THE BADAREA CELL ARRAY
    headingBadArea = [
        headingPoint1.lon * -1, headingPoint1.lat;
        headingPoint2.lon * -1, headingPoint2.lat;
        headingPoint3.lon * -1, headingPoint3.lat;
        headingPoint4.lon * -1, headingPoint4.lat
        ];
    
    badAreas{maxBadAreaInd + 1} = headingBadArea;
    
    % Set all bad area point arrays to be Waypoint arrays
    badWaypoints = {};
    for i = 1:length(badAreas)
        badArea = badAreas{i};
        badWaypointSet = [];
        for j = 1:length(badArea(:,1))
            badWaypointSet = [badWaypointSet; Waypoint(badArea(j,2), badArea(j,1) * -1, 0.0, 0, referencePoint)];
        end
        badWaypoints{i} = badWaypointSet;
    end
    
    % Initialize variable for total path distance.
    totalDistance = 0.0;

    disp('Route planning algo started.')
    
    %% Begin iterating.
    for iteration = 1:numIterations
        % Initialize trigger and switch variables to determine when path has been
        % forged between start and stop points (pathCreated), and to set which
        % point will be used as the base for the next iteration (startOrStop).
        pathCreated = 0;
        startOrStop = 0; % start == 0, stop == 1

        % Initialize variables for waypoints from start and waypoints from stop.
        % These will be joined to create the final path plan after iterations
        % complete.
        Branches_start = [Branch([startPoint])];
        Branches_stop = [Branch([stopPoint])];
    
    %% Begin generating trees.
        while ~pathCreated
            % Generate random X and Y values.
            randE = freespace(1,1) + (freespace(1,2) - freespace(1,1)) * rand;
            randN = freespace(2,1) + (freespace(2,2) - freespace(2,1)) * rand;
            newRandPoint = Waypoint(randN, randE, 0.0, 1, referencePoint);

            % Grow the tree from the start or stop depending on the iteration.
            if ~startOrStop
                [closestNode, branchInd, pointInd, distance] = newRandPoint.findClosestBranchPoint(Branches_start);
            elseif startOrStop
                [closestNode, branchInd, pointInd, distance] = newRandPoint.findClosestBranchPoint(Branches_stop);
            end

            % Do not want to extend past point.
            extendDistance = delta;
%             if distance < delta
%                 extendDistance = distance;
%             else
%                 extendDistance = delta;
%             end

            % Extend out by set distance
            newEast  = closestNode.east + ((newRandPoint.east - closestNode.east) / distance) * extendDistance;
            newNorth = closestNode.north + ((newRandPoint.north - closestNode.north) / distance) * extendDistance;

            % Create new waypoint and check if new point is a valid point in
            % the freespace. Add point to appropriate waypoint array if it is
            % valid and check if a complete path has been forged.
            newPoint = Waypoint(newNorth, newEast, 0.0, 1, referencePoint);
            if ~validPoint(newPoint, closestNode, badWaypoints)
                continue
            end

            % Create new line segment to test for path completion.
            newSegment = LineSegment(closestNode, newPoint);

            %% Check which iteration is being performed.

            if ~startOrStop

                % Variable showing whether a point has been added or not
                % for this iteration; defaults to false and goes to true if
                % a point is added.
                
                pointAdded = false;
                
                % Add point to the correct branch if the closest node is the
                % final node of an existing branch.            
                % If no branches are created yet, create a new one with the
                % provided points.
                if length(Branches_start) == 0
                    Branches_start = [Branch([closestNode, newPoint])];
                    pointAdded = true;

                % Otherwise, check if the closest node is the final point on an
                % existing branch; if it is, you can add the new point as the
                % final point for that branch.
                elseif length(Branches_start(branchInd).waypoints) == pointInd
                    if length(Branches_start(branchInd).segments) == 0
                        Branches_start(branchInd) = Branches_start(branchInd).addWaypoint(newPoint);
                        pointAdded = true;
                    elseif angleCheck(newSegment, Branches_start(branchInd).segments(end), minAngle)
                        Branches_start(branchInd) = Branches_start(branchInd).addWaypoint(newPoint);
                        pointAdded = true;
                    end

                % Otherwise, the closest node must be a body point of an
                % existing branch, in which case it would create a new branch
                % from the points up to and including the closest node, with
                % the new point as the final point.
                else
                    if pointInd == 1
                        Branches_start = [Branches_start, Branch([Branches_start(branchInd).waypoints(1:pointInd), newPoint])];
                        pointAdded = true;
                    elseif angleCheck(newSegment, Branches_start(branchInd).segments(pointInd - 1), minAngle)
                        Branches_start = [Branches_start, Branch([Branches_start(branchInd).waypoints(1:pointInd), newPoint])];
                        pointAdded = true;
                    end
                end

                % If +1 points exist in opposing array, check new segment to
                % see if it intersects any segment in opposing array. If it
                % does, get intersection, add new point to both arrays, and
                % claim success. Only do this check if a point has been
                % added for this iteration.

                if pointAdded
                    
                    [goodIntersect, checkPoint, interBranchInd, interWaypointInd] = checkBranchIntersect(Branches_stop, newSegment, referencePoint, minAngle);

                    if goodIntersect
                        finalPoint = checkPoint;
                        % Trim intersected branch to point of intersection.
                        finalStopBranch = Branch(Branches_stop(interBranchInd).waypoints(1:interWaypointInd));
                        finalStopBranch = finalStopBranch.removeLastWaypoint();
                        finalStopBranch = finalStopBranch.addWaypoint(finalPoint);
                        % Remove last waypoint of intersecting branch and
                        % add new final waypoint.
                        finalStartBranch = Branches_start(branchInd);
                        finalStartBranch = finalStartBranch.removeLastWaypoint();
                        finalStartBranch = finalStartBranch.addWaypoint(finalPoint);
                        pathCreated = true;
                    end
                end
                
% %                 for i = 1:length(Branches_stop)
% %                     for j = 2:length(Branches_stop(i).waypoints)
% %                         if newSegment.doIntersect(Branches_stop(i).segments(j-1))
% %                             % If intersecting segments can be found, check
% %                             % if the resulting trimmed segments meet the
% %                             % minimum angle requirement
% %                             [checkEast, checkNorth] = newSegment.getIntersection(Branches_stop(i).segments(j-1));
% %                             checkPoint = Waypoint(checkNorth, checkEast, 0.0, 1, referencePoint);
% %                             newSegment_Trimmed = LineSegment(checkPoint, newSegment.startPoint);
% %                             existingSegment_Trimmed = LineSegment(Branches_stop(i).segments(j-1).startPoint, checkPoint);
% %                             if angleCheck(newSegment_Trimmed, existingSegment_Trimmed, minAngle)
% %                                 finalPoint = checkPoint;
% %                                 % Trim intersected branch to point of intersection.
% %                                 finalStopBranch = Branch(Branches_stop(i).waypoints(1:j));
% %                                 finalStopBranch = finalStopBranch.removeLastWaypoint();
% %                                 finalStopBranch = finalStopBranch.addWaypoint(finalPoint);
% %                                 % Remove last waypoint of intersecting branch and
% %                                 % add new final waypoint.
% %                                 finalStartBranch = Branches_start(branchInd);
% %                                 finalStartBranch = finalStartBranch.removeLastWaypoint();
% %                                 finalStartBranch = finalStartBranch.addWaypoint(finalPoint);
% %                                 pathCreated = true;
% %                             end
% %                         end
% %                     end
% %                 end

                %% 
                
                % If path still not complete, grow opposing tree in the same
                % direction. Will use similar process as above.
                if pathCreated == false
                    
                    [closestNode, branchInd, pointInd, distance] = newRandPoint.findClosestBranchPoint(Branches_stop);

                    % Do not want to extend past point.
                    extendDistance = delta;
%                     if distance < delta
%                         extendDistance = distance;
%                     else
%                         extendDistance = delta;
%                     end

                    % Extend out by set distance
                    newEast  = closestNode.east + ((newRandPoint.east - closestNode.east) / distance) * extendDistance;
                    newNorth = closestNode.north + ((newRandPoint.north - closestNode.north) / distance) * extendDistance;

                    % Create new waypoint and check if new point is a valid point in
                    % the freespace. Add point to appropriate waypoint array if it is
                    % valid and check if a complete path has been forged.
                    newPoint = Waypoint(newNorth, newEast, 0.0, 1, referencePoint);
                    if ~validPoint(newPoint, closestNode, badWaypoints)
                        continue
                    end

                    % Create new line segment to test for path completion.
                    newSegment = LineSegment(closestNode, newPoint);

                    % Variable showing whether a point has been added or not
                    % for this iteration; defaults to false and goes to true if
                    % a point is added.

                    pointAdded = false;

                    % Add point to the correct branch if the closest node is the
                    % final node of an existing branch.            
                    % If no branches are created yet, create a new one with the
                    % provided points.
                    if length(Branches_stop) == 0
                        Branches_stop = [Branch([closestNode, newPoint])];
                        pointAdded = true;

                    % Otherwise, check if the closest node is the final point on an
                    % existing branch; if it is, you can add the new point as the
                    % final point for that branch.
                    elseif length(Branches_stop(branchInd).waypoints) == pointInd
                        if length(Branches_stop(branchInd).segments) == 0
                            Branches_stop(branchInd) = Branches_stop(branchInd).addWaypoint(newPoint);
                            pointAdded = true;
                        elseif angleCheck(newSegment, Branches_stop(branchInd).segments(end), minAngle)
                            Branches_stop(branchInd) = Branches_stop(branchInd).addWaypoint(newPoint);
                            pointAdded = true;
                        end

                    % Otherwise, the closest node must be a body point of an
                    % existing branch, in which case it would create a new branch
                    % from the points up to and including the closest node, with
                    % the new point as the final point.
                    else
                        if pointInd == 1
                            Branches_stop = [Branches_stop, Branch([Branches_stop(branchInd).waypoints(1:pointInd), newPoint])];
                            pointAdded = true;
                        elseif angleCheck(newSegment, Branches_stop(branchInd).segments(pointInd - 1), minAngle)
                            Branches_stop = [Branches_stop, Branch([Branches_stop(branchInd).waypoints(1:pointInd), newPoint])];
                            pointAdded = true;
                        end
                    end

                    % If +1 points exist in opposing array, check new segment to
                    % see if it intersects any segment in opposing array. If it
                    % does, get intersection, add new point to both arrays, and
                    % claim success.

                    if pointAdded

                        [goodIntersect, checkPoint, interBranchInd, interWaypointInd] = checkBranchIntersect(Branches_start, newSegment, referencePoint, minAngle);

                        if goodIntersect
                            finalPoint = checkPoint;
                            % Trim intersected branch to point of intersection.
                            finalStartBranch = Branch(Branches_start(interBranchInd).waypoints(1:interWaypointInd));
                            finalStartBranch = finalStartBranch.removeLastWaypoint();
                            finalStartBranch = finalStartBranch.addWaypoint(finalPoint);
                            % Remove last waypoint of intersecting branch and
                            % add new final waypoint.
                            finalStopBranch = Branches_stop(branchInd);
                            finalStopBranch = finalStopBranch.removeLastWaypoint();
                            finalStopBranch = finalStopBranch.addWaypoint(finalPoint);
                            pathCreated = true;
                        end
                    end

    % %                 for i = 1:length(Branches_start)
    % %                     for j = 2:length(Branches_start(i).waypoints)
    % %                         if newSegment.doIntersect(Branches_start(i).segments(j-1))
    % %                             % If intersecting segments can be found, check
    % %                             % if the resulting trimmed segments meet the
    % %                             % minimum angle requirement
    % %                             [checkEast, checkNorth] = newSegment.getIntersection(Branches_start(i).segments(j-1));
    % %                             checkPoint = Waypoint(checkNorth, checkEast, 0.0, 1, referencePoint);
    % %                             newSegment_Trimmed = LineSegment(checkPoint, newSegment.startPoint);
    % %                             existingSegment_Trimmed = LineSegment(Branches_start(i).segments(j-1).startPoint, checkPoint);
    % %                             if angleCheck(newSegment_Trimmed, existingSegment_Trimmed, minAngle)
    % %                                 finalPoint = checkPoint;
    % %                                 % Trim intersected branch to point of intersection.
    % %                                 finalStartBranch = Branch(Branches_start(i).waypoints(1:j));
    % %                                 finalStartBranch = finalStartBranch.removeLastWaypoint();
    % %                                 finalStartBranch = finalStartBranch.addWaypoint(finalPoint);
    % %                                 % Remove last waypoint of intersecting branch and
    % %                                 % add new final waypoint.
    % %                                 finalStopBranch = Branches_stop(branchInd);
    % %                                 finalStopBranch = finalStopBranch.removeLastWaypoint();
    % %                                 finalStopBranch = finalStopBranch.addWaypoint(finalPoint);
    % %                                 pathCreated = true;
    % %                             end
    % %                         end
    % %                     end
    % %                 end
                end
            elseif startOrStop
                % Add point to waypoint array
                
                % Variable showing whether a point has been added or not
                % for this iteration; defaults to false and goes to true if
                % a point is added.
                
                pointAdded = false;

                % Add point to the correct branch if the closest node is the
                % final node of an existing branch.            
                % If no branches are created yet, create a new one with the
                % provided points.
                if length(Branches_stop) == 0
                    Branches_stop = [Branch([closestNode, newPoint])];
                    pointAdded = true;

                % Otherwise, check if the closest node is the final point on an
                % existing branch; if it is, you can add the new point as the
                % final point for that branch.
                elseif length(Branches_stop(branchInd).waypoints) == pointInd
                    if length(Branches_stop(branchInd).segments) == 0
                        Branches_stop(branchInd) = Branches_stop(branchInd).addWaypoint(newPoint);
                        pointAdded = true;
                    elseif angleCheck(newSegment, Branches_stop(branchInd).segments(end), minAngle)
                        Branches_stop(branchInd) = Branches_stop(branchInd).addWaypoint(newPoint);
                        pointAdded = true;
                    end

                % Otherwise, the closest node must be a body point of an
                % existing branch, in which case it would create a new branch
                % from the points up to and including the closest node, with
                % the new point as the final point.
                else
                    if pointInd == 1
                        Branches_stop = [Branches_stop, Branch([Branches_stop(branchInd).waypoints(1:pointInd), newPoint])];
                        pointAdded = true;
                    elseif angleCheck(newSegment, Branches_stop(branchInd).segments(pointInd - 1), minAngle)
                        Branches_stop = [Branches_stop, Branch([Branches_stop(branchInd).waypoints(1:pointInd), newPoint])];
                        pointAdded = true;
                    end
                end

                % If +1 points exist in opposing array, check new segment to
                % see if it intersects any segment in opposing array. If it
                % does, get intersection, add new point to both arrays, and
                % claim success.
                
                if pointAdded
                
                    [goodIntersect, checkPoint, interBranchInd, interWaypointInd] = checkBranchIntersect(Branches_start, newSegment, referencePoint, minAngle);

                    if goodIntersect   
                        finalPoint = checkPoint;
                        % Trim intersected branch to point of intersection.
                        finalStartBranch = Branch(Branches_start(interBranchInd).waypoints(1:interWaypointInd));
                        finalStartBranch = finalStartBranch.removeLastWaypoint();
                        finalStartBranch = finalStartBranch.addWaypoint(finalPoint);
                        % Remove last waypoint of intersecting branch and
                        % add new final waypoint.
                        finalStopBranch = Branches_stop(branchInd);
                        finalStopBranch = finalStopBranch.removeLastWaypoint();
                        finalStopBranch = finalStopBranch.addWaypoint(finalPoint);
                        pathCreated = true;
                    end
                end

% %                 for i = 1:length(Branches_start)
% %                     for j = 2:length(Branches_start(i).waypoints)
% %                         if newSegment.doIntersect(Branches_start(i).segments(j-1))
% %                             % If intersecting segments can be found, check
% %                             % if the resulting trimmed segments meet the
% %                             % minimum angle requirement
% %                             [checkEast, checkNorth] = newSegment.getIntersection(Branches_start(i).segments(j-1));
% %                             checkPoint = Waypoint(checkNorth, checkEast, 0.0, 1, referencePoint);
% %                             newSegment_Trimmed = LineSegment(checkPoint, newSegment.startPoint);
% %                             existingSegment_Trimmed = LineSegment(Branches_start(i).segments(j-1).startPoint, checkPoint);
% %                             if angleCheck(newSegment_Trimmed, existingSegment_Trimmed, minAngle)
% %                                 finalPoint = checkPoint;
% %                                 % Trim intersected branch to point of intersection.
% %                                 finalStartBranch = Branch(Branches_start(i).waypoints(1:j));
% %                                 finalStartBranch = finalStartBranch.removeLastWaypoint();
% %                                 finalStartBranch = finalStartBranch.addWaypoint(finalPoint);
% %                                 % Remove last waypoint of intersecting branch and
% %                                 % add new final waypoint.
% %                                 finalStopBranch = Branches_stop(branchInd);
% %                                 finalStopBranch = finalStopBranch.removeLastWaypoint();
% %                                 finalStopBranch = finalStopBranch.addWaypoint(finalPoint);
% %                                 pathCreated = true;
% %                             end
% %                         end
% %                     end
% %                 end

                %% 
                % If path still not complete, grow opposing tree in the same
                % direction. Will use similar process as above.

                if pathCreated == false

                    [closestNode, branchInd, pointInd, distance] = newRandPoint.findClosestBranchPoint(Branches_start);

                    % Do not want to extend past point.
                    extendDistance = delta;
%                     if distance < delta
%                         extendDistance = distance;
%                     else
%                         extendDistance = delta;
%                     end

                    % Extend out by set distance
                    newEast  = closestNode.east + ((newRandPoint.east - closestNode.east) / distance) * extendDistance;
                    newNorth = closestNode.north + ((newRandPoint.north - closestNode.north) / distance) * extendDistance;

                    % Create new waypoint and check if new point is a valid point in
                    % the freespace. Add point to appropriate waypoint array if it is
                    % valid and check if a complete path has been forged.
                    newPoint = Waypoint(newNorth, newEast, 0.0, 1, referencePoint);
                    if ~validPoint(newPoint, closestNode, badWaypoints)
                        continue
                    end

                    % Create new line segment to test for path completion.
                    newSegment = LineSegment(closestNode, newPoint);

                    % Variable showing whether a point has been added or not
                    % for this iteration; defaults to false and goes to true if
                    % a point is added.

                    pointAdded = false;

                    % Add point to the correct branch if the closest node is the
                    % final node of an existing branch.            
                    % If no branches are created yet, create a new one with the
                    % provided points.
                    if length(Branches_start) == 0
                        Branches_start = [Branch([closestNode, newPoint])];
                        pointAdded = true;

                    % Otherwise, check if the closest node is the final point on an
                    % existing branch; if it is, you can add the new point as the
                    % final point for that branch.
                    elseif length(Branches_start(branchInd).waypoints) == pointInd
                        if length(Branches_start(branchInd).segments) == 0
                            Branches_start(branchInd) = Branches_start(branchInd).addWaypoint(newPoint);
                            pointAdded = true;
                        elseif angleCheck(newSegment, Branches_start(branchInd).segments(end), minAngle)
                            Branches_start(branchInd) = Branches_start(branchInd).addWaypoint(newPoint);
                            pointAdded = true;
                        end

                    % Otherwise, the closest node must be a body point of an
                    % existing branch, in which case it would create a new branch
                    % from the points up to and including the closest node, with
                    % the new point as the final point.
                    else
                        if pointInd == 1
                            Branches_start = [Branches_start, Branch([Branches_start(branchInd).waypoints(1:pointInd), newPoint])];
                            pointAdded = true;
                        elseif angleCheck(newSegment, Branches_start(branchInd).segments(pointInd - 1), minAngle)
                            Branches_start = [Branches_start, Branch([Branches_start(branchInd).waypoints(1:pointInd), newPoint])];
                            pointAdded = true;
                        end
                    end

                    % If +1 points exist in opposing array, check new segment to
                    % see if it intersects any segment in opposing array. If it
                    % does, get intersection, add new point to both arrays, and
                    % claim success.

                    if pointAdded

                        [goodIntersect, checkPoint, interBranchInd, interWaypointInd] = checkBranchIntersect(Branches_stop, newSegment, referencePoint, minAngle);

                        if goodIntersect
                            finalPoint = checkPoint;
                            % Trim intersected branch to point of intersection.
                            finalStopBranch = Branch(Branches_stop(interBranchInd).waypoints(1:interWaypointInd));
                            finalStopBranch = finalStopBranch.removeLastWaypoint();
                            finalStopBranch = finalStopBranch.addWaypoint(finalPoint);
                            % Remove last waypoint of intersecting branch and
                            % add new final waypoint.
                            finalStartBranch = Branches_start(branchInd);
                            finalStartBranch = finalStartBranch.removeLastWaypoint();
                            finalStartBranch = finalStartBranch.addWaypoint(finalPoint);
                            pathCreated = true;
                        end
                    end

    % %                 for i = 1:length(Branches_stop)
    % %                     for j = 2:length(Branches_stop(i).waypoints)
    % %                         if newSegment.doIntersect(Branches_stop(i).segments(j-1))
    % %                             % If intersecting segments can be found, check
    % %                             % if the resulting trimmed segments meet the
    % %                             % minimum angle requirement
    % %                             [checkEast, checkNorth] = newSegment.getIntersection(Branches_stop(i).segments(j-1));
    % %                             checkPoint = Waypoint(checkNorth, checkEast, 0.0, 1, referencePoint);
    % %                             newSegment_Trimmed = LineSegment(checkPoint, newSegment.startPoint);
    % %                             existingSegment_Trimmed = LineSegment(Branches_stop(i).segments(j-1).startPoint, checkPoint);
    % %                             if angleCheck(newSegment_Trimmed, existingSegment_Trimmed, minAngle)
    % %                                 finalPoint = checkPoint;
    % %                                 % Trim intersected branch to point of intersection.
    % %                                 finalStopBranch = Branch(Branches_stop(i).waypoints(1:j));
    % %                                 finalStopBranch = finalStopBranch.removeLastWaypoint();
    % %                                 finalStopBranch = finalStopBranch.addWaypoint(finalPoint);
    % %                                 % Remove last waypoint of intersecting branch and
    % %                                 % add new final waypoint.
    % %                                 finalStartBranch = Branches_start(branchInd);
    % %                                 finalStartBranch = finalStartBranch.removeLastWaypoint();
    % %                                 finalStartBranch = finalStartBranch.addWaypoint(finalPoint);
    % %                                 pathCreated = true;
    % %                             end
    % %                         end
    % %                     end
    % %                 end

                end
            end

            % Flip startOrStop boolean.
            startOrStop = not(startOrStop);
            pointAdded = false;

        end

        % Generate final output products and check if they are better than
        % prior runs.
        
        %newTotalDistance = finalStartBranch.distance + finalStopBranch.distance;
        
        newWaypoints = [finalStartBranch.waypoints];
        for j = length(finalStopBranch.waypoints)-1:-1:1
            newWaypoints = [newWaypoints finalStopBranch.waypoints(j)];
        end

        % Perform a simple one-time pruning on the selected points
        if length(newWaypoints) > 4
            reducedWaypoints = [newWaypoints(1)];
            pointInd = 3;
            startInd = 1;
            newStartInd = 1;
            startNode = newWaypoints(startInd);
            while pointInd < length(newWaypoints)
                if startInd == 1
                    while validPoint(newWaypoints(startInd), newWaypoints(pointInd), badWaypoints) && ...
                            angleCheck(LineSegment(newWaypoints(pointInd), newWaypoints(pointInd + 1)), LineSegment(newWaypoints(startInd), newWaypoints(pointInd)), minAngle) &&...
                            pointInd < (length(newWaypoints) - 1)
                            pointInd = pointInd + 1;
                    end
                    reducedWaypoints = [reducedWaypoints, newWaypoints(pointInd - 1)];
                    startInd = pointInd - 1;
                    pointInd = pointInd + 1;
                    newStartInd = newStartInd + 1;
                else
                    while validPoint(newWaypoints(startInd), newWaypoints(pointInd), badWaypoints) && ...
                            angleCheck(LineSegment(newWaypoints(pointInd), newWaypoints(pointInd + 1)), LineSegment(newWaypoints(startInd), newWaypoints(pointInd)), minAngle) &&...
                            angleCheck(LineSegment(reducedWaypoints(newStartInd - 1), reducedWaypoints(newStartInd)), LineSegment(newWaypoints(startInd), newWaypoints(pointInd)), minAngle) &&...
                            pointInd < (length(newWaypoints) - 1)
                            pointInd = pointInd + 1;
                    end
                    reducedWaypoints = [reducedWaypoints, newWaypoints(pointInd - 1)];
                    startInd = pointInd - 1;
                    pointInd = pointInd + 1;
                    newStartInd = newStartInd + 1;
                end
            end
            if validPoint(newWaypoints(startInd), newWaypoints(end), badWaypoints) && ...
                            angleCheck(LineSegment(reducedWaypoints(newStartInd - 1), reducedWaypoints(newStartInd)), LineSegment(newWaypoints(startInd), newWaypoints(end)), minAngle)
                reducedWaypoints = [reducedWaypoints, newWaypoints(end)];
            else
                reducedWaypoints = [reducedWaypoints, newWaypoints(end - 1), newWaypoints(end)];
            end
        end
        
        newTotalDistance = 0.0;
        for j = 2:length(reducedWaypoints)
            newTotalDistance = newTotalDistance + LineSegment(reducedWaypoints(j - 1), reducedWaypoints(j)).distance;
        end
        
        if newTotalDistance < totalDistance || iteration == 1
            totalDistance = newTotalDistance;
            Waypoints = reducedWaypoints;
        end
        
        % Display progress text
        donePercent = (double(iteration)/double(numIterations)) * 100;
        statusText = ['Algo ', num2str(donePercent), '% complete.'];
        disp(statusText)
    end
    
    % Plot new full path.
    figure(2)
    eastPlot = [];
    northPlot = [];
    for j = 1:length(Waypoints)
        eastPlot = [eastPlot Waypoints(j).east];
        northPlot = [northPlot Waypoints(j).north];
    end
    plot(eastPlot, northPlot, '--xk')
    hold on
    for i = 1:length(badWaypoints)
        badAreaBoundary = [];
        badWaypointSet = badWaypoints{i};
        for j = 1:length(badWaypointSet)
            badAreaBoundary = [badAreaBoundary; badWaypointSet(j).east, badWaypointSet(j).north];
        end
        if i ~= length(badWaypoints)
            badAreaBoundary = [badAreaBoundary; badWaypointSet(1).east, badWaypointSet(1).north];
        end
        plot(badAreaBoundary(:,1), badAreaBoundary(:,2), '--r')
    end
    plot(startPoint.east, startPoint.north, '^b');
    plot(stopPoint.east, stopPoint.north, 'xg');
    xlabel('East (m)')
    ylabel('North (m)')
    title('Final Waypoint Path NED Plot')
    % Set axes to be equivalent length to preserve visual ratio
    axesDeltaE = abs(freespace(1,1) - freespace(1,2));
    axesDeltaN = abs(freespace(2,1) - freespace(2,2));
    axesDelta = abs(axesDeltaE - axesDeltaN) * 0.5;
    if axesDeltaE > axesDeltaN
        freespace(2,2) = freespace(2,2) + axesDelta;
        freespace(1,2) = freespace(1,2) - axesDelta;
    elseif axesDeltaN > axesDeltaE
        freespace(1,2) = freespace(1,2) + axesDelta;
        freespace(1,1) = freespace(1,1) - axesDelta;
    end
    axis([freespace(1,1) freespace(1,2) freespace(2,1) freespace(2,2)])
    
    h = zeros(4, 1);
    h(1) = plot(NaN,NaN,'--xk');
    h(2) = plot(NaN,NaN,'--r');
    h(3) = plot(NaN,NaN,'^b');
    h(4) = plot(NaN,NaN,'xg');
    legend(h, 'Planned Path','Restricted Regions','UAV Current Position','Landing Strip Center', 'Location', 'best');
    
    % Plot informational text
    text1 = ['Delta Distance: ', num2str(delta), ' m'];
    text2 = ['Minimum Turn Angle: ', num2str(minAngle), ' deg'];
    text3 = ['UAV Heading: ', num2str(heading), ' deg'];
    text4 = ['Number of Iterations: ', num2str(numIterations)];
    plotText = {text1, text2, text3, text4};
    text(freespace(1,1) + 50,freespace(2,1) + 400,plotText)
    hold off
    
    
    % Plot new full path.
    figure(3)
    eastPlot = [];
    northPlot = [];
    for j = 1:length(newWaypoints)
        eastPlot = [eastPlot newWaypoints(j).east];
        northPlot = [northPlot newWaypoints(j).north];
    end
    plot(eastPlot, northPlot, '--xk')
    hold on
    for i = 1:length(badWaypoints)
        badAreaBoundary = [];
        badWaypointSet = badWaypoints{i};
        for j = 1:length(badWaypointSet)
            badAreaBoundary = [badAreaBoundary; badWaypointSet(j).east, badWaypointSet(j).north];
        end
        if i ~= length(badWaypoints)
            badAreaBoundary = [badAreaBoundary; badWaypointSet(1).east, badWaypointSet(1).north];
        end
        plot(badAreaBoundary(:,1), badAreaBoundary(:,2), '--r')
    end
    plot(startPoint.east, startPoint.north, '^b');
    plot(stopPoint.east, stopPoint.north, 'xg');
    xlabel('East (m)')
    ylabel('North (m)')
    title('Final Waypoint Path NED Plot')
    % Set axes to be equivalent length to preserve visual ratio
    axesDeltaE = abs(freespace(1,1) - freespace(1,2));
    axesDeltaN = abs(freespace(2,1) - freespace(2,2));
    axesDelta = abs(axesDeltaE - axesDeltaN) * 0.5;
    if axesDeltaE > axesDeltaN
        freespace(2,2) = freespace(2,2) + axesDelta;
        freespace(1,2) = freespace(1,2) - axesDelta;
    elseif axesDeltaN > axesDeltaE
        freespace(1,2) = freespace(1,2) + axesDelta;
        freespace(1,1) = freespace(1,1) - axesDelta;
    end
    axis([freespace(1,1) freespace(1,2) freespace(2,1) freespace(2,2)])
    
    h = zeros(4, 1);
    h(1) = plot(NaN,NaN,'--xk');
    h(2) = plot(NaN,NaN,'--r');
    h(3) = plot(NaN,NaN,'^b');
    h(4) = plot(NaN,NaN,'xg');
    legend(h, 'Planned Path','Restricted Regions','UAV Current Position','Landing Strip Center', 'Location', 'best');
    
    % Plot informational text
    text1 = ['Delta Distance: ', num2str(delta), ' m'];
    text2 = ['Minimum Turn Angle: ', num2str(minAngle), ' deg'];
    text3 = ['UAV Heading: ', num2str(heading), ' deg'];
    text4 = ['Number of Iterations: ', num2str(numIterations)];
    plotText = {text1, text2, text3, text4};
    text(freespace(1,1) + 50,freespace(2,1) + 400,plotText)
    hold off
    
    Waypoints = newWaypoints;
    
    % Revert all longitudes
    for i = 1:length(Waypoints)
        Waypoints(i).lon = Waypoints(i).lon * -1;
    end
end

