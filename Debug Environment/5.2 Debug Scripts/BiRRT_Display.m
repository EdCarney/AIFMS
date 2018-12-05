function [Waypoints, totalDistance] = BiRRT_Display(LandingStrip, currentPosition, alt, delta, numIterations, badAreas)
%BIRRT Builds a path between the LzStrip and UAV using a Bi-directional RRT
%algoirthm.
%   Takes three arguments, an LandingStrip object, the current UAV posiiton
%   as a 1 x 2 vector of [Lon Lat], the altitude as the current UAV
%   altitude, and the delta value (in meters) for the growth of the tree.
%   Outputs an n x 2 matrix, where n is the number of waypoints between the
%   UAV position and LZ strip centriod, and the 2 columns corresponding to
%   the waypoint Lon and Lat, respectively.

    % Will define the centriod of the LzStrip as the reference point for all
    % conversions to NED; defined as Lat, Lon, Alt.
    referencePoint = [LandingStrip.centroid(2),LandingStrip.centroid(1), LandingStrip.avgDEM];

    % Define start and stop points for the algorithm.
    startPoint = Waypoint(currentPosition(2), currentPosition(1), alt, 0, referencePoint);
    stopPoint  = Waypoint(LandingStrip.centroid(2), LandingStrip.centroid(1), LandingStrip.avgDEM, 0, referencePoint);

    % Set domain of algorithm freespace as [minE maxE ; minN maxN].
    minE = Waypoint(referencePoint(1), min(min(LandingStrip.LonLat(:,1), startPoint.lon)), 0.0, 0, referencePoint).east;
    maxE = Waypoint(referencePoint(1), max(max(LandingStrip.LonLat(:,1), startPoint.lon)), 0.0, 0, referencePoint).east;
    minN = Waypoint(min(min(LandingStrip.LonLat(:,2), startPoint.lat)), referencePoint(2), 0.0, 0, referencePoint).north;
    maxN = Waypoint(max(max(LandingStrip.LonLat(:,2), startPoint.lat)), referencePoint(2), 0.0, 0, referencePoint).north;
    freespace = [minE maxE ; minN maxN];
    
    % Set all bad area point arrays to be Waypoint arrays
    badWaypoints = {};
    for i = 1:length(badAreas)
        badArea = badAreas{i};
        badWaypointSet = [];
        for j = 1:length(badArea(:,1))
            badWaypointSet = [badWaypointSet; Waypoint(badArea(j,2), badArea(j,1), 0.0, 0, referencePoint)];
        end
        badWaypoints{i} = badWaypointSet;
    end
    
    % Initialize variable for total path distance.
    totalDistance = 0.0;

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

         % Initiate plot for points.
         close all
         figure(1)
         hold on
    
    %% Begin generating trees.
        while ~pathCreated
            % Generate random X and Y values.
            randE = minE + (maxE-minE)*rand;
            randN = minN + (maxN-minN)*rand;
            newRandPoint = Waypoint(randN, randE, 0.0, 1, referencePoint);

            % Grow the tree from the start or stop depending on the iteration.
            if ~startOrStop
                [closestNode, branchInd, pointInd, distance] = newRandPoint.findClosestBranchPoint(Branches_start);
            elseif startOrStop
                [closestNode, branchInd, pointInd, distance] = newRandPoint.findClosestBranchPoint(Branches_stop);
            end

            % Do not want to extend past point.
            if distance < delta
                extendDistance = distance;
            else
                extendDistance = delta;
            end

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

                % Add point to the correct branch if the closest node is the
                % final node of an existing branch.            
                % If no branches are created yet, create a new one with the
                % provided points.
                if length(Branches_start) == 0
                    Branches_start = [Branch([closestNode, newPoint])];

                % Otherwise, check if the closest node is the final point on an
                % existing branch; if it is, you can add the new point as the
                % final point for that branch.
                elseif length(Branches_start(branchInd).waypoints) == pointInd
                    Branches_start(branchInd) = Branches_start(branchInd).addWaypoint(newPoint);

                % Otherwise, the closest node must be a body point of an
                % existing branch, in which case it would create a new branch
                % from the points up to and including the closest node, with
                % the new point as the final point.
                else
                    Branches_start = [Branches_start, Branch([Branches_start(branchInd).waypoints(1:pointInd), newPoint])];
                end

                % If +1 points exist in opposing array, check new segment to
                % see if it intersects any segment in opposing array. If it
                % does, get intersection, add new point to both arrays, and
                % claim success.

                for i = 1:length(Branches_stop)
                    for j = 2:length(Branches_stop(i).waypoints)
                        if newSegment.doIntersect(Branches_stop(i).segments(j-1))
                            [finalEast, finalNorth] = newSegment.getIntersection(Branches_stop(i).segments(j-1));
                            finalPoint = Waypoint(finalNorth, finalEast, 0.0, 1, referencePoint);
                            % Trim intersected branch to point of intersection.
                            finalStopBranch = Branch(Branches_stop(i).waypoints(1:j));
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
                end

                 % Plot new point.
                 plot(newPoint.east, newPoint.north, 'bx');
                 hold on;

                %% 
                % If path still not complete, grow opposing tree in the same
                % direction. Will use similar process as above.

                [closestNode, branchInd, pointInd, distance] = newRandPoint.findClosestBranchPoint(Branches_stop);

                % Do not want to extend past point.
                if distance < delta
                    extendDistance = distance;
                else
                    extendDistance = delta;
                end

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

                % Add point to the correct branch if the closest node is the
                % final node of an existing branch.            
                % If no branches are created yet, create a new one with the
                % provided points.
                if length(Branches_stop) == 0
                    Branches_stop = [Branch([closestNode, newPoint])];

                % Otherwise, check if the closest node is the final point on an
                % existing branch; if it is, you can add the new point as the
                % final point for that branch.
                elseif length(Branches_stop(branchInd).waypoints) == pointInd
                    Branches_stop(branchInd) = Branches_stop(branchInd).addWaypoint(newPoint);

                % Otherwise, the closest node must be a body point of an
                % existing branch, in which case it would create a new branch
                % from the points up to and including the closest node, with
                % the new point as the final point.
                else
                    Branches_stop = [Branches_stop, Branch([Branches_stop(branchInd).waypoints(1:pointInd), newPoint])];
                end

                % If +1 points exist in opposing array, check new segment to
                % see if it intersects any segment in opposing array. If it
                % does, get intersection, add new point to both arrays, and
                % claim success.

                for i = 1:length(Branches_start)
                    for j = 2:length(Branches_start(i).waypoints)
                        if newSegment.doIntersect(Branches_start(i).segments(j-1))
                            [finalEast, finalNorth] = newSegment.getIntersection(Branches_start(i).segments(j-1));
                            finalPoint = Waypoint(finalNorth, finalEast, 0.0, 1, referencePoint);
                            % Trim intersected branch to point of intersection.
                            finalStartBranch = Branch(Branches_start(i).waypoints(1:j));
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
                end

                 % Plot new point.
                 plot(newPoint.east, newPoint.north, 'gx');
                 hold on;

            elseif startOrStop
                % Add point to waypoint array

                % Add point to the correct branch if the closest node is the
                % final node of an existing branch.            
                % If no branches are created yet, create a new one with the
                % provided points.
                if length(Branches_stop) == 0
                    Branches_stop = [Branch([closestNode, newPoint])];

                % Otherwise, check if the closest node is the final point on an
                % existing branch; if it is, you can add the new point as the
                % final point for that branch.
                elseif length(Branches_stop(branchInd).waypoints) == pointInd
                    Branches_stop(branchInd) = Branches_stop(branchInd).addWaypoint(newPoint);

                % Otherwise, the closest node must be a body point of an
                % existing branch, in which case it would create a new branch
                % from the points up to and including the closest node, with
                % the new point as the final point.
                else
                    Branches_stop = [Branches_stop, Branch([Branches_stop(branchInd).waypoints(1:pointInd), newPoint])];
                end

                % If +1 points exist in opposing array, check new segment to
                % see if it intersects any segment in opposing array. If it
                % does, get intersection, add new point to both arrays, and
                % claim success.

                for i = 1:length(Branches_start)
                    for j = 2:length(Branches_start(i).waypoints)
                        if newSegment.doIntersect(Branches_start(i).segments(j-1))
                            [finalEast, finalNorth] = newSegment.getIntersection(Branches_start(i).segments(j-1));
                            finalPoint = Waypoint(finalNorth, finalEast, 0.0, 1, referencePoint);
                            % Trim intersected branch to point of intersection.
                            finalStartBranch = Branch(Branches_start(i).waypoints(1:j));
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
                end

                 % Plot new point.
                 plot(newPoint.east, newPoint.north, 'gx');
                 hold on;

                %% 
                % If path still not complete, grow opposing tree in the same
                % direction. Will use similar process as above.

                [closestNode, branchInd, pointInd, distance] = newRandPoint.findClosestBranchPoint(Branches_start);

                % Do not want to extend past point.
                if distance < delta
                    extendDistance = distance;
                else
                    extendDistance = delta;
                end

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

                % Add point to the correct branch if the closest node is the
                % final node of an existing branch.            
                % If no branches are created yet, create a new one with the
                % provided points.
                if length(Branches_start) == 0
                    Branches_start = [Branch([closestNode, newPoint])];

                % Otherwise, check if the closest node is the final point on an
                % existing branch; if it is, you can add the new point as the
                % final point for that branch.
                elseif length(Branches_start(branchInd).waypoints) == pointInd
                    Branches_start(branchInd) = Branches_start(branchInd).addWaypoint(newPoint);

                % Otherwise, the closest node must be a body point of an
                % existing branch, in which case it would create a new branch
                % from the points up to and including the closest node, with
                % the new point as the final point.
                else
                    Branches_start = [Branches_start, Branch([Branches_start(branchInd).waypoints(1:pointInd), newPoint])];
                end

                % If +1 points exist in opposing array, check new segment to
                % see if it intersects any segment in opposing array. If it
                % does, get intersection, add new point to both arrays, and
                % claim success.

                for i = 1:length(Branches_stop)
                    for j = 2:length(Branches_stop(i).waypoints)
                        if newSegment.doIntersect(Branches_stop(i).segments(j-1))
                            [finalEast, finalNorth] = newSegment.getIntersection(Branches_stop(i).segments(j-1));
                            finalPoint = Waypoint(finalNorth, finalEast, 0.0, 1, referencePoint);
                            % Trim intersected branch to point of intersection.
                            finalStopBranch = Branch(Branches_stop(i).waypoints(1:j));
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
                end

                 % Plot new point.
                 plot(newPoint.east, newPoint.north, 'bx');
                 hold on;
            end

             for i = 1:length(Branches_stop)
                 eastPlot = [];
                 northPlot = [];
                 for j = 1:length(Branches_stop(i).waypoints)
                     eastPlot = [eastPlot, Branches_stop(i).waypoints(j).east];
                     northPlot = [northPlot, Branches_stop(i).waypoints(j).north];
                 end
                 plot(eastPlot, northPlot, 'g--')
             end
 
             for i = 1:length(Branches_start)
                 eastPlot = [];
                 northPlot = [];
                 for j = 1:length(Branches_start(i).waypoints)
                     eastPlot = [eastPlot, Branches_start(i).waypoints(j).east];
                     northPlot = [northPlot, Branches_start(i).waypoints(j).north];
                 end
                 plot(eastPlot, northPlot, 'b--')
             end

            % Flip startOrStop boolean.
            startOrStop = not(startOrStop);

        % Display specific stuff
        hold on
        for i = 1:length(badWaypoints)
            badAreaBoundary = [];
            badWaypointSet = badWaypoints{i};
            for j = 1:length(badArea)
                badAreaBoundary = [badAreaBoundary; badWaypointSet(j).east, badWaypointSet(j).north];
            end
            badAreaBoundary = [badAreaBoundary; badWaypointSet(1).east, badWaypointSet(1).north];
            plot(badAreaBoundary(:,1), badAreaBoundary(:,2), '--r')
        end
        plot(startPoint.east, startPoint.north, '^k');
        plot(stopPoint.east, stopPoint.north, 'xk');
        xlabel('East (m)')
        ylabel('North (m)')
            
        end

         xlabel('East (m)')
         ylabel('North (m)')
         title('Waypoint NED Plot')
         axis([minE maxE minN maxN])
         hold off

        % Generate final output products and check if they are better than
        % prior runs.
        newTotalDistance = finalStartBranch.distance + finalStopBranch.distance;
        newWaypoints = [finalStartBranch.waypoints];
        for j = length(finalStopBranch.waypoints)-1:-1:1
            newWaypoints = [newWaypoints finalStopBranch.waypoints(j)];
        end

        if newTotalDistance < totalDistance || iteration == 1
            totalDistance = newTotalDistance;
            Waypoints = newWaypoints;
        end
        
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
        for j = 1:length(badArea)
            badAreaBoundary = [badAreaBoundary; badWaypointSet(j).east, badWaypointSet(j).north];
        end
        badAreaBoundary = [badAreaBoundary; badWaypointSet(1).east, badWaypointSet(1).north];
        plot(badAreaBoundary(:,1), badAreaBoundary(:,2), '--r')
    end
    plot(startPoint.east, startPoint.north, '^k');
    plot(stopPoint.east, stopPoint.north, 'xk');
    xlabel('East (m)')
    ylabel('North (m)')
    title('Final Waypoint Path NED Plot')
    axis([minE maxE minN maxN])
    hold off
end

