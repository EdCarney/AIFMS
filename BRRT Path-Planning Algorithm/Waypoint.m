classdef Waypoint
    %WAYPOINT Basic class to hold data on waypoints.
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
    %   Holds NED and LLA information on datapoint. Initialized via
    %   Waypoint(arg1, arg2, arg3, LLAorNED, referencePoint), where
    %   arguments arg1, ag2, arg3 define the waypoint values either
    %   in NED or LLA depending on the value or argument LLAorNED.
    %   LLAorNED == 0 indicates args will be Lat, Lon, Alt; LLAorNED == 1
    %   indicates args will be North, East, Down. The referencePoint sets
    %   the reference point for the LLA-NED conversion; given in LLA.
    
    properties
        north
        east
        down
        lat
        lon
        alt
    end
    
    methods
        
        function obj = Waypoint(arg1, arg2, arg3, LLAorNED, referencePoint)
            %WAYPOINT Construct an instance of this class.
            %   Arguments arg1, ag2, arg3 define the waypoint values either
            %   in NED or LLA depending on the value or argument LLAorNED.
            %   LLAorNED == 0 indicates args will be LLA; LLAorNED == 1
            %   indicated args will be NED. The referencePoint sets the
            %   reference point for the LLA-NED conversion; given in LLA.
            if LLAorNED
                obj.north = arg1;
                obj.east  = arg2;
                obj.down  = arg3;
                [obj.lat, obj.lon, obj.alt] = ned2geodetic(arg1, arg2,...
                    arg3, referencePoint(1), referencePoint(2), ...
                    referencePoint(3), wgs84Ellipsoid);
            else
                obj.lat = arg1;
                obj.lon = arg2;
                obj.alt = arg3;
                [obj.north, obj.east, obj.down] = geodetic2ned(arg1, arg2,...
                    arg3, referencePoint(1), referencePoint(2), ...
                    referencePoint(3), wgs84Ellipsoid);
            end
        end
        
        function  distance = compute2dDistance(obj,Waypoint)
            %computeDistance Method to compute distance between this point
            %and another (returns meter value).
            distance = sqrt((obj.north - Waypoint.north)^2 + (obj.east - Waypoint.east)^2);
        end
        
        function [closetPoint, pointInd, minDist] = findClosestPoint(obj,Waypoints)
            %computeDistance Method to find the closest point to the
            %current point given a list of points to check against.
            
            % Consider trival case.
            if length(Waypoints) == 1
                closetPoint = Waypoints(1);
                pointInd = 1;
                minDist = obj.compute2dDistance(Waypoints(1));
                
            % Otherwise loop through list of points and find closest.
            else
                minDist = 0;
                for i = 1:length(Waypoints)
                    newDistance = obj.compute2dDistance(Waypoints(i));
                    if newDistance < minDist || i == 1
                        closetPoint = Waypoints(i);
                        pointInd = i;
                        minDist = newDistance;
                    end
                end
            end
        end
        
        function [closestPoint, branchInd, pointInd, minDist] = findClosestBranchPoint(obj,Branches)
            
            % Consider trivial case when there are no branches.
            if length(Branches) == 0
                returnWaypint = obj;
                closestPoint = returnWaypint;
                branchInd = 0;
                pointInd = 0;
            % Iterate through each branch and get the closest point.
            else            
                minDist = 0.0;
                for i = 1:length(Branches)
                    [newClosestPoint, newPointInd, newMinDist] = obj.findClosestPoint(Branches(i).waypoints);
                    if newMinDist < minDist || i == 1
                        minDist = newMinDist;
                        closestPoint = newClosestPoint;
                        pointInd = newPointInd;
                        branchInd = i;
                    end
                end
            end
        end
    end
end

