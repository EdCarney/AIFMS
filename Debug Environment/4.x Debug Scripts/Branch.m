classdef Branch
    %BRANCH Class to represent multiple points forming a full branch.
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
    %   Properties include an array of n waypoint objects forming the
    %   branch, and an array of n-1 line segments constructed by the same
    %   points. Accessible branch attributes include the total straight line
    %   distance along the branch path.
    
    properties
        waypoints
        segments
        distance
    end
    
    methods
        function obj = Branch(waypointArray)
            %BRANCH Construct an instance of this class
            %   Branches are instantiated with an array of waypoint.
            
            obj.distance = 0.0;
            obj.waypoints = waypointArray;
            if length(waypointArray) > 1
                for i = 1:length(waypointArray) - 1
                    tempSegments = obj.segments;
                    obj.segments = [tempSegments,...
                        LineSegment(waypointArray(i), waypointArray(i + 1))];
                end
            end
            
            for i = 1:length(obj.segments)
                tempDistance = obj.distance;
                obj.distance = tempDistance + obj.segments(i).distance;
            end
        end
        
        function obj = addWaypoint(obj,waypoint)
            %ADDWAYPOINT Method to add a waypoint to the branch.
            %   Automatically adds a new line segment and recomputes the
            %   total branch distance.
            tempWaypoints = obj.waypoints;
            obj.waypoints = [tempWaypoints waypoint];
            
            tempSegments = obj.segments;
            obj.segments = [tempSegments,...
                LineSegment(obj.waypoints(end - 1), obj.waypoints(end))];
            
            obj.distance = obj.distance + obj.segments(end).distance;
            
        end
        
        function obj = removeLastWaypoint(obj)
            %REMOVELASTWAYPOINT Method to remove the last waypoint of a
            %branch.
            %   Automatically removes old line segment and recomputes the
            %   total branch distance.
            obj.distance = obj.distance - obj.segments(end).distance;
            obj.segments = obj.segments(1:end-1);
            obj.waypoints = obj.waypoints(1:end-1);
            
        end
        
    end
end

