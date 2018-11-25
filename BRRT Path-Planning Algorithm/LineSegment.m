classdef LineSegment
    %LINESEGMENT Basic class to hold data on line segments formed by
    %waypoint classes.
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
    %   Takes two waypoints upon instantiation and can be used to determine
    %   line segment properties and intersections with other line segments.
    
    properties
        startPoint
        stopPoint
        distance
    end
    
    methods (Access = public)
        function obj = LineSegment(startPoint, stopPoint)
            %LINESEGMENT Construct an instance of this class using two
            %waypoint objects.
            obj.startPoint = startPoint;
            obj.stopPoint = stopPoint;
            obj.distance = startPoint.compute2dDistance(stopPoint);
        end
        
        function intersect = doIntersect(obj, LineSeg)
            %DOINTERSECT Determines whether a given line segment intersects
            %with this line segment. Outputs a true/false boolean depending
            %on the answer.
            
            % Get four orientations needed for all cases.
            orientation1 = obj.getOrientation(obj.startPoint, obj.stopPoint, LineSeg.startPoint);
            orientation2 = obj.getOrientation(obj.startPoint, obj.stopPoint, LineSeg.stopPoint);
            orientation3 = obj.getOrientation(LineSeg.startPoint, LineSeg.stopPoint, obj.startPoint);
            orientation4 = obj.getOrientation(LineSeg.startPoint, LineSeg.stopPoint, obj.stopPoint);
        
            % Check if general intersection case is true.
            if orientation1 ~= orientation2 && orientation3 ~= orientation4
                intersect = true;
            elseif orientation1 == 0 && obj.onSegment(obj.startPoint, obj.stopPoint, LineSeg.startPoint)
                intersect = true;
            elseif orientation2 == 0 && obj.onSegment(obj.startPoint, obj.stopPoint, LineSeg.stopPoint)
                intersect = true;
            elseif orientation3 == 0 && obj.onSegment(LineSeg.startPoint, LineSeg.stopPoint, obj.startPoint)
                intersect = true;
            elseif orientation4 == 0 && obj.onSegment(LineSeg.startPoint, LineSeg.stopPoint, obj.stopPoint)
                intersect = true;
            else
                intersect = false;
            end
        end
        
        function [eastIntersect, northIntersect] = getIntersection(obj, LineSeg)
            %GETINTERSECTION Calculates intersection point of two line
            %segments. Returns waypoint object.
            
            % Parent line segment as A1*x + B1*y = C1.
            A1 = obj.stopPoint.north - obj.startPoint.north;
            B1 = obj.startPoint.east - obj.stopPoint.east;
            C1 = A1 * obj.startPoint.east + B1 * obj.startPoint.north;
            
            % Comparision line segment as A2*x + B2*y = C2.
            A2 = LineSeg.stopPoint.north - LineSeg.startPoint.north;
            B2 = LineSeg.startPoint.east - LineSeg.stopPoint.east;
            C2 = A2 * LineSeg.startPoint.east + B2 * LineSeg.startPoint.north;
            
            % Calculate intersection.
            determinant = A1*B2 - A2*B1;
            eastIntersect  = (B2 * C1 - B1 * C2) / determinant;
            northIntersect = (A1 * C2 - A2 * C1) / determinant;
            
        end
    end
    methods (Access = protected)
        function onSegment = onSegment(obj, startPoint, stopPoint, Waypoint)
            %ONSEGMENT Determines if a given waypoing lies on the line
            %segment.
            if Waypoint.east <= max([startPoint.east,stopPoint.east])...
                    && Waypoint.east >= min([startPoint.east,stopPoint.east])...
                    && Waypoint.north <= max([startPoint.north,stopPoint.north])...
                    && Waypoint.north >= min([startPoint.north,stopPoint.north]);
                onSegment = true;
            else
                onSegment = false;
            end
        end
        
        function orientation = getOrientation(obj, startPoint, stopPoint, Waypoint)
            %GETORIENTATION Determines the orientation of three coplanar
            %points. Orientation are set as 0 == coplanar, 1 == clockwise,
            %2 == counterclockwise.
            val = (Waypoint.north - startPoint.north) * (stopPoint.east...
                - Waypoint.east) - (Waypoint.east - startPoint.east) * ...
                (stopPoint.north - Waypoint.north);
            if val == 0
                orientation = 0;
            elseif val > 0
                orientation = 1;
            else
                orientation = 2;
            end
        end
    end
end

