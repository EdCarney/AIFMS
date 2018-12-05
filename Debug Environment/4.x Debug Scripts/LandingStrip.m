classdef LandingStrip
    %LANDINGSTRIP A class to hold all LZ data.
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
    
    properties
        Indices
        LonLat
        DEM
        TRI
        LC
        avgDEM
        stdDEM
        maxDiffDEM
        avgTRI
        avgLC
        centroid
        distFromOtherStrips
        distFromUAV
        score
        deg2met
        validStrip
    end
    
    % Will define all weightings as constant properties.
    properties (Constant)
        W_avgTri = -15; % Weighting for average ruggedness (lower values preferred)
        W_avgDem = -15; % Weighting for average elevation (lower altitude is preferred)
        W_LcType = -50; % Weighting for land cover type (prefer bare earth, vegetation, and farm land in that order)
        W_maxDiffDEM = -15; % Weighting for max elevation difference from average (lower values preferred)
        W_distFromOtherStrips = -0.01; % Weighting for proximity to other strips (higher values preferred)
        W_distFromUAV = -0.1; % Weighting for ditance to UAV (lower values preferred)
    end
        
    
    methods
        function obj = LandingStrip(Indices, LonLat, DEM, TRI, LC, x_UAV, y_UAV, deg2met)
            %LANDINGSTRIP Construct an instance of this class
            if nargin > 0
                obj.Indices = Indices;
                obj.LonLat = LonLat;
                obj.DEM = DEM;
                obj.TRI = TRI;
                obj.LC = LC;
                obj.avgDEM = mean(DEM);
                obj.stdDEM = std(DEM);
                obj.maxDiffDEM = max(DEM) - obj.avgDEM;
                obj.avgTRI = mean(TRI);
                obj.avgLC = mean(LC);
                obj.centroid = [mean(LonLat(:,1)) mean(LonLat(:,2))];
                obj.distFromOtherStrips = 0;
                obj.distFromUAV = sqrt((obj.centroid(1) - x_UAV)^2 + (obj.centroid(2) - y_UAV)^2) * deg2met; % meter units
                obj.score = 0;
                obj.deg2met = deg2met;
                obj.validStrip = 1;
            end
        end
        
        function  distance = computeDistance(obj,LandingStrip)
            %computeDistance Method to compute distance between this strip
            %and another (returns degree value).
            distance = sqrt((obj.centroid(1) - LandingStrip.centroid(1))^2 + (obj.centroid(2) - LandingStrip.centroid(2))^2);
            distance = distance * obj.deg2met;
        end
                
        function score = computeScore(obj)
            %computeScore Method to compute the overall strip score after
            %setting strip non-specific values.
            score = obj.avgDEM * obj.W_avgDem + obj.maxDiffDEM * obj.W_maxDiffDEM + ...
                obj.avgTRI * obj.W_avgTri + obj.distFromOtherStrips * obj.W_distFromOtherStrips + ...
                obj.distFromUAV * obj.W_distFromUAV + obj.avgLC * obj.W_LcType;
        end
    end
end

