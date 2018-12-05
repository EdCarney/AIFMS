%%
% Created Based on Concepts From:
%   Files:
%       LZ_Calculation_FGIF.m
%       LandingStrip.m
%   Author:
%       Edward Carney
%       Systems Engineer, Orbit Logic Inc.
%       Graduate Student, Aerospace Engineering, University of Maryland
%   Creation Date:
%       24 Sep 2018

function bestStrip = GetLandingZones(N_new, x_earth, y_earth, deg2met, Alt)

    tic()

    % Define minimum LZ length and width and buffer size.
    minLength = 3;
    minWidth = 1;
    buffer = 1;

    % Set N_new so that populated zones are 1 and unpopulated are 0.
    N_new(:,:,1) = N_new(:,:,1)==10000;
    
    % Find the size of N_new; will be used to limit the for-loop iterations
    [row_size,col_size] = size(N_new(:,:,1));

    % Get the indices for areas of FGIF that are unpopulated and use them
    % to create the SRM matrix
    [row,col] = find(~N_new(:,:,1));
    SRM = [row,col];
    len = length(SRM);

    % Initialize variables in for-loop; will set size of LZs to be the
    % maximum possible and then trim after (more efficient than dynamic
    % allocation)
    LZs = zeros((minWidth + 2*buffer) * (minLength + 2*buffer),2,row_size * col_size * (minLength + 2*buffer));
    numLZs = 0;
    
    for k=1:len

        % Set entry indices to check.
        i_row = SRM(k,1);
        i_col = SRM(k,2);

        % ROW CHECK:
        % Check for any horizontal LZs (check along row); will check for
        % length along column and width along row
        for j = -(minLength + 2*buffer - 1):1:0
            
            % Set minimums and maximums for this iteration
            col_min = i_col + j;
            col_max = i_col + j + minLength + 2*buffer - 1;
            row_min = i_row - (buffer + minWidth - 1);
            row_max = i_row + (buffer + minWidth - 1);
            
            % Skip over iterations that are out of range of N_new
            if (col_min) < 1 || (col_max) > col_size
                continue
            end
            if (row_min) < 1 || (row_max) > row_size
                continue
            end
            
            % Check if any entries are populated; if they are, exclude that
            % region from consideration; else, add that region to the row
            % LZ list.
            if any(any(N_new(row_min:row_max,col_min:col_max,1)))
                continue
            else
                
                % Check number of LZs to determine where to add matrix of
                % LZ points
                numLZs = ceil(nnz(LZs)/(((minWidth + 2*buffer) * (minLength + 2*buffer))*2));
                
                % Initialize new LZ matrix and row and column arrays
                LZ = zeros((minWidth + 2*buffer) * (minLength + 2*buffer),2,1);
                row_array = [0 0];
                col_array = [0 0];
                
                % Set row and column arrays to be the same as those for the
                % ANY(ANY()) check above.
                row_array = row_min:row_max;
                col_array = col_min:col_max;
                
                % Populate the entries of the new LZ                
                for row_entry = 1:length(row_array)
                    for col_entry = 1:length(col_array)
                        LZ((row_entry - 1) * length(col_array) + col_entry,:) = [row_array(row_entry),col_array(col_entry)];
                    end
                end
                
                % Add new LZ to 3D LZs matrix
                LZs(:,:,numLZs + 1) = LZ;
                
            end
        end
        
        
        % COLUMN CHECK:
        % Check for any vertical LZs (check along row); will check for
        % length along row and width along column
        for j = -(minLength + 2*buffer - 1):1:0
            
            % Set minimums and maximums for this iteration
            row_min = i_row + j;
            row_max = i_row + j + minLength + 2*buffer - 1;
            col_min = i_col - (buffer + minWidth - 1);
            col_max = i_col + (buffer + minWidth - 1);
            
            % Skip over iterations that are out of range of N_new
            if (col_min) < 1 || (col_max) > col_size
                continue
            end
            if (row_min) < 1 || (row_max) > row_size
                continue
            end
            
            % Check if any entries are populated; if they are, exclude that
            % region from consideration; else, add that region to the row
            % LZ list.
            if any(any(N_new(row_min:row_max,col_min:col_max,1)))
                continue
            else
                
                % Check number of LZs to determine where to add matrix of
                % LZ points
                numLZs = ceil(nnz(LZs)/(((minWidth + 2*buffer) * (minLength + 2*buffer))*2));
                
                % Initialize new LZ matrix and row and column arrays
                LZ = zeros((minWidth + 2*buffer) * (minLength + 2*buffer),2,1);
                row_array = [0 0];
                col_array = [0 0];
                
                % Set row and column arrays to be the same as those for the
                % ANY(ANY()) check above.
                row_array = row_min:row_max;
                col_array = col_min:col_max;
                
                % Populate the entries of the new LZ                
                for row_entry = 1:length(row_array)
                    for col_entry = 1:length(col_array)
                        LZ((row_entry - 1) * length(col_array) + col_entry,:) = [row_array(row_entry),col_array(col_entry)];
                    end
                end
                                       
                % Add new LZ to 3D LZs matrix
                LZs(:,:,numLZs + 1) = LZ;
                
            end
        end
        
    end
    
    toc()
    
    % Check if there were any possible contiguous regions for landing, exit
    % function if not.
    if numLZs == 0
        disp('No possible landing strips.')
        return
    end
    
    % Remove all zero entries from LZs matrix
    numLZs = numLZs + 1;
    LZs = LZs(:,:,1:numLZs);
    
    % Remove non-unique entries from LZs matrix
    [n,m,p]=size(LZs);
    a=reshape(LZs,n,[],1);
    b=reshape(a(:),n*m,[])';
    c=unique(b,'rows','stable')';
    LZs = reshape(c,n,m,[]);
    [n,m,p]=size(LZs);
    
    % Inititalize variables to use to populate array of LandingStrip
    % objects
    LzStrips(1,p) = LandingStrip();
    LonLat = zeros(n,2);
    DEM = zeros(n,1);
    TRI = zeros(n,1);
    LC = zeros(n,1);
        
    % Create array of LZ strips and calculate values
    for k = 1:p
        LzStrip = LZs(:,:,k);
        for coord = 1:length(LzStrip)
            LonLat(coord,:) = [N_new(LzStrip(coord,1),LzStrip(coord,2),5) N_new(LzStrip(coord,1),LzStrip(coord,2),6)];
            DEM(coord,:) = N_new(LzStrip(coord,1),LzStrip(coord,2),2);
            TRI(coord,:) = N_new(LzStrip(coord,1),LzStrip(coord,2),3);
            LC(coord,:) = N_new(LzStrip(coord,1),LzStrip(coord,2),4);
        end
        LzStrips(k) = LandingStrip(LzStrip, LonLat, DEM, TRI, LC, x_earth, y_earth, deg2met);
    end
       
    toc()
    
    % Determine the distance of each strip from all others and add to each
    % strip's distFromOtherStrips property.
%     distance = single(0);
%     for k = 1:p
%         for i = [1:k-1,k+1:p]
%             distance = LzStrips(k).computeDistance(LzStrips(i));
%             LzStrips(k).distFromOtherStrips = LzStrips(k).distFromOtherStrips + distance;
%         end
%     end
    
    toc()
    
    % ADD VALIDITY CHECKS HERE
    % Set validity to zero for strips with Land Cover type not 1, 2, or 3
    % (bare earth, vegetation, farm land), and/or with TRI > 2
    for k = 1:p
        %if LzStrips(k).avgTRI > 2 || LzStrips(k).avgDEM > Alt*0.5 || any([LzStrips(k).LC] == 6) || any([LzStrips(k).LC] == 1)
        if LzStrips(k).avgTRI > 1 || LzStrips(k).avgDEM > Alt*0.25 || any([LzStrips(k).LC] == 6) || any([LzStrips(k).LC] == 1)
            LzStrips(k).validStrip = 0;
        end
        LzStrips(k).score = LzStrips(k).computeScore();
    end
    
    toc()
    
    % Check if any valid landing strips remain, exits function if not
    numTotalValidLzStrips = nnz([LzStrips.validStrip]);
    if numTotalValidLzStrips == 0
        disp('No valid landing strips remain after pruning.');
        return
    end
        
    % Create a new array for the remaining valid LZs
    ValidLzStrips(1,numTotalValidLzStrips) = LandingStrip();
    numValidLzStrips = 0;
    for k = 1:p
        if LzStrips(k).validStrip
            numValidLzStrips = numValidLzStrips + 1;
            ValidLzStrips(numValidLzStrips) = LzStrips(k);
        end
    end
    
    % Find valid strip with maximum score and output parameters
    [maxScore, maxScoreIndex] = max([ValidLzStrips.score]);
    [minScore, minScoreIndex] = min([ValidLzStrips.score]);
    bestStrip = ValidLzStrips(maxScoreIndex);
    
    scoreGradient = [minScore:(maxScore - minScore)/3:maxScore];
    
    % Plot the resultant strips, valid strips, and chosen landing strip.
    figure(9)
    plot(-N_new(:,:,5),N_new(:,:,6),'s','MarkerSize',8,'MarkerFaceColor','red','MarkerEdgeColor','k')
    hold on
    for k = 1:p
        plot(-LzStrips(k).LonLat(:,1),LzStrips(k).LonLat(:,2),'s','MarkerSize',8,'MarkerFaceColor','blue','MarkerEdgeColor','k');
    end
    for k = 1:numValidLzStrips
        plot(-ValidLzStrips(k).LonLat(:,1),ValidLzStrips(k).LonLat(:,2),'s','MarkerSize',8,'MarkerFaceColor','green','MarkerEdgeColor','k');
    end
    plot(-bestStrip.LonLat(:,1),bestStrip.LonLat(:,2),'s','MarkerSize',8,'MarkerFaceColor','m','MarkerEdgeColor','k');
    plot(-x_earth, y_earth, '^','MarkerSize',8,'MarkerFaceColor','cyan','MarkerEdgeColor','red')
    xlabel('Longitude (deg)')
    ylabel('Latitude (deg)')
    title('Landing Strip Identification')
    
    h = zeros(5, 1);
    h(1) = plot(NaN,NaN,'^','MarkerSize',8,'MarkerFaceColor','c','MarkerEdgeColor','r');
    h(2) = plot(NaN,NaN,'s','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k');
    h(3) = plot(NaN,NaN,'s','MarkerSize',8,'MarkerFaceColor','b','MarkerEdgeColor','k');
    h(4) = plot(NaN,NaN,'s','MarkerSize',8,'MarkerFaceColor','g','MarkerEdgeColor','k');
    h(5) = plot(NaN,NaN,'s','MarkerSize',8,'MarkerFaceColor','m','MarkerEdgeColor','k');
    legend(h, 'Current UAV Position','High-Population','Low-Population; Invalid LZs','Low-Population; Valid LZs','Chosen Landing Strip', 'Location', 'bestoutside');
    
    hold off
    
    toc()
    
end