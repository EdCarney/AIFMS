clear LZs LZs_xy LZs_all_plot LZs_horz LZs_horiz LZs_horiz_plot LZs_vert LZs_vert_plot LZs_safe_pos LZs_safe_horiz_pos LZs_safe_vert_pos LZs_safe_horiz LZs_safe_vert


SF = 0.003;
n = 1;
len = length(SRM);

for k=1:len

    i = SRM(k,1);
    j = SRM(k,2);
    
		try
            
			% check if any consideration along the row is true
		
			if all(all((N_new(i-1:i+3,j-1:j+1,1)<= SF)))
                LZs_vert(k,:) = [i-1:i+3,j-1:j+1];
				%disp('CASE 1');
            elseif all(all((N_new(i-2:i+2,j-1:j+1,1)<= SF)))
                LZs_vert(k,:) = [i-2:i+2,j-1:j+1];
				%disp('CASE 2');
            elseif all(all((N_new(i-3:i+1,j-1:j+1,1)<= SF)))
                LZs_vert(k,:) = [i-3:i+1,j-1:j+1];
				%disp('CASE 3');
                
			% check if any consideration along the column is true
			
            elseif all(all((N_new(i-1:i+1,j-1:j+3,1)<= SF)))
                LZs_horiz(k,:) = [i-1:i+1,j-1:j+3];
				%disp('CASE 4');
            elseif all(all((N_new(i-1:i+1,j-2:j+2,1)<= SF)))
                LZs_horiz(k,:) = [i-1:i+1,j-2:j+2];
				%disp('CASE 5');
            elseif all(all((N_new(i-1:i+1,j-3:j+1,1)<= SF)))
                LZs_horiz(k,:) = [i-1:i+1,j-3:j+1];
				%disp('CASE 6');
				
			% continue if no LZ possibilities
            
            else
                continue
            end
				
		catch me
			if strcmp(me.identifier,'MATLAB:badsubscript') | strcmp(me.identifier,'MATLAB:subsassigndimmismatch');
				continue
			else
				disp('BAD STUFF ERROR!');
				disp(me.identifier);
				return
			end
			
        end
end

% Remove zero rows

LZs_vert(~any(LZs_vert,2),:) = [];
LZs_horiz(~any(LZs_horiz,2),:) = [];

% Use row and column values to get all 15 points for each vertical
% and horizontal LZ

for k = 1:length(LZs_vert)
    for i = 1:5
        for j = 1:3
           LZs_vert_plot((j+(i-1)*3)+15*(k-1),:) = [LZs_vert(k,i),LZs_vert(k,5+j)];
        end
    end
end

for k = 1:length(LZs_horiz)
    for i = 1:3
        for j = 1:5
            LZs_horiz_plot((j+(i-1)*5)+15*(k-1),:) = [LZs_horiz(k,i),LZs_horiz(k,3+j)];
        end
    end
end

% Concatenate arrays and remove duplicate values

LZs_all_plot = [LZs_vert_plot;LZs_horiz_plot];
LZs_all_plot = unique(LZs_all_plot,'rows');

for i = 1:length(SRM)
    SRM_longs(i) = N_new(SRM(i,1),SRM(i,2),2);
    SRM_lats(i) = N_new(SRM(i,1),SRM(i,2),3);
end

for i = 1:length(LZs_all_plot)
    LZs_all_plot_longs(i) = N_new(LZs_all_plot(i,1),LZs_all_plot(i,2),2);
    LZs_all_plot_lats(i) = N_new(LZs_all_plot(i,1),LZs_all_plot(i,2),3);
end

% Check if the potential landing strips are safe for landing; here "safe"
% is quantified as an 30mx30m square whose area does not have a standard
% deviation in elevation greater than 2 meters.

LZs_all_plot_longs = -1*LZs_all_plot_longs';
LZs_all_plot_lats = LZs_all_plot_lats';

maxStdDev = 6; % m

for i = 1:length(LZs_all_plot_lats)
    [minLatDistance, indexOfLatMin] = min(abs(dLats-LZs_all_plot_lats(i)));
    [minLongDistance, indexOfLongMin] = min(abs(dLongs-LZs_all_plot_longs(i)));
    if std_all(indexOfLatMin,indexOfLongMin) > maxStdDev
        continue
    else
        safeLZs(i,1) = LZs_all_plot_lats(i);
        safeLZs(i,2) = LZs_all_plot_longs(i);
    end
end

safeLZs = unique(safeLZs,'rows');
safeLZs(1,:) = [];

% Find if remaining safe strips provide a contiguous landing region. Use
% N_new matrix and safeLZ lat/longs to find unique landing strips.

longVec = -N_new(1,:,2)';
latVec = N_new(:,1,3);

findMatrix = zeros(numel(N_new(:,1,1)));
findLZs = zeros(numel(safeLZs(:,1)),2); 

for i = 1:numel(safeLZs(:,1))
    findLat = find(latVec == safeLZs(i,1));
    findLong = find(longVec == safeLZs(i,2));
    findLZs(i,:) = [findLat,findLong];
    findMatrix(findLong,findLat) = 1;
end


for k=1:length(safeLZs)

    j = findLZs(k,1);
    i = findLZs(k,2);
    
		try
            
			% check if any consideration along the column is true
		
			if all(all((findMatrix(i-1:i+3,j-1:j+1,1))))
                LZs_safe_horiz(k,:) = [i-1:i+3,j-1:j+1];
				%disp('CASE 1');
            elseif all(all((findMatrix(i-2:i+2,j-1:j+1,1))))
                LZs_safe_horiz(k,:) = [i-2:i+2,j-1:j+1];
				%disp('CASE 2');
            elseif all(all((findMatrix(i-3:i+1,j-1:j+1,1))))
                LZs_safe_horiz(k,:) = [i-3:i+1,j-1:j+1];
				%disp('CASE 3');
                
			% check if any consideration along the row is true
			
            elseif all(all((findMatrix(i-1:i+1,j-1:j+3,1))))
                LZs_safe_vert(k,:) = [i-1:i+1,j-1:j+3];
				%disp('CASE 4');
            elseif all(all((findMatrix(i-1:i+1,j-2:j+2,1))))
                LZs_safe_vert(k,:) = [i-1:i+1,j-2:j+2];
				%disp('CASE 5');
            elseif all(all((findMatrix(i-1:i+1,j-3:j+1,1))))
                LZs_safe_vert(k,:) = [i-1:i+1,j-3:j+1];
				%disp('CASE 6');
				
			% continue if no LZ possibilities
            
            else
                continue
            end
				
		catch me
			if strcmp(me.identifier,'MATLAB:badsubscript') | strcmp(me.identifier,'MATLAB:subsassigndimmismatch');
				continue
			else
				disp('BAD STUFF ERROR!');
				disp(me.identifier);
				return
			end
			
        end
end

% Remove zero and duplicate rows; calculate the posiiton of the centriod of
% each landing strip to determine which is closest.

if exist('LZs_safe_vert','var')
    LZs_safe_vert(~any(LZs_safe_vert,2),:) = [];
    LZs_safe_vert = unique(LZs_safe_vert,'rows');
    for k = 1:length(LZs_safe_vert(:,1))
        LZs_safe_vert_pos(k,:) = [mean(N_new(LZs_safe_vert(k,4:8),1,3)),-1*mean(N_new(1,LZs_safe_vert(k,1:3),2))];
        LZs_safe_vert_coords(k,:) = [N_new(LZs_safe_vert(k,4:8),1,3)',N_new(1,LZs_safe_vert(k,1:3),2)];
    end
end
if exist('LZs_safe_horiz','var')
    LZs_safe_horiz(~any(LZs_safe_horiz,2),:) = [];
    LZs_safe_horiz = unique(LZs_safe_horiz,'rows');
    for k = 1:length(LZs_safe_horiz(:,1))
        LZs_safe_horiz_pos(k,:) = [mean(N_new(LZs_safe_horiz(k,6:8),1,3)),-1*mean(N_new(1,LZs_safe_horiz(k,1:5),2))];
        LZs_safe_horiz_coords(k,:) = [N_new(LZs_safe_horiz(k,6:8),1,3)',N_new(1,LZs_safe_horiz(k,1:5),2)];
    end
end

% Check if there are safe LZs in both or just one direction and use to
% create one matrix of safe LZs.

if exist('LZs_safe_vert','var') && exist('LZs_safe_horiz','var')
    LZs_safe_pos = [LZs_safe_vert_pos;LZs_safe_horiz_pos];
    LZs_safe_coords = [LZs_safe_vert_coords;LZs_safe_horiz_coords];
    LZs_safe = [LZs_safe_vert;LZs_safe_horiz];
elseif exist('LZs_safe_vert','var')
    LZs_safe_pos = LZs_safe_vert_pos;
    LZs_safe_coords = [LZs_safe_vert_coords];
    LZs_safe = [LZs_safe_vert];
elseif exist('LZs_safe_horiz','var')
    LZs_safe_pos = LZs_safe_horiz_pos;
    LZs_safe_coords = [LZs_safe_horiz_coords];
    LZs_safe = [LZs_safe_horiz];
end
    
% Check which safe LZ is closest to the current AC position.

x_earth = -76.94695; % Simulated AC x posititon
y_earth = 38.99479; % Simulated AC y posititon.

for i = 1:length(LZs_safe_pos)
    dist(i) = sqrt((LZs_safe_pos(i,1)-y_earth)^2 + (LZs_safe_pos(i,2)-x_earth)^2);
end

[minDist,I] = min(dist);

% Get resultant LZ points for plottng purposes. Will check if solution is
% vertical or horizontal.

if LZs_safe(I,5)-LZs_safe(I,1) > 4
    for i = 1:5
        for j = 1:3
            responseCoords((i-1)*3+j,:) = [LZs_safe_coords(I,i),LZs_safe_coords(I,5+j)];
        end
    end
else
    for i = 1:3
        for j = 1:5
            responseCoords((i-1)*5+j,:) = [LZs_safe_coords(I,i),LZs_safe_coords(I,3+j)];
        end
    end
end

response = LZs_safe_pos(I,:)

figure(1)
plot(-SRM_longs,SRM_lats,'rx')
hold on
plot(LZs_all_plot_longs,LZs_all_plot_lats,'gx')
plot(safeLZs(:,2),safeLZs(:,1),'bx')
plot(-responseCoords(:,2),responseCoords(:,1),'sm')
xlabel('Latitude')
ylabel('Longitude')
legend('PD < 1','Continuous Region - PD < 1','PD < 1 and ESD < 6','Selcted Landing Zone')