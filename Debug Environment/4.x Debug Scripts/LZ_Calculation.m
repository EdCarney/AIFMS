clear LZs LZs_xy rows cols i j newr newc

SF = 0.003;

[rows,cols] = size(N_new(:,:,1));

LZs = zeros(rows,cols);

for i=1:rows

	for j=1:cols
    
		try
            
            %N_new(i-1:i+3,j-1:j+1,1)<= SF
            
            %all(all((N_new(i-1:i+3,j-1:j+1,1)<= SF)))
            
			% check if any consideration along the row is true
		
			if all(all((N_new(i-1:i+3,j-1:j+1,1)<= SF)))
                LZs(i-1:i+3,j-1:j+1) = N_new(i-1:i+3,j-1:j+1,1);
				disp('CASE 1');
            elseif all(all((N_new(i-2:i+2,j-1:j+1,1)<= SF)))
                LZs(i-2:i+2,j-1:j+1) = N_new(i-2:i+2,j-1:j+1,1);
				disp('CASE 2');
            elseif all(all((N_new(i-3:i+1,j-1:j+1,1)<= SF)))
                LZs(i-3:i+1,j-1:j+1) = N_new(i-3:i+1,j-1:j+1,1);
				disp('CASE 3');
                
			% check if any consideration along the column is true
			
            elseif all(all((N_new(i-1:i+1,j-1:j+3,1)<= SF)))
                LZs(i-1:i+1,j-1:j+3) = N_new(i-1:i+1,j-1:j+3,1);
				disp('CASE 4');
            elseif all(all((N_new(i-1:i+1,j-2:j+2,1)<= SF)))
                LZs(i-1:i+1,j-2:j+2) = N_new(i-1:i+1,j-2:j+2,1);
				disp('CASE 5');
            elseif all(all((N_new(i-1:i+1,j-3:j+1,1)<= SF)))
                LZs(i-1:i+1,j-3:j+1) = N_new(i-1:i+1,j-3:j+1,1);
				disp('CASE 6');
				
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
end

I = find(LZs);
[newr,newc] = ind2sub(size(N_new(:,:,1)),I);
for i=1:numel(newr)
    LZs_xy(i,:) = [N_new(newr(i),newc(i),2),N_new(newr(i),newc(i),3)];
end