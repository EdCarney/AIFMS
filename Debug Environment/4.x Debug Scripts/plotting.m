clear lats longs
figure(1)
plot(N_new(:,:,2),N_new(:,:,3),'g.')
hold on
%plot(FGIF(:,1),FGIF(:,2),'.b')
for i = 1:length(SRM(:,1))
lats(i) = N_new(SRM(i,1),SRM(i,2),2);
longs(i) = N_new(SRM(i,1),SRM(i,2),3);
end
longs = longs';
lats = lats';
plot(lats(:,1),longs(:,1),'rx')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SF = 0.0035;

[rows,cols] = size(N_new(:,:,1));

n = 1;

casevar = 0;

for i=1:rows

	for j=1:cols
		
		try
		
			% check if any consideration along the row is true
		
			if all(all(N_new(i-1:i+3,j-1:j+1,1))) <= SF
			  LZs(n) = N_new(i-1:i+3,j-1:j+1,1);
				n = n + 1;
				disp('CASE 1');
			elseif all(all(N_new(i-2:i+2,j-1:j+1,1))) <= SF
				  LZs(n) = N_new(i-2:i+2,j-1:j+1,1);
				n = n + 1;
				disp('CASE 2');
			elseif all(all(N_new(i-3:i+1,j-1:j+1,1))) <= SF
				  LZs(n) = N_new(i-3:i+1,j-1:j+1,1);
				n = n + 1;
				disp('CASE 3');
				
			% check if any consideration along the column is true
			
			elseif all(all(N_new(i-1:i+1,j-1:j+3,1))) <= SF
				  LZs(n) = N_new(i-1:i+1,j-1:j+3,1);
				n = n + 1;
				disp('CASE 4');
			elseif all(all(N_new(i-1:i+1,j-2:j+2,1))) <= SF
				  LZs(n) = N_new(i-1:i+1,j-2:j+2,1);
				n = n + 1;
				disp('CASE 5');
			elseif all(all(N_new(i-1:i+1,j-3:j+1,1))) <= SF
				  LZs(n) = N_new(i-1:i+1,j-3:j+1,1);
				n = n + 1;
				disp('CASE 6');
				
			% continue if no LZ possibilities
				
			else
				continue
				
			end
				
		catch me
			if strcmp(me.identifier,'MATLAB:badsubscript');
				continue
			else
				disp('BAD STUFF ERROR!');
				disp(me.identifier);
				return
			end
			
		end
	end
end