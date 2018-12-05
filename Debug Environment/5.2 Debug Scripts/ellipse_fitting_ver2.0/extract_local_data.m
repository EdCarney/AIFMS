function [ map_local ] = extract_local_data( map, y_min, y_max, x_min, x_max )
%extract_data: Get heatmap data into array for user-input latitude and
%longitude range.
   ii_local = 1;
   for ii=y_min:y_max
       jj_local = 1;
       for jj=x_min:x_max
           map_local(ii_local,jj_local) = map(ii,jj);
           jj_local = jj_local+1;
       end
       ii_local = ii_local+1;
   end
end

