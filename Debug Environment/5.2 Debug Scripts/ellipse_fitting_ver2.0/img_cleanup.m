function [ img_processed ] = img_cleanup( img, spursize, colorval )
%img_cleanup: Fill singleton pixels (30.87meters) or spursize blocks with colorval
   sz = size(img);
   img_processed = img;
   
   for itn=1:10
   for idx=spursize:-1:1
       for ii=1+idx:sz(1)-idx
           for jj=1+idx:sz(2)-idx
               if ( img(ii,jj)~=colorval )
                   %fill in horizontal or vertical gaps
                   if ( (img(ii,jj+idx)==colorval && img(ii,jj-idx)==colorval) || (img(ii-idx,jj)==colorval && img(ii+idx,jj)==colorval))
                       img_processed(ii,jj) = colorval;
                   end
               end
           end
       end
       img = img_processed; %reset input image for next pass to output of current pass.
   end
   end
   
   spursize=1;
   %cut off diagonal connections, 1 iteration only
   for idx=spursize:-1:1
        for ii=1+idx:sz(1)-idx
           for jj=1+idx:sz(2)-idx
               if ( img(ii,jj)~=colorval )
                   %fill in diagonal gaps
                   if ( (img(ii+idx,jj+idx)==colorval && img(ii-idx,jj-idx)==colorval) || (img(ii-idx,jj+idx)==colorval && img(ii+idx,jj-idx)==colorval))
                       img_processed(ii,jj) = colorval;
                   end
               end
           end
        end
       img = img_processed;
   end
end

