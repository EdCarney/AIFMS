function [ img_array_thresh ] = img_thresh( img_array, threshold_lo, threshold_hi, flip_colors )
   if ~exist('threshold_hi','var')
       threshold_hi = 9999;
   end
   if ~exist('flip_colors','var')
       flip_colors = 0;
   end
   sz = size(img_array);
   sz_1d = prod(sz); %in case of a many-dimensional array
   tmp = reshape(img_array, [sz_1d 1]);
   hi_val = 100;
   lo_val = -9999;
   for idx=1:sz_1d
       if (tmp(idx) > threshold_lo && tmp(idx) < threshold_hi)
           if ~flip_colors
               tmp(idx) = hi_val;
           else
               tmp(idx) = lo_val;
           end
       else
           if ~flip_colors
               tmp(idx) = lo_val;
           else
               tmp(idx) = hi_val;
           end
       end
   end
   img_array_thresh = reshape(tmp, sz);
end
         