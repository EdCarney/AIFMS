function [perimeter, numpixels] = find_perimeter( pixel_block, use_orig_dims, dont_use_window_edges )
   if ~exist('use_orig_dims','var')
       use_orig_dims = 0; %dark blue
   end
   if ~exist('dont_use_window_edges','var')
       dont_use_window_edges = 0;
   end
   debug = 0;
   hi_pop = 100;
   szp = size(pixel_block);
   perimeter = zeros(szp(1), szp(2));
   %enclose the pixel-block in a box
   szi = max(pixel_block, [], 1);
   img = hi_pop*ones(szi(1),szi(2));
   numnonzero = nnz(pixel_block(:,1));
   for blknum=1:numnonzero
      img(pixel_block(blknum,1), pixel_block(blknum,2)) = -9999;
   end
   
   ii=1; jj=1; px_idx=1; pxb_idx=1;
   while (ii && jj) && (px_idx <= szp(1))
      ii = pixel_block(px_idx, 1);
      jj = pixel_block(px_idx, 2);
      if (ii && jj)
         if dont_use_window_edges
            crit_1 = 0;
         else
            crit_1 = ii==1 || jj==1 || ii==szi(1) || jj==szi(2);
         end
         crit_2 = img(ii,min(jj+1, szi(2))) == hi_pop || img(ii,max(jj-1, 1)) == hi_pop || img(max(ii-1, 1),jj) == hi_pop || img(min(ii+1, szi(1)),jj) == hi_pop;
         if (crit_1 || crit_2)
            perimeter(pxb_idx, 1) = ii;
            perimeter(pxb_idx, 2) = jj;
            pxb_idx = pxb_idx + 1;
            img(ii,jj) = 20;
         end
      end
      px_idx = px_idx + 1;
   end
   numpixels = [px_idx-1 pxb_idx-1];
   if ~use_orig_dims
      perimeter = perimeter(1:numpixels(2),:);
   end
   if debug
      figure; %DEBUG ONLY
      image(img); %DEBUG ONLY
   end
end
