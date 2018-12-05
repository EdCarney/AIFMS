function [ A_thin ] = thin_bin_img( A, foregrndcol, backgrndcol )
%thin_img: Perform morphological thinning of a binary image

   f = foregrndcol; b = backgrndcol;
   struct_elem_1_rot0 = [b b b; -1 f -1; f f f];
   struct_elem_2_rot0 = [-1 b b; f f b; -1 f -1];
   
   struct_elem_1_rot90 = [f -1 b; f f b; f -1 b];
   struct_elem_2_rot90 = [-1 f -1; f f b; -1 b b];
   
   struct_elem_1_rot180 = [f f f; -1 f -1; b b b];
   struct_elem_2_rot180 = [-1 f -1; b f f; b b -1];
   
   struct_elem_1_rot270 = [b -1 f; b f f; b -1 f];
   struct_elem_2_rot270 = [b b -1; b f f; -1 f -1];
   
   struct_elem_1 = [struct_elem_1_rot0 struct_elem_1_rot90 struct_elem_1_rot180 struct_elem_1_rot270]; % 3x12 array
   struct_elem_2 = [struct_elem_2_rot0 struct_elem_2_rot90 struct_elem_2_rot180 struct_elem_2_rot270]; % 3x12 array
   
   imax = size(A,1);
   jmax = size(A,2);
   A_thin = A; %initialization
   
   for itn=1:30
       for rot=1:4
           arr_lo = (rot-1)*3 + 1;
           arr_hi = (rot-1)*3 + 3;
           for se=1:2 %apply thinning using each structuring element, for a particular 90deg rotation.
               if se==1
                   se_array = struct_elem_1(:,arr_lo:arr_hi);
               else
                   se_array = struct_elem_2(:,arr_lo:arr_hi);
               end               
               for ii=2:imax-1 %apply structuring element to entire image.
                   for jj=2:jmax-1               
                       remove_pixel = 1;
                       for se_i=1:3 %compare with 'structuring element'
                           offset_ii = se_i - 2;
                           for se_j=1:3
                               offset_jj = se_j - 2;
                               if ((se_array(se_i,se_j) ~= A_thin(ii+offset_ii,jj+offset_jj)) && (se_array(se_i,se_j) ~= -1))
                                   remove_pixel = 0; %if mismatch, don't modify pixel.
                               end
                           end
                       end
                       if remove_pixel == 1
                           A_thin(ii,jj) = backgrndcol; %if applicable, change foreground pixel to background pixel.
                       end
                   end
               end
           end
       end
   end
end

