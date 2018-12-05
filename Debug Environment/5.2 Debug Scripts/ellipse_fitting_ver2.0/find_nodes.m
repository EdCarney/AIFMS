function [ nodes, endpts, A_nodes ] = find_nodes( A, foregrnd_col, backgrnd_col, node_col, endpt_col )
%find_nodes: Find intersection points in a pixel-wide skeletonized image.

    nodes=[]; endpts=[]; A_nodes=[];
    f = foregrnd_col;
    b = backgrnd_col;
    pattern_1a = [b f -1; f f -1; -1 -1 -1];
    pattern_1b = [-1 f b; -1 f f; -1 -1 -1];
    pattern_1c = [-1 -1 -1; -1 f f; -1 f b];
    pattern_1d = [-1 -1 -1; f f -1; b f -1];
    
    pattern_2a = [f b f; -1 f -1; -1 f -1];
    pattern_2b = [-1 -1 f; f f b; -1 -1 f];
    pattern_2c = [-1 f -1; -1 f -1; f b f];
    pattern_2d = [f -1 -1; b f f; f -1 -1];

    pattern_3a = [b b b; b f b; b f b];
    pattern_3b = [b b b; f f b; b b b];
    pattern_3c = [b f b; b f b; b b b];
    pattern_3d = [b b b; b f f; b b b];
    
    pattern_4a = [b b b; b f b; b b f];
    pattern_4b = [b b b; b f b; f b b];
    pattern_4c = [f b b; b f b; b b b];
    pattern_4d = [b b f; b f b; b b b];
    
    pattern_1 = [pattern_1a pattern_1b pattern_1c pattern_1d]; %90deg angles
    pattern_2 = [pattern_2a pattern_2b pattern_2c pattern_2d]; %acute angles...y shape
    pattern_3 = [pattern_3a pattern_3b pattern_3c pattern_3d]; %endpoint type 1
    pattern_4 = [pattern_4a pattern_4b pattern_4c pattern_4d]; %endpoint type 2
    
    imax = size(A,1);
    jmax = size(A,2);
    A_nodes = A; %initialization
    idx = 1;
    idx2 = 1;
   
    for rot=1:4
        arr_lo = (rot-1)*3 + 1;
        arr_hi = (rot-1)*3 + 3;
        for se=1:4
            if se==1
                se_array = pattern_1(:,arr_lo:arr_hi);
            elseif se==2
                se_array = pattern_2(:,arr_lo:arr_hi);
            elseif se==3
                se_array = pattern_3(:,arr_lo:arr_hi);
            else
                se_array = pattern_4(:,arr_lo:arr_hi);
            end
            for ii=2:imax-1 %apply structuring element to entire image.
                for jj=2:jmax-1
                    if se<=2
                        node_pixel = 1;
                        endpt_pixel = 0;
                    else
                        endpt_pixel = 1;
                        node_pixel = 0;
                    end
                    for se_i=1:3 %compare with 'structuring element'
                        offset_ii = se_i - 2;
                        for se_j=1:3
                            offset_jj = se_j - 2;
                            if ((se_array(se_i,se_j) ~= A(ii+offset_ii,jj+offset_jj)) && (se_array(se_i,se_j) ~= -1))
                                if se<=2
                                    node_pixel = 0; %if any mismatch, pixel is not a node.
                                else
                                    endpt_pixel = 0; %if any mismatch, pixel is not an endpoint.
                                end
                            end
                        end
                    end
                    if node_pixel == 1
                        %check for duplicates
                        numnodes = idx-1; %max(size(nodes));
                        dup = 0;
                        for idx_n=1:numnodes
                            if (ii==nodes(idx_n,1) && jj==nodes(idx_n,2))
                                dup = 1;
                            end
                            distvec = [ii jj] - [nodes(idx_n,1) nodes(idx_n,2)];
                            dist_lim = 5;
                            if norm(distvec, 2)<dist_lim
                                dup = 1;
                            end
                        end
                        if dup==0
                            A_nodes(ii,jj) = node_col;
                            nodes(idx,:) = [ii jj];
                            idx = idx + 1;
                        end
                    end
                    if endpt_pixel == 1
                        A_nodes(ii,jj) = endpt_col;
                        endpts(idx2,:) = [ii jj];
                        idx2 = idx2 + 1;
                    end
                end
            end
        end
    end
    if isempty(endpts) == 1
        img_thin = A;
        img_thin(img_thin == -9999) = 1;
        img_thin(img_thin == 100) = 0;
        endpts_mat = bwmorph(img_thin, 'endpoints');
        [xendpts yendpts] = find(endpts_mat == 1);
        endpts = [yendpts(:) xendpts(:)];
    end
end

