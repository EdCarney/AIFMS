function A = decode_data(loadfilename)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % decode_data.m
    %   This function loads saved compressed map data consisting of pixels
    %   that represent above- or below-threshold population values.
    %
    % Revision history
    %   6-26-18 Implemented by yuts, based on original concept by briggsmm
    %           and UMD requirements.
    %   7-3-18  Made into a function and renamed.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    gridsize = 10;
    strt=2;
    %% algorithm
    load(loadfilename);
    size_i = 1515; size_j = 1165;
    szM = size(M); j_offset = 0;
    crust_i_hi = mod(size_i-strt+1,gridsize); crust_j_hi = mod(size_j-strt+1,gridsize);
    for ii=1:szM(1)
       for jj=1:szM(2)
            i_start = (ii-1)*gridsize+strt; i_end = ii*gridsize+strt-1; j_start = (jj-1)*gridsize+strt; j_end = jj*gridsize+strt-1;
            if M(ii,jj) %binary
                A(i_start:i_end,j_start:j_end) = BM(1:gridsize,j_offset+1:j_offset+gridsize);
                j_offset = j_offset + gridsize;
            else %run-length encoding
                A(i_start:i_end,j_start:j_end) = R(ii,jj)*ones(gridsize,gridsize);
            end
        end
    end
    A(1:strt-1,1:size_j) = BM_crust_i_lo(1:strt-1,1:size_j);
    A(1:size_i,1:strt-1) = BM_crust_j_lo(1:size_i,1:strt-1);
    if crust_i_hi
        A(size_i-crust_i_hi+1:size_i,1:size_j) = BM_crust_i_hi(1:crust_i_hi,1:size_j);
    end
    if crust_j_hi
        A(1:size_i,size_j-crust_j_hi+1:size_j) = BM_crust_j_hi(1:size_i,1:crust_j_hi);
    end
    elapsed = toc;
end