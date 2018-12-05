function endpts = gen_endpts(FGIF_p, blob, img_45m, hi_pop)
endpts_tmp = [round(FGIF_p(1,2),0) round(FGIF_p(1,1),0); round(FGIF_p(end,2),0) round(FGIF_p(end,1),0)];
endpts_tmp(~endpts_tmp) = 1;
% Choose endpts within given blob
check = ismember(endpts_tmp, blob, 'rows');
cnt = 0;
cnt1 = 0;
cnt2 = 0;
cnt3 = 0;
crit_flag = 0;
while check(1) == 0 || check(2) == 0
    check = ismember(endpts_tmp, blob, 'rows');
    cnt = cnt + 1;
    if cnt >= (size(FGIF_p,1) - 1)
        endpts_tmp = [];
        crit_flag = 1;
        disp('broke gen_endpts loop 1')
        break
    end
    if check(1) == 0 && check(2) == 0
        cnt1 = cnt1 + 1;
        endpts_tmp = [round(FGIF_p(1+cnt1,2),0) round(FGIF_p(1+cnt1,1),0); round(FGIF_p(end-cnt1,2),0) round(FGIF_p(end-cnt1,1),0)];
        endpts_tmp(~endpts_tmp) = 1;
    elseif check(1) == 0 && check(2) == 1
        cnt2 = cnt2 + 1;
        endpts_tmp = [round(FGIF_p(1+cnt1+cnt2,2),0) round(FGIF_p(1+cnt1+cnt2,1),0); endpts_tmp(2,1) endpts_tmp(2,2)];
        endpts_tmp(~endpts_tmp) = 1;
    elseif check(1) == 1 && check(2) == 0
        cnt3 = cnt3 + 1;
        endpts_tmp = [endpts_tmp(1,1) endpts_tmp(1,2); round(FGIF_p(end-cnt1-cnt3,2),0) round(FGIF_p(end-cnt1-cnt3,1),0)];
        endpts_tmp(~endpts_tmp) = 1;
    end
end
if crit_flag == 1
    endpts = [];
    return
end

FGIF_round = round(FGIF_p(:,:), 0);
FGIF_round(~FGIF_round) = 1;
FGIF_round(:,[1 2]) = FGIF_round(:, [2 1]);
[~,~,idx_min] = intersect(endpts_tmp(1,:), FGIF_round(:,:), 'rows');
idx_min = idx_min(1);
[~,~,idx_max] = intersect(endpts_tmp(2,:), FGIF_round(:,:), 'rows');
idx_max = idx_max(1);
%Get list of FGIF pts within blob
endpts_inblob = [];
for i = 1:(idx_max - idx_min + 1)
    idx = idx_min + i - 1;
    endpts_add = [round(FGIF_p(idx,2),0) round(FGIF_p(idx,1),0)];
    endpts_add(~endpts_add) = 1;
    endpts_inblob = [endpts_inblob; endpts_add];
end
endpts_inblob = unique(endpts_inblob, 'rows');

% Determine if endpoints intersect any high pop
pts_bad = [];
for idx_pts = 1:size(endpts_inblob,1)
    if img_45m(endpts_inblob(idx_pts,1), endpts_inblob(idx_pts,2)) == hi_pop
        pts_bad = [pts_bad; idx_pts];
    else 
        endpts = endpts_inblob(:,:);
    end
end

iter = 0;
endpts_pot = [endpts_inblob(1,:); endpts_inblob(end,:)];
while isempty(pts_bad) == 0
    pts_bad;
    iter = iter + 1;
    if iter >= (size(endpts_inblob,1) - 1)
        endpts = [];
        disp('broke gen_endpts loop 2')
        break
    end
    % Choose different endpts to check
    val = zeros(size(pts_bad,1)+1,1);
    val(1) = pts_bad(1) - 1; 
    if size(pts_bad,1) > 1
        for i = 2:size(pts_bad,1)
            val(i) = pts_bad(i) - pts_bad(i-1); 
        end
        val(size(pts_bad,1)+1,1) = size(endpts_inblob,1) - pts_bad(end);
        [~,idx_max] = max(val(:));
        if idx_max == 1 
            idx1 = 0;
            idx2 = pts_bad(1);
        elseif idx_max == size(val,1)
            idx1 = pts_bad(idx_max - 1);
            idx2 = size(endpts_inblob,1) + 1;
        else
            idx1 = pts_bad(idx_max - 1); 
            idx2 = pts_bad(idx_max); 
        end
    else
        val(2) = size(endpts_inblob,1) - pts_bad(1);
        [~,idx_max] = max(val(:));
        if idx_max == 1
            idx1 = 0
            idx2 = pts_bad(1);
        elseif idx_max == 2 
            idx1 = pts_bad(1)
            idx2 = size(endpts_inblob,1);
        end
    end    
    endpts_pot = endpts_inblob((idx1+1):(idx2-1),:); % need to modify to account for if first or last endpts are used 
    
    pts_bad = [];
    for idx_pts = 1:size(endpts_pot,1)
        if img_45m(endpts_pot(idx_pts,1), endpts_pot(idx_pts,2)) == hi_pop
            pts_bad = [pts_bad; idx_pts];
        else 
            endpts = endpts_pot;
        end
    end
end   
end