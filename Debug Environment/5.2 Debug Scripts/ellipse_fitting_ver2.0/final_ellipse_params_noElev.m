function [ra_, rb_, centroid_, ang_, area_, high_low_rat_, dist_, xnew_, ynew_] = final_ellipse_params_noElev(ra, rb, centroid, ang, img, x_earth, y_earth, imin, jmin, LS_USA, deg2met, FGIF)
IM = img < 150;
[y_bel, x_bel] = find(IM);
[y_abv, x_abv] = find(~IM);

if isempty(ra) == 0
    ra_ = ra;
    rb_ = rb;
    centroid_ = centroid(:,:);
    x_centroid = centroid(:,2);
    y_centroid = centroid(:,1);
    ang_ = ang;
    dist_ = [];

    centroids = centroid'; % y,x
    centroids([1 2],:) = centroids([2 1],:);
    P_in_bel = zeros(2, length(x_bel), length(ra));
    P_in_abv = zeros(2, length(x_bel), length(ra));
 
    for k = 1:length(ra)
        ab = [ra(k) rb(k)];
        r(k) = wrapTo2Pi(ang(k));
        R = [cos(r(k)) -sin(r(k)); sin(r(k)) cos(r(k))];
        P_bel = [x_bel'; y_bel'];
        P_abv = [x_abv'; y_abv'];
        dP_bel = P_bel - centroids(:,k);
        dP_bel = R'*dP_bel; 
        dP_bel = diag(1./ab)*dP_bel;
        dP_abv = P_abv - centroids(:,k);
        dP_abv = R'*dP_abv; 
        dP_abv = diag(1./ab)*dP_abv;
        id_in_bel = sum(dP_bel.^2, 1) <= 1;
        id_in_abv = sum(dP_abv.^2, 1) <= 1;
        X = sum(id_in_bel(:) == 1);
        Y = sum(id_in_abv(:) == 1);

        X = sum(id_in_bel(:) == 1);
        Y = sum(id_in_abv(:) == 1);
        P_in_bel(1:2,1:X,k) = P_bel(:, id_in_bel);
        P_in_abv(1:2,1:Y,k) = P_abv(:, id_in_abv);
        
        [~,idx_FGIF,~] = intersect(FGIF, P_in_bel(:,:,k)', 'rows');
        if isempty(idx_FGIF) == 0
            Num_below(k) = nnz(P_in_bel(:,:,k))/2; % number of cells under pop thresh
            Num_above(k) = nnz(P_in_abv(:,:,k))/2;
            high_low_rat_(k) = Num_above(k)/(Num_above(k) + Num_below(k));
            area_(k) = pi*ra_(k)*rb_(k);   
            yrow = FGIF(idx_FGIF,2) + min(imin) - 1;
            xcol = FGIF(idx_FGIF,1) + min(jmin) - 1;
            ME2 = LS_USA(:,:,2);
            ME3 = LS_USA(:,:,3);
            xnew1 = ME2(sub2ind(size(ME2), yrow, xcol));
            xnew_(k) = xnew1(1);
            ynew1 = ME3(sub2ind(size(ME3), yrow, xcol));
            ynew_(k) = ynew1(1);
            dist_(k) = sqrt((x_earth - xnew_(k))^2 + (y_earth - ynew_(k))^2)*deg2met;
        else
            high_low_rat_(k) = 0;
            area_(k) = 0;
            xnew_(k) = 0;
            ynew_(k) = 0;
            dist_(k) = 0;
        end
    end
    % Extract pts inside each ellipse and check if ellipses are within
    % others
    for k = 1:size(P_in_bel, 3)
        var{k} = P_in_bel(:,:,k)';
    end
    pts_in = [];
    if numel(var) > 1
        for i = 1:numel(var)
            pts_in_tmp = var{i};
            pts_in = [pts_in; var{i}];
        end
    end

    row_del = [];
    if numel(var) > 1
        for i = 1:numel(var)
            var1_tmp = var{i};
            var1_tmp = var1_tmp(any(var1_tmp, 2), :);
            var2_tmp = var;
            var2_tmp{i} = [];
            var2_tmp(~cellfun('isempty', var2_tmp));
            var3_tmp = cat(1, var2_tmp{:});
            var3_tmp = var3_tmp(any(var3_tmp, 2), :);
            [~,idx_var1, ~] = intersect(var1_tmp, var3_tmp, 'rows');
            if size(idx_var1, 1) == size(var1_tmp, 1)
                eqlcheck = [];
                P_in_bel_del = P_in_bel(:,:,:);
                P_in_bel_del(:,:,i) = [];
                for j = 1:size(P_in_bel_del, 3)
                    pts_non0 = P_in_bel_del(:,:,j)';
                    pts_non0(~any(pts_non0,2),:) = [];
                    eql = isequal(var1_tmp, pts_non0);
                    eqlcheck = [eqlcheck; eql];
                end
                check = find(double(eqlcheck) == 1);
                if isempty(check) == 1              
                    row_del = [row_del; i];
                end
            elseif xnew_(1,i) == 0 && ynew_(1,i) == 0
                row_del = [row_del; i];
            end
        end
    end

    pts_in_del = [];
    if size(row_del,1) > 1
        for i = 1:size(row_del,1)
            pts_tmp = P_in_bel(:,:,row_del(i))';
            pts_in_del = [pts_in_del; pts_tmp];
        end
        ra_([row_del]) = [];
        rb_([row_del]) = [];
        ang_([row_del]) = [];
        centroid_([row_del], :) = [];
        area_([row_del]) = [];
        high_low_rat_([row_del]) = [];
        P_in_bel(:,:,[row_del]) = [];
        P_in_abv(:,:,[row_del]) = [];
        dist_([row_del]) = [];
        xnew_([row_del]) = [];
        ynew_([row_del]) = [];
        disp(['Removed ' num2str(size(row_del,1))  ' Ellipses'])
    end      
else
    ra_ = [];
    rb_ = [];
    centroid_ = [];
    ang_ = [];
    xnew_ = [];
    ynew_ = [];
    dist_ = [];
    area_ = [];
    high_low_rat_ = [];
end
end