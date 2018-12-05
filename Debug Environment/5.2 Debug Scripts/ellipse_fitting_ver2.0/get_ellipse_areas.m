function [area_ high_low_rat_] = get_ellipse_areas(ra_, rb_, centroid_, ang_, img)
IM = img < 100;
[y_bel, x_bel] = find(IM);
[y_abv, x_abv] = find(~IM);

centroids = centroid_'; % y,x
centroids([1 2],:) = centroids([2 1],:);
for k = 1:length(ra_)
    ab = [ra_(k) rb_(k)];
    r(k) = wrapTo2Pi(ang_(k));
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
    
    P_in_bel = P_bel(:, id_in_bel); 
    P_in_abv = P_abv(:, id_in_abv);
    
    Num_below(k) = length(P_in_bel); % number of cells under pop thresh
    Num_above(k) = length(P_in_abv);
    high_low_rat_(k) = Num_above(k)/(Num_above(k) + Num_below(k));
    area_(k) = pi*ra_(k)*rb_(k);
    
    % plotting 
    plot(P_in_bel(1,:), P_in_bel(2,:), '.g', 'MarkerSize', 10, 'Linewidth', 2);
    plot(P_in_abv(1,:), P_in_abv(2,:), '.r', 'MarkerSize', 10, 'Linewidth', 2);
end
end