function endpts = gen_endpts_noElev(centroid_opts, blob, img_45m, hi_pop, psi0)
grow_len = 2;
ang = psi0 + centroid_opts(1,3);
angle = [ang; ang];
endpts_y = [centroid_opts(1,1) - grow_len*sin(ang); centroid_opts(1,1) + grow_len*sin(ang)];
endpts_x = [centroid_opts(1,2) - grow_len*cos(ang); centroid_opts(1,2) + grow_len*cos(ang)];   
endpts_tmp = [round(endpts_y,0) round(endpts_x,0) angle];
endpts_tmp = endpts_tmp(:,1:2);
% Choose endpts within given blob
check = ismember(endpts_tmp, blob, 'rows');

if check(1) == 0 && check(2) == 0
    endpts = [];
    return
else
    if check(1) == 0 && check(2) == 1
        endpts_y2 = [centroid_opts(1,1); centroid_opts(1,1) + grow_len*sin(ang)];
        endpts_x2 = [centroid_opts(1,2); centroid_opts(1,2) + grow_len*cos(ang)];
        endpts_tmp = [round(endpts_y2,0) round(endpts_x2,0)];
    elseif check(1) == 1 && check(2) == 0
        endpts_y2 = [centroid_opts(1,1) - grow_len*sin(ang); centroid_opts(1,1)];
        endpts_x2 = [centroid_opts(1,2) - grow_len*cos(ang); centroid_opts(1,2)];
        endpts_tmp = [round(endpts_y2,0) round(endpts_x2,0)];
    end
    check = ismember(endpts_tmp, blob, 'rows');
    if check(1) == 1 && check(2) == 1
        endpts_pot = endpts_tmp;
    else
        endpts = [];
        return
    end
end
endpts = endpts_pot;

end