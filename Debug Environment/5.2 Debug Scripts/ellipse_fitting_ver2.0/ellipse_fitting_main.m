%Get the border pixels, pixel-blocks, centroids, and sizes!
function [xnew ynew] = ellipse_fitting_main(img, img_filt, img_45m, img_filt_45m, LS_USA, LS_USA_45m, x_earth, y_earth, imin, jmin, imin_45m, jmin_45m, deg2met, FM, FGIF, d_deg, minlat, minlon, psi0)
% Convert FGIF to pixels and extract endpoints
FGIF_xp = (FGIF(:,1) - minlon)./d_deg;
FGIF_yp = (FGIF(:,2) - minlat)./d_deg;
FGIF_p = [FGIF_xp FGIF_yp];
FGIF_start = [round((x_earth - minlon)/d_deg*2+1,0) round((y_earth - minlat)/d_deg*2+1,0)];
FGIF_p45 = [FGIF_xp.*2 + 1 FGIF_yp.*2 + 1]; % [x, y]
FGIF_p45i = [round(FGIF_p45(:,1),0) round(FGIF_p45(:,2),0)];
high_pop = 10000;
format short g
tic
min_border_pix = 4;
IM = img < 150;
IM_45m = img_45m < 150;
IM_filt = img_filt < 150;
IM_filt_45m = img_filt_45m < 150;

if FM == 0 || FM == 3 || FM == 4 || FM == 5 || FM == 6
    figure(1)
    imagesc(IM_filt)
    hold on 
    title('Filtered Map')
    set(gca,'YDir','normal')
    cmap = [0.5843 0.8157 0.9882; 1 1 1];
    colormap(cmap);
    hold off
    
    if ismember(0,IM) == 0
        disp('All low population zones!')
        land_pts = [];
    elseif ismember(1,IM) == 0
        disp('No low population zones!')
        land_pts = [];
    else
        [blobs, blob_numpix, crash_options, flag] = find_contiguous_areas_bwmorph(img_filt, LS_USA, min_border_pix, x_earth, y_earth, imin, jmin, deg2met);
        if size(crash_options,1) > 0
            crash_pts = sortrows(crash_options, [4 3], {'ascend', 'descend'});
        else
            [blobs, blob_numpix, crash_options, flag] = find_contiguous_areas_bwmorph(img, LS_USA, min_border_pix, x_earth, y_earth, imin, jmin, deg2met);
            if size(crash_options,1) > 0
                crash_pts = sortrows(crash_options, [4 3], {'ascend', 'descend'});
            else
                crash_pts = [];
            end
        end

            %Obtain the ellipses...
            ellipses = []; 
            Skel_mat = 100*ones(size(img_filt));
            row_skel = [];
            col_skel = [];
            skel_brnchpts = [];
            skel_endpts = [];
            for idx=1:size(blobs,1)
                tmpblob = blobs(idx,1:blob_numpix(idx,1),:);
                thisblob = reshape(tmpblob,size(tmpblob,2),size(tmpblob,3));
                [ra, rb, ang, centroid, nodes, endpts, area, high_low_rat, dist, A_thin, xnew_opt, ynew_opt] = get_pixelcore_fit_ellipses ( thisblob, 0.5, img_filt, x_earth, y_earth, imin, jmin, LS_USA, deg2met);

                %Format skeleton matrix
                [row_skel1 col_skel1] = find(A_thin == -9999);
                row_skel = [row_skel; row_skel1];
                col_skel = [col_skel; col_skel1];
                skel_brnchpts = [skel_brnchpts; nodes];
                skel_endpts = [skel_endpts; endpts];
                this_ellipse_set = [ra rb ang centroid area high_low_rat dist xnew_opt ynew_opt];
                ellipses = [ellipses; this_ellipse_set];
            end

            if size(ellipses, 1) > 0
                land_pts = [ellipses(:,9) ellipses(:,10) ellipses(:,3).*180/pi ellipses(:,6) ellipses(:,7) ellipses(:,8)]; %[x, y, orientation, area, high:low pop, dist]
                land_pts = sortrows(land_pts, [5 4], {'ascend' 'descend'});
                land_pts_wcent = [ellipses(:,9) ellipses(:,10) ellipses(:,3).*180/pi ellipses(:,6) ellipses(:,7) ellipses(:,8) ellipses(:,4) ellipses(:,5)]; 
                land_pts_wcent = sortrows(land_pts_wcent, [5 4], {'ascend' 'descend'});

                xnew_pix = land_pts_wcent(1,8);
                ynew_pix = land_pts_wcent(1,7);

            else 
                land_pts = [];
                xnew_pix = [];
                ynew_pix = [];
            end
            time = toc

            %plot skeleton image 
            if size(row_skel,1) > 0
                for i = 1:length(row_skel)
                    Skel_mat(row_skel(i), col_skel(i)) = -9999;
                end
                figure(2)
                display_pixel_array(Skel_mat);
                hold on 
                title('Skeletonized Image')
                if size(skel_brnchpts,1) > 0 
                plot(skel_brnchpts(:,2), skel_brnchpts(:,1), '.c', 'Markersize', 20)
                end
                if size(skel_endpts,1) > 0
                plot(skel_endpts(:,2), skel_endpts(:,1), '.r', 'Markersize', 20)
                end
                hold off
            end

            if size(ellipses,1) > 0
                figure(3)
                imagesc(IM)
                hold on 
                title('Fitted Ellipses')
                set(gca,'YDir','normal')
                cmap = [0.5843 0.8157 0.9882; 1 1 1];
                colormap(cmap);
               % ...and display them
                x=[]; y=[];
                for idx=1:size(ellipses,1)
                    if ellipses(idx,1)>0
                        [x_cur,y_cur]=ellipse(ellipses(idx,1),ellipses(idx,2),ellipses(idx,3),ellipses(idx,5),ellipses(idx,4),'c');
                        tmp_x = x; tmp_y = y;
                        x = [tmp_x; x_cur']; y = [tmp_y; y_cur'];
                    end
                end
                plot(x,y,'r.','LineWidth',1);
                plot(xnew_pix, ynew_pix, '.k', 'Markersize', 15)
                legend('Fitted Ellipses', 'Chosen Landing wpt')
                hold off
            end
        end
    
elseif FM == 1 || FM == 3 || FM == 8
    figure(1)
    imagesc(IM_filt_45m)
    hold on 
    title('High Res Filtered')
    set(gca,'YDir','normal')
    plot(FGIF_p45(:,1), FGIF_p45(:,2))
    cmap = [0.5843 0.8157 0.9882; 1 1 1];
    colormap(cmap);
    hold off
    
    figure(2)
    imagesc(IM_45m)
    hold on 
    title('High Res Map')
    set(gca,'YDir','normal')
    cmap = [0.5843 0.8157 0.9882; 1 1 1];
    colormap(cmap);
    hold off
    
    if ismember(0,IM) == 0
       disp('All low population zones!')
       land_pts = [];
       [blobs, blob_numpix, crash_options, flag] = find_contiguous_areas_noturn(img_45m, LS_USA_45m, min_border_pix, x_earth, y_earth, imin_45m, jmin_45m, deg2met, FGIF_p45, high_pop);
        if size(crash_options,1) > 0
            crash_pts = sortrows(crash_options, [4 3], {'ascend', 'descend'});
        else
            crash_pts = [];
        end
    elseif ismember(1,IM) == 0
       disp('No low population zones!')
       land_pts = [];    
       crash_pts = [];
    else
        [blobs, blob_numpix, crash_options, flag] = find_contiguous_areas_noturn(img_45m, LS_USA_45m, min_border_pix, x_earth, y_earth, imin_45m, jmin_45m, deg2met, FGIF_p45, high_pop);
        if size(crash_options,1) > 0
            crash_pts = sortrows(crash_options, [4 3], {'ascend', 'descend'});
        else
            crash_pts = [];
        end
    end  
        %Obtain the ellipses...
        ellipses = []; 
        for idx=1:size(blobs,1)
            tmpblob = blobs(idx,1:blob_numpix(idx,1),:);
            thisblob = reshape(tmpblob,size(tmpblob,2),size(tmpblob,3));
            [ra, rb, ang, centroid, nodes, endpts, area, high_low_rat, dist, A_thin, xnew_opt, ynew_opt] = get_pixelcore_fit_ellipses_noturn ( thisblob, 0.5, img_45m, x_earth, y_earth, imin_45m, jmin_45m, LS_USA_45m, deg2met, FGIF_p45, d_deg);
            this_ellipse_set = [ra rb ang centroid area high_low_rat dist xnew_opt ynew_opt];
            ellipses = [ellipses; this_ellipse_set];
        end
        if size(ellipses, 1) > 0 && nnz(ellipses(:,9)) ~= 0 && nnz(ellipses(:,10)) ~= 0
            land_pts = [ellipses(:,9) ellipses(:,10) ellipses(:,3).*180/pi ellipses(:,6) ellipses(:,7) ellipses(:,8)]; %[x, y, orientation, area, high:low pop, dist]
            land_pts = sortrows(land_pts, [5 4], {'ascend' 'descend'});
        else 
            land_pts = [];
        end
        time = toc

        if size(ellipses,1) > 0
            figure(2)
            imagesc(IM_45m)
            hold on 
            title('Fitted High Resolution Ellipse')
            set(gca,'YDir','normal')
            cmap = [0.5843 0.8157 0.9882; 1 1 1];
            colormap(cmap);
           % ...and display them
            x=[]; y=[];
            for idx=1:size(ellipses,1)
                if ellipses(idx,1)>0
                    [x_cur,y_cur]=ellipse(ellipses(idx,1),ellipses(idx,2),ellipses(idx,3),ellipses(idx,5),ellipses(idx,4),'c');
                    tmp_x = x; tmp_y = y;
                    x = [tmp_x; x_cur']; y = [tmp_y; y_cur'];
                end
            end
            plot(x,y,'r.','LineWidth',1);
            hold off

            figure(3)
            imagesc(IM)
            hold on 
            title('Fitted Ellipse')
            set(gca,'YDir','normal')
            cmap = [0.5843 0.8157 0.9882; 1 1 1];
            colormap(cmap);
           % ...and display them
            x=[]; y=[];
            for idx=1:size(ellipses,1)
                if ellipses(idx,1)>0
                    [x_cur,y_cur]=ellipse(ellipses(idx,1)/2,ellipses(idx,2)/2,ellipses(idx,3),ellipses(idx,5)/2,ellipses(idx,4)/2,'c');
                    tmp_x = x; tmp_y = y;
                    x = [tmp_x; x_cur']; y = [tmp_y; y_cur'];
                end
            end
            plot(x+1,y,'r.','LineWidth',1);

            if isempty(land_pts) == 0
                plot(ellipses(1,5)/2+1, ellipses(1,4)/2, '.k', 'Markersize', 15)
            end
            legend('Fitted Ellipse', 'Chosen Landing wpt')
            hold off
        end
            
elseif FM == 2 || FM == 7 
    figure(1)
    imagesc(IM_filt_45m)
    hold on 
    title('High Res Filtered')
    set(gca,'YDir','normal')
    cmap = [0.5843 0.8157 0.9882; 1 1 1];
    plot(FGIF_p45(:,1), FGIF_p45(:,2), '.k', 'Markersize', 10)
    legend('FGIF')
    colormap(cmap);
    hold off
    
    figure(2)
    imagesc(IM_45m)
    hold on 
    title('High Res Map')
    set(gca,'YDir','normal')
    plot(FGIF_p45(:,1), FGIF_p45(:,2), '.k', 'Markersize', 10)
    cmap = [0.5843 0.8157 0.9882; 1 1 1];
    legend('FGIF')
    colormap(cmap);
    hold off
    
    flg_img = 0;
    flg_filt = 0;
    if ismember(0,IM) == 0
        disp('All low population zones!')
        land_pts = [];
    elseif ismember(1,IM) == 0
        disp('No low population zones!')
        land_pts = [];
    else
        [blobs, blob_numpix, crash_options, flag] = find_contiguous_areas_noElev(img_filt_45m, LS_USA_45m, min_border_pix, x_earth, y_earth, imin_45m, jmin_45m, deg2met, FGIF_p45, high_pop);
        if size(crash_options,1) > 0
            crash_pts = sortrows(crash_options, [4 3], {'ascend', 'descend'});
            flg_filt = 1;
        else
            [blobs, blob_numpix, crash_options, flag] = find_contiguous_areas_noElev(img_45m, LS_USA_45m, min_border_pix, x_earth, y_earth, imin_45m, jmin_45m, deg2met, FGIF_p45, high_pop);
            flg_img = 1;
            if size(crash_options,1) > 0
                crash_pts = sortrows(crash_options, [4 3], {'ascend', 'descend'});
            else
                crash_pts = [];
            end
        end
            %Obtain the ellipses...
            ellipses = []; 
            for idx=1:size(blobs,1)
                tmpblob = blobs(idx,1:blob_numpix(idx,1),:);
                thisblob = reshape(tmpblob,size(tmpblob,2),size(tmpblob,3));
                if flg_img == 1
                    [ra, rb, ang, centroid, nodes, endpts, area, high_low_rat, dist, A_thin, xnew_opt, ynew_opt] = get_pixelcore_fit_ellipses_noElev ( thisblob, 0.5, img_45m, x_earth, y_earth, imin_45m, jmin_45m, LS_USA_45m, deg2met, FGIF_p45, FGIF_start, d_deg, psi0);
                elseif flg_filt == 1
                    [ra, rb, ang, centroid, nodes, endpts, area, high_low_rat, dist, A_thin, xnew_opt, ynew_opt] = get_pixelcore_fit_ellipses_noElev ( thisblob, 0.5, img_filt_45m, x_earth, y_earth, imin_45m, jmin_45m, LS_USA_45m, deg2met, FGIF_p45, FGIF_start, d_deg, psi0);
                end
                this_ellipse_set = [ra rb ang centroid area high_low_rat dist xnew_opt ynew_opt];
                ellipses = [ellipses; this_ellipse_set];
            end
            if size(ellipses, 1) > 0
                land_pts = [ellipses(:,9) ellipses(:,10) ellipses(:,3).*180/pi ellipses(:,6) ellipses(:,7) ellipses(:,8)]; %[x, y, orientation, area, high:low pop, dist]
                land_pts = sortrows(land_pts, [5 6 4], {'ascend' 'descend' 'descend'});
            else 
                land_pts = [];
            end
            time = toc

            if size(ellipses,1) > 0
                figure(3)
                imagesc(IM_45m)
                hold on 
                title('Fitted High Resolution Ellipses')
                plot(FGIF_p45(:,1), FGIF_p45(:,2), '.k', 'Markersize', 10)
                set(gca,'YDir','normal')
                cmap = [0.5843 0.8157 0.9882; 1 1 1];
                colormap(cmap);
               % ...and display them
                x=[]; y=[];

                for idx=1:size(ellipses,1)
                    if ellipses(idx,1)>0
                        [x_cur,y_cur]=ellipse(ellipses(idx,1),ellipses(idx,2),ellipses(idx,3),ellipses(idx,5),ellipses(idx,4),'c');
                        tmp_x = x; tmp_y = y;
                        x = [tmp_x; x_cur']; y = [tmp_y; y_cur'];
                    end
                end
                plot(x,y,'r.','LineWidth',1);
                plot(ellipses(1,5), ellipses(1,4), '.g', 'Markersize', 15)
                legend('FGIF', 'Fitted Ellipses', 'Chosen Landing wpt')
                hold off

                figure(6)
                imagesc(IM)
                hold on 
                title('Fitted Ellipses')
                plot(FGIF_p45(:,1)./2, FGIF_p45(:,2)./2, '.k', 'Markersize', 10)
                set(gca,'YDir','normal')
                cmap = [0.5843 0.8157 0.9882; 1 1 1];
                colormap(cmap);
               % ...and display them
                x=[]; y=[];
                for idx=1:size(ellipses,1)
                    if ellipses(idx,1)>0
                        [x_cur,y_cur]=ellipse(ellipses(idx,1)/2,ellipses(idx,2)/2,ellipses(idx,3),ellipses(idx,5)/2,ellipses(idx,4)/2,'c');
                        tmp_x = x; tmp_y = y;
                        x = [tmp_x; x_cur']; y = [tmp_y; y_cur'];
                    end
                end
                plot(x,y,'r.','LineWidth',1);    
                plot(ellipses(1,5)/2, ellipses(1,4)/2, '.g', 'Markersize', 15)
                legend('FGIF', 'Fitted Ellipses', 'Chosen Landing wpt')
                hold off
            end
        end
end
crash_pts
land_pts
if size(land_pts,1) > 0
    xnew = land_pts(1,1);
    ynew = land_pts(1,2);
elseif size(crash_pts,1) > 0
    xnew = crash_pts(1,1);
    ynew = crash_pts(1,2);
else
    xnew = x_earth;
    ynew = y_earth;
end
