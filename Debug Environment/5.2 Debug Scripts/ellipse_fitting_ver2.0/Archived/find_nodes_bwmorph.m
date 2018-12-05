function [nodes, endpts] = find_nodes_bwmorph(img_skel)
branchpts = bwmorph(img_skel, 'branchpoints');
end_idx = bwmorph(img_skel, 'endpoints');
[row_branch col_branch] = find(double(branchpts) == 1);
[row_endpts col_endpts] = find(double(end_idx) == 1);

nodes = [col_branch(:) row_branch(:)];
endpts = [col_endpts(:) row_endpts(:)];
end
