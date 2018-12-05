function [A4] = thin_skeleton(img)
IM = img < 100;
A4 =  bwmorph(IM,'skel',inf);
A3 = bwmorph(A4, 'thin', inf);
A4 = double(A3);
A4(A4 == 1) = -9999;
A4(A4 == 0) = 100;
end
