function A3 = thin_skeleton(img)
IM = img < 100;
A3 =  bwmorph(IM,'skel',inf);
A3 = double(A3);
A3(A3 == 1) = -9999;
A3(A3 == 0) = 100;
end
