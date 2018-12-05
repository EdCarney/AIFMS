function [ stat ] = diplay_pixel_array( img )
%display_pixel_array: given an array of pixels, display in consistent
%format, with 1:1 aspect ratio and gca
    image(img);
    set(gca,'YDir','normal');
    daspect([1 1 1]);
    stat = 1;
end

