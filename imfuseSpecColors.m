function imageOut = imfuseSpecColors(chan1, chan2, color1, color2, whiteOrBlack)
%imfuseSpecColors - generate fused image of chan1 + chan2 given colors
%   color1 and color2
%   Output is uint8 RGB image w/ max scaled to 255 over white or black bkgd

c1 = (chan1 - min(chan1(:)))/(max(chan1(:)) - min(chan1(:)));
c2 = (chan2 - min(chan2(:)))/(max(chan2(:)) - min(chan2(:)));

if strcmp(whiteOrBlack, 'white');
    color1 = 1-color1;
    color2 = 1-color2;

    c1c = 1 - cat(3, c1*color1(1), c1*color1(2), c1*color1(3));
    c2c = 1 - cat(3, c2*color2(1), c2*color2(2), c2*color2(3));

    imageOut = (c1c + c2c - 1);
    
elseif strcmp(whiteOrBlack, 'black')

    c1c = cat(3, c1*color1(1), c1*color1(2), c1*color1(3));
    c2c = cat(3, c2*color2(1), c2*color2(2), c2*color2(3));

    imageOut = c1c + c2c;

end


end

