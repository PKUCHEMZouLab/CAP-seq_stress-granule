% For convertion of input wavelength (nm) value to a RGB array
% Peng Zou lab, 2018
function [RGB]=getrgb(WL);
if (WL>=380)&&(WL<440)
    R = -1*(WL-440)/(440-380)*255/256;
    G = 0;
    B = 255/256;
elseif (WL>=440)&&(WL<490)
    R = 0;
    G = (WL-440)/(490-440);
    B = 255/256;
elseif (WL>=490)&&(WL<510)
    R = 0;
    G = 255/256;
    B = -1*(WL-510)/(510-490)*255/256;
elseif(WL>=510)&&(WL<580)
    R = (WL-510)/(580-510)*255/256;
    G = 255/256;
    B = 0;
elseif(WL>=580)&&(WL<645)
    R = 255/256;
    G = -1*(WL-645)/(645-580)*255/256;
    B = 0;
else(WL>=645)&&(WL<=780)
    R = 255/256;
    G = 0;
    B = 0;
end
RGB = [R G B];
