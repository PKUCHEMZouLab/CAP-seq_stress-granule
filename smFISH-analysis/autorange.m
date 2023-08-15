% For batch process of contrast improvement in RGB converstion
% The cmin and cmax, who are defined as the intensity values of the pixels
% in the original image that rank at the saturated% and 1-saturated% positions, respectively. 
% The default value 0.0035 is from the "auto" selection in ImageJ
% Peng Zou lab, 2022
function [cmin, cmax, autoimg]=autorange(img,saturated,colormap)
if nargin == 1
    saturated = 0.0035;
    colormap = repmat(((0:255)./255)',1,3);
end
if size(colormap,1) ~=256
    colormap = [(interp1([0 255],[0 colormap(1)],0:255))' (interp1([0 255],[0 colormap(2)],0:255))' (interp1([0 255],[0 colormap(3)],0:255))'];
end
temp = sort(reshape(img,1,numel(img)));
cmin = temp(round(saturated/2*numel(img)));
cmax = temp(round((1-saturated/2)*numel(img)));
temp1 = img;
temp1(temp1<cmin) = cmin;
temp1(temp1>cmax) = cmax;
autoimg = ind2rgb(round((temp1-cmin)./(cmax-cmin)*255), colormap);

