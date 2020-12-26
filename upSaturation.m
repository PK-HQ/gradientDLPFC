function [colorOut]=upSaturation(colorRGB,nColors)
color(1:9,1,1:3)=colorRGB;
colorHSV = rgb2hsv(color);
% "20% more" saturation:
colorHSV(:, :, 2) = colorHSV(:, :, 2) * 1.5;
% or add:
% HSV(:, :, 2) = HSV(:, :, 2) + 0.2;
colorHSV(colorHSV > 1) = 1;  % Limit values
colorRGB = hsv2rgb(colorHSV);
colorOut(1:9,1:3)=colorRGB;
colorOut=colorOut(end-nColors+1:end,:);
end