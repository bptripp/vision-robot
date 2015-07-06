% This script calculates lens distortion needed to emulate macaque retinal 
% ganglion cell density, taken from Wassle et al. (1989) with 
% http://arohatgi.info/WebPlotDigitizer/app/. 
% 
% It produces a lens distortion plot and an example image with the
% desired distortion. 

% run from the directory with this file 
wassle = csvread('wassle data.csv');

[~,I] = sort(wassle(:,1));
wassle = wassle(I,:);
offset = wassle(:,1);
linearDensity = wassle(:,2).^.5;

% We will go out to 16mm in the nasal (higher density) direction
ind = find(abs(offset) <= 16 & offset <= .15);
offset = abs(offset(ind));
offset(end) = 0; % I think non-zero was an error in figure extraction
linearDensity = linearDensity(ind);
offset = offset(end:-1:1); % distance from fovea centre
linearDensity = linearDensity(end:-1:1);

% plot Wassle data
figure
scatter(offset, linearDensity)

% To be conservative about motion sensitivity and lower optical performance
% in periphery, we can increase density in the periphery. 
if 0
    scaleFactor = 1 + 0.5*offset/max(offset);
    linearDensity = linearDensity .* scaleFactor;
end

% plot(offset, linearDensity, '.'), set(gca, 'YLim', [0 1000])

% Assume 16mm from fovea corresponds to 90 degrees (~full range in nasal 
% direction from Wassle et al. (1989) and ~full macaque field of view from 
% Van Essen et al. (1984) Vision Research. 
FOV = 45; % note: this is half the full FOV
max_e = 16; % mm
max_a = 90; % degrees

% trapezoidal integration of linear density to find cumulative cells 

ecc = FOV/max_a *max_e;
[~, idx] = min( abs(offset-ecc) );
    
cells = (linearDensity(1:idx-1) + linearDensity(2:idx)) / 2 .* diff(offset(1:idx));
cumulativeCells = cumsum(cells);

relativeOffset = offset(1:idx) / offset(idx);
relativeCumulativeCells = [0; cumulativeCells] / cumulativeCells(idx-1);

angle = offset(1:idx) * FOV*pi/180 / ecc;

% cells = (linearDensity(1:end-1) + linearDensity(2:end)) / 2 .* diff(offset);
% cumulativeCells = cumsum(cells);
% 
% relativeOffset = offset / offset(end);
% relativeCumulativeCells = [0; cumulativeCells] / cumulativeCells(end);

%%%% figure for lens designer %%%%
figure
subplot(1,2,1), hold on
plot(relativeOffset, relativeCumulativeCells, 'k')
plot([0 1], [0 1], 'k--')
legend('Desired','Linear Reference')
xlabel('Fraction of FOV from Centre')
ylabel('Fraction of Sensor from Centre')
title_str = sprintf('%d degree FOV', FOV*2);
title(title_str);
% set(gca, 'FontSize', 18)
subplot(1,2,2)
plot(relativeOffset, relativeCumulativeCells./relativeOffset*100-100, 'k')
xlabel('Fraction of FOV from Centre')
ylabel('% Distortion')
title_str = sprintf('%d degree FOV', FOV*2);
title(title_str);
% set(gca, 'FontSize', 18)
set(gcf, 'Position', [440 440 881 358])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot a similarly warped example figure ... 
sCam1 = 1048; % Flea3 1.3MP resolution is 1328 x 1048
sCam2 = 1328;
cameraPixel = (1:sCam2/2)-1/2; 
maxCameraAngle = max(angle);
cameraAngle = interp1(relativeCumulativeCells, angle, cameraPixel/cameraPixel(end));

% A = imread('cars.jpg', 'jpg');
A = imread('mountains.jpg', 'jpg');
sIm1 = size(A,1);
sIm2 = size(A,2);
imagePixel = (1:sIm2/2)-1/2;
maxImageAngle = 90*pi/180;
imageAngle = imagePixel * maxImageAngle / imagePixel(end);
mappedPixel = interp1(imageAngle, imagePixel, cameraAngle);

cameraPixelsHor = repmat((1:sCam2)-(sCam2/2)-1/2, sCam1, 1);
cameraPixelsVer = repmat((1:sCam1)'-(sCam1/2)-1/2, 1, sCam2);
cameraRadialPixel = (cameraPixelsHor.^2 + cameraPixelsVer.^2).^.5;
cameraRadialAngle = interp1(cameraPixel, cameraAngle, cameraRadialPixel);
mappedRadialPixel = interp1(imageAngle, imagePixel, cameraRadialAngle);
angles = atan(cameraPixelsVer ./ cameraPixelsHor);
imagePixelsHor = mappedRadialPixel .* cos(angles);
imagePixelsHor(:,1:sCam2/2) = -imagePixelsHor(:,1:sCam2/2);
imagePixelsVer = mappedRadialPixel .* sin(angles);
imagePixelsVer(:,1:sCam2/2) = -imagePixelsVer(:,1:sCam2/2);

lookupVer = round(imagePixelsVer+size(A,1)/2+1/2);
lookupVer(lookupVer < 1) = NaN;
lookupVer(lookupVer > size(A,1)) = NaN;
lookupVer(isnan(lookupVer)) = 1;
lookupHor = round(imagePixelsHor+size(A,2)/2+1/2);
lookupHor(lookupHor < 1) = NaN;
lookupHor(lookupHor > size(A,2)) = NaN;
lookupHor(isnan(lookupHor)) = 1;

remapIndices = (lookupHor(:)-1)*size(A,1) + lookupVer(:);
flatA = reshape(A, size(A,1)*size(A,2), 3);
remapped = flatA(remapIndices, :);
remapped = reshape(remapped, sCam1, sCam2, 3);

figure, imshow(A)
title('Original image')
figure, imshow(uint8(remapped))
title_str = sprintf('Lens Distorted Image, %d degree FOV', FOV*2);
title(title_str) 

