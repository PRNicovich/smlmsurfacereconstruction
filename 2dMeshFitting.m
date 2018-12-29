%% Generate test pattern ground truth
% Drawing a circle

circRadius = 100; % microns
circCenter = [-25, -40];

th = 0:pi/100:2*pi;
xunit = circRadius * cos(th) + circCenter(1);
yunit = circRadius * sin(th) + circCenter(2);
h = plot(xunit, yunit);

figure(1)
plot(xunit, yunit); 



% Simulate points around ground truth

sizeX = 25.6; % microns
sizeY = 25.6; % microns

pixelSize = 0.1; % microns

nPoints = 100;
nDetectsPerPoint = 20;
locPrecision = 20; % nm



% Render points into image

% Find ridge line

% Sub-sample at some reasonable spacing

% Resample to spacing appropriate for density

% Fit resampled points to data cloud

% Display