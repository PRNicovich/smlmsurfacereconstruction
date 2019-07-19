
%% Simulate points around ground truth
% This section can be substituted with other means to generate SMLM
% dataset. Output should be N x 3 matrix of SMLM points named 'pts'.

% Double-check units here
sizeX = 10000; % nm
sizeY = 10000; % nm
sizeZ = 10000; % nm

nPoints = 40000;
nDetectsPerPoint = 12; % mean, poisson distributed

exampleDataFile = 'C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\testData\CF647 KRas_8.txt';

template.Format = 'stlFile';
% template.File = 'C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\High_Resolution_Stanford_Bunny\FLATFOOT_StanfordBunny_jmil_HIGH_RES_Smoothed.stl';
template.File = 'C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\parametricSurface.stl';
% template.Scale = 0.08;
template.Scale = 0.0005;
template.Center = [5000, 5000, 3000];

xSlice = 60; % Thickness of single slice, in nm
subsampleStep = 10; % Subsample of fitted curve, in points
onPoint = 5; % Start curve on this point.  Avoids deviation from center of point clound from using min span tree to find ridge


rotationAroundX = pi/10; % radians
rotationAroundZ = pi/4; % radians

zSlice = xSlice/2;



%%
[z, zComb, polygonReturn] = meshFitting3DFcn(sizeX, sizeY, sizeZ, nPoints, exampleDataFile, nDetectsPerPoint, ...
                                                        template, xSlice, subsampleStep, onPoint, ...
                                                        rotationAroundX, rotationAroundZ, zSlice);
                                                    
save('parametricSurfaceOutput.mat');


template.Format = 'stlFile';
template.File = 'C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\High_Resolution_Stanford_Bunny\FLATFOOT_StanfordBunny_jmil_HIGH_RES_Smoothed.stl';
template.Scale = 0.08;
template.Center = [5000, 5000, 5000];
[z, zComb, polygonReturn] = meshFitting3DFcn(sizeX, sizeY, sizeZ, nPoints, exampleDataFile, nDetectsPerPoint, ...
                                                        template, xSlice, subsampleStep, onPoint, ...
                                                        rotationAroundX, rotationAroundZ, zSlice);
                                                    
save('bunnySurfaceOutput.mat');


% %% Plot things on calculated surface
% figure(4)
% plot3(zComb(:,1), zComb(:,2), zComb(:,3), 'k.', 'markersize', 1)


