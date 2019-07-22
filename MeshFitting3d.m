
%% Input parameters

% Size of volume bounding box
sizeX = 10000; % nm
sizeY = 10000; % nm
sizeZ = 10000; % nm

nPoints = 40000; % Number of emitters to generate
nDetectsPerPoint = 12; % mean number of detection events per emitter, poisson distributed

% File containing real SMLM data.  Used to generate localization precision
% distribution
exampleDataFile = 'C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\testData\CF647 KRas_8.txt';

% 3D template parameters
template.Format = 'stlFile'; % Type.  'stlFile' supported. 
% template.File = 'C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\High_Resolution_Stanford_Bunny\FLATFOOT_StanfordBunny_jmil_HIGH_RES_Smoothed.stl';
template.File = 'C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\parametricSurface.stl'; % Reference STL file
% template.Scale = 0.08;
template.Scale = 0.0005; % Scaling factor for native STL file size -> size for this simulation. 
template.Center = [5000, 5000, 3000]; % Location of center of mass of STL file in simulated space

% Meridian cutting and 2D point cloud approximation
xSlice = 60; % Thickness of single slice, in nm
subsampleStep = 10; % Subsample of fitted curve, in points
onPoint = 5; % Start curve on this point.  Avoids deviation from center of point clound from using min span tree to find ridge

% Rotation of point cloud in analysis 
rotationAroundX = pi/10; % Rotation around X axis to avoid bad analysis of area parallel to XY plane; in radians
rotationAroundZ = pi/4; % Rotation around Z to re-slice second pass of ponit cloud fitting; in radians

% Step size between resampled Z slices. 
zSlice = xSlice/2;



%% Example 1 - Open surface of sine wave
[z, zComb, polygonReturn] = meshFitting3DFcn(sizeX, sizeY, sizeZ, nPoints, exampleDataFile, nDetectsPerPoint, ...
                                                        template, xSlice, subsampleStep, onPoint, ...
                                                        rotationAroundX, rotationAroundZ, zSlice);
                                                    
save('parametricSurfaceOutput.mat');

%% Example 2 - Stanford bunny
% template.Format = 'stlFile';
% template.File = 'C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\High_Resolution_Stanford_Bunny\FLATFOOT_StanfordBunny_jmil_HIGH_RES_Smoothed.stl';
% template.Scale = 0.08;
% template.Center = [5000, 5000, 5000];
% [z, zComb, polygonReturn] = meshFitting3DFcn(sizeX, sizeY, sizeZ, nPoints, exampleDataFile, nDetectsPerPoint, ...
%                                                         template, xSlice, subsampleStep, onPoint, ...
%                                                         rotationAroundX, rotationAroundZ, zSlice);
%                                                     
% save('bunnySurfaceOutput.mat');


% %% Plot things on calculated surface
% figure(4)
% plot3(zComb(:,1), zComb(:,2), zComb(:,3), 'k.', 'markersize', 1)


