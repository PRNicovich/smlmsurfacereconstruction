% Run with current working directory as .\smlmmeshfitting\

% Requires:
% stlread
% stlwrite
% ICP_finite

% Hausdorff distance
% Zachary Danziger (2020). 
% Hausdorff Distance (https://www.mathworks.com/matlabcentral/fileexchange/26738-hausdorff-distance), 
% MATLAB Central File Exchange. Retrieved July 8, 2020.
% MATLAB 2019a

%% Input parameters

meshObj = meshfitter();

%% Example 2 - Stanford bunny

% Transform Stanford bunny model to scale, position desired

meshObj = meshFitting3DFcn(meshObj);

%%

meshObj.analysis = analyzeMeshFittingResults(meshObj.inputs, meshObj.results, meshObj.analysis);

%%
save('bunnySurfaceOutput.mat', 'meshObj');


% %% Plot things on calculated surface
% figure(4)
% plot3(zComb(:,1), zComb(:,2), zComb(:,3), 'k.', 'markersize', 1)


