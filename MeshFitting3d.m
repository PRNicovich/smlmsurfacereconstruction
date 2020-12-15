% Demo script for SMLM surface reconstruction
% Run with current working directory as .\smlmmeshfitting\

%% Input parameters

meshObj = meshfitter();

%% Example - Stanford bunny recon from STL file

meshObj = meshFitting3DFcn(meshObj);

%% Analysis on resulting mesh

meshObj.analysis = analyzeMeshFittingResults(meshObj.inputs, meshObj.results, meshObj.analysis);

%% Save output
save('bunnySurfaceOutput.mat', 'meshObj');

%% Plots
% Input in black. Point cloud in blue. Output in red.


% Plot input model
figure(1);
title("Input model");
tri = triangulation(meshObj.inputs.groundTruth.faces, ...
    meshObj.inputs.groundTruth.vertices(:,1), ...
    meshObj.inputs.groundTruth.vertices(:,2), ...
    meshObj.inputs.groundTruth.vertices(:,3));
triPlot = trisurf(tri, 'edgecolor', 'none');
set(triPlot, 'faceColor', [0.6 0.6 0.6]);
material('dull');
camlight('headlight');
set(triPlot, ...
    'FaceLighting',    'gouraud', ...
    'AmbientStrength', 0.001)
lightHand = findobj('parent', gca, 'type', 'light');
axis('image');
title("Input model");
xlabel('X position (nm)');
ylabel('Y position (nm)');
zlabel('Z position (nm)');

figure(2);
ptsCld = plot3(meshObj.results.z(:,1), ...
               meshObj.results.z(:,2), ...
               meshObj.results.z(:,3), '.', ...
               'markersize', 1, 'markeredgecolor', [0.2, 0.4, 1]);
axis('image');
grid('on');
title("Point cloud");
xlabel('X position (nm)');
ylabel('Y position (nm)');
zlabel('Z position (nm)');


figure(3);
tri = triangulation(meshObj.results.polygonReturn.faces, ...
    meshObj.results.polygonReturn.vertices(:,1), ...
    meshObj.results.polygonReturn.vertices(:,2), ...
    meshObj.results.polygonReturn.vertices(:,3));
triPlot = trisurf(tri, 'edgecolor', 'none');
set(triPlot, 'faceColor', [1 0.6 0.4]);
material('dull');
camlight('headlight');
set(triPlot, ...
    'FaceLighting',    'gouraud', ...
    'AmbientStrength', 0.001)
lightHand = findobj('parent', gca, 'type', 'light');
axis('image');
title("Surface reconstruction");
xlabel('X position (nm)');
ylabel('Y position (nm)');
zlabel('Z position (nm)');


