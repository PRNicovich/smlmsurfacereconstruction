
function [z, zComb, polygonReturn, pwOut, groundTruth] = meshFitting3DFcn(sizeX, sizeY, sizeZ, nPoints, exampleDataFile, nDetectsPerPoint, ...
                                                        template, xSlice, subsampleStep, onPoint, ...
                                                        rotationAroundX, rotationAroundZ, zSlice)
%----------------------------------------------------------%

    ImportData = Import1File(exampleDataFile);
    data = ImportData.Data;

    % Make a sphere
    % circRadius = 700; % nm
    circCenter = template.Center;
    % [xunit, yunit, zunit] = sampleSpherePoints(nPoints, circRadius, circCenter);

    % Make a teapot
    % [xunit, yunit, zunit] = GenerateTestDataFromSTLFile('C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\classic-teapot\TEAPOT.stl', ...
    %     nPoints, circCenter, .07);

    if strcmp(template.Format, 'stlFile')
        [xunit, yunit, zunit] = GenerateTestDataFromSTLFile(template.File, ...
            nPoints, circCenter, template.Scale);
    end
    
    groundTruth = {xunit, yunit, zunit};

    % From zero-noise points sampled on STL file, make some proxy SMLM data
    z = pointsToSMLMPointCloud(xunit, yunit, zunit, nPoints, nDetectsPerPoint,...
                                data(:,7:8));

    % Filter points to remain in ROI bounds
    z = enforceROIBounds(z, sizeX, sizeY, sizeZ);

    %% Rotate input point cloud to avoid having one big plane parallel to the z axis
    % First step in 3d reconstruction analysis.  Required if object has a large
    % area that is parallel to one of the spatial planes.  Adding a tilt to the
    % object works around this problem. 

    zRot = rotateCloudAroundAxis(z, rotationAroundX, 'x');

    %% Loop over slices in x
    % Main fitting function along single axis

    ptRange = [max(zRot(:,(1:3))); min(zRot(:,1:3))];
    
    xPiecewise = fitPointSetLoop(zRot, ptRange(2,1), ptRange(1,1), xSlice, 2, subsampleStep, onPoint);
    
    
    %% Loop over slices in y
    % Repeat above section with same parameters but on perpendicular axis.
    yPiecewise = fitPointSetLoop(zRot, ptRange(2,2), ptRange(1,2), xSlice, 1, subsampleStep, onPoint);
    
    %% Loop over slices in z
    % Repeat above section with same parameters but on perpendicular axis.
    zPiecewise = fitPointSetLoop(zRot, ptRange(2,2), ptRange(1,2), xSlice, 3, subsampleStep, onPoint);


    %% Loop over slices between X and Y axes
    % Rotate point cloud by 45 degrees and go again.  This is done to add
    % points to next step and generate smoother cloud. 
    % Rotation amount should be sufficient to avoid having much of object faces
    % to be close to parallel to either X or Y axis. 

    rotPts = rotateCloudAroundAxis(zRot, rotationAroundZ, 'z');

    ptRange = [max(rotPts(:,(1:3))); min(rotPts(:,1:3))];
    bPiecewise = fitPointSetLoop(rotPts, ptRange(2,1), ptRange(1,1), xSlice, 2, subsampleStep, onPoint);
    cPiecewise = fitPointSetLoop(rotPts, ptRange(2,2), ptRange(1,2), xSlice, 1, subsampleStep, onPoint);
    dPiecewise = fitPointSetLoop(rotPts, ptRange(2,2), ptRange(1,2), xSlice, 3, subsampleStep, onPoint);

    %% Double-check for erroneously long links and split if needed
    % Determine threshold for splitting
    % Can feed eithe xPiecewise or yPiecewise values
    % Do both and then take mean as best choice
    [threshX, splitDistX] = rosinThreshold(xPiecewise, 'doPlot', true);
    [threshY, splitDistY] = rosinThreshold(yPiecewise, 'doPlot', true);
    [threshZ, splitDistZ] = rosinThreshold(zPiecewise, 'doPlot', true);

    [threshB, splitDistB] = rosinThreshold(bPiecewise, 'doPlot', true);
    [threshC, splitDistC] = rosinThreshold(cPiecewise, 'doPlot', true);
    [threshD, splitDistD] = rosinThreshold(dPiecewise, 'doPlot', true);

    threshVal = (threshX + threshY + threshZ + threshB + threshC + threshD)/6;

    %%

    % Split up cells at bad links
    % Generates extra points that may or may not matter. 
    % Def an error, but don't care enough to find out right now.

    xPiecewise = splitByLargeDistance(xPiecewise, splitDistX, threshVal);
    yPiecewise = splitByLargeDistance(yPiecewise, splitDistY, threshVal);
    zPiecewise = splitByLargeDistance(zPiecewise, splitDistZ, threshVal);
    
    bPiecewise = splitByLargeDistance(bPiecewise, splitDistB, threshVal);
    cPiecewise = splitByLargeDistance(cPiecewise, splitDistC, threshVal);
    dPiecewise = splitByLargeDistance(dPiecewise, splitDistD, threshVal);
    
    pwOut = {xPiecewise, yPiecewise, zPiecewise, bPiecewise, cPiecewise, dPiecewise};

%     %% Z slice fitted curves
%     % Generate regularly-spaced points that hopefully don't have too much trash
%     % in the set
% 
%     % Make all connected point pairs into a single variable
%     % No way this isn't done better in a graph class, but don't know how to
%     % implement that....

%     zVect = ptRange(2,3) : zSlice : ptRange(1,3);
% 
%     zPiecewise = cell(numel(zVect), 2);
%     zP = [];
%     zP2 = [];
% 
%     for k = 1:numel(zVect)
% 
%         zP = [zP; generateZSlicePoints(xPiecewise, zVect(k), 'x')];
%         zP = [zP; generateZSlicePoints(yPiecewise, zVect(k), 'y')];
%         zP2 = [zP2; generateZSlicePoints(bPiecewise, zVect(k), 'x')];
%         zP2 = [zP2; generateZSlicePoints(cPiecewise, zVect(k), 'y')];    
% 
%     end
try    
    zP = [vertcat(xPiecewise{:}); vertcat(yPiecewise{:}); vertcat(zPiecewise{:})];
    zP2 = [vertcat(bPiecewise{:}); vertcat(cPiecewise{:}); vertcat(dPiecewise{:})];

    %% Un-rotate point cloud
    % Rotate back zP2
    zP2 = rotateCloudAroundAxis(zP2, -rotationAroundZ, 'z');
    zComb = [zP; zP2];
    % Rotate back rest
    zComb = rotateCloudAroundAxis(zComb, -rotationAroundX, 'x');
catch
    assignin('base', 'xPiecewise', xPiecewise);
    assignin('base', 'yPiecewise', yPiecewise);
    assignin('base', 'zPiecewise', zPiecewise);
end
    %% Unskew original point cloud
    % Undo previous rotation to return to original point cloud orientation. 

    % zP = rotateCloudAroundAxis(zP, -rotationAroundX, 'x');

    %% Display
    figure(1)
    clf(1);
    % plot3(pts(:,1), pts(:,2), pts(:,3), '.', 'markersize', 1, 'markeredgecolor', [0.8, 0.8, 0.8]);
    hold on
    % plot(minSet(:,1), minSet(:,2), 'bx');
    % for k = 1:length(bPiecewise)
    %     
    %     if ~isempty(bPiecewise{k, 1})
    %         plot3(bPiecewise{k,2}*ones(size(bPiecewise{k}, 1), 1), ...
    %             bPiecewise{k,1}(:,1), ...
    %             bPiecewise{k,1}(:,2), ...    
    %             'r', 'linewidth', 1);
    %     end
    % end
    % 
    % for k = 1:length(cPiecewise)
    %     
    %     if ~isempty(cPiecewise{k, 1})
    %         plot3(cPiecewise{k,1}(:,1), cPiecewise{k,2}*ones(size(cPiecewise{k}, 1), 1), ...
    %             cPiecewise{k,1}(:,2), ...    
    %             'b', 'linewidth', 1);
    %     end
    % end

    plot3(zComb(:,1), zComb(:,2), zComb(:,3), 'k.', 'markersize', 1)

    % 
    % plot(shp, 'facecolor', [.9, .8, .1], 'edgecolor', [.8, .7 ,.5], 'FaceAlpha', 0.6);
    % 
    hold off
    % set(gca, 'xlim', [0, sizeX], 'ylim', [0, sizeY]);

    % figure(2)
    % plot(slicePts(:,1), slicePts(:,2), '.', 'markersize', 3, 'markeredgecolor', [0.6, 0.6, 0.6])
    % hold on
    % plot(piecewisePoints(:,1), piecewisePoints(:,2), 'r-x', 'linewidth', 2);
    % hold off
    % set(gca, 'xlim', [0, sizeX], 'ylim', [0, sizeY]);

    % %% 
    % 
    % numPoints = 10;
    % 
    % ptCloud = pointCloud(zComb);
    % % Calculate normals
    % normals = pcnormals(ptCloud, numPoints);
    % ptCloud = pointCloud(zComb, 'normal', normals);
    % % Save as ply file
    % pcwrite(ptCloud, 'ptCloud20190527.ply');
    % 


    %% Generate a fitted mesh to the extracted points
    % Can pass off to MeshLab to do the crunching here

    dlmwrite('ptsOut20200618.xyz', zComb, 'delimiter', '\t');

    %%
    % In MeshLab
    % Load file
    % compute normals w/ defaults
    % reconstruct surface w/ Poisson + num pts 15
    % Discrete curvature w/ mean curvature
    % smooth: laplace vertex color x 3

    % To add : call out modifications to script as inputs to this function.
    % 
    scriptFile = modifyMeshLabScript('poissonRecon.mlx', 'saveInPlace', false, 'poissonDepth', 9);

    polygonReturn = processWithMeshLab('ptsOut20200618.xyz', scriptFile); % Works!
    delete(scriptFile);

end

