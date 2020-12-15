
function varargout = meshFitting3DFcn(varargin)
%% Inputs

if nargin == 14
    sizeX = varargin{1};
    sizeY = varargin{2};
    sizeZ = varargin{3};
    nPoints = varargin{4};
    nDetectsPerPoint = varargin{5};
    template = varargin{6};
    xSlice = varargin{7};
    subsampleStep = varargin{8};
    onPoint = varargin{9};
    rotationAroundX = varargin{10};
    rotationAroundZ = varargin{11};
    tempFolder = varargin{12};
    reconScript = varargin{13};
    reconParams = varargin{14};
elseif nargin == 1
    
    meshObj = varargin{1};
    
    % meshfitter object
    sizeX = meshObj.inputs.sizeX;
    sizeY = meshObj.inputs.sizeY;
    sizeZ = meshObj.inputs.sizeZ;
    nPoints = meshObj.inputs.nPoints;
    nDetectsPerPoint = meshObj.inputs.nDetectsPerPoint;
    template = meshObj.inputs.template;
    xSlice = meshObj.inputs.xSlice;
    subsampleStep = meshObj.inputs.subsampleStep;
    onPoint = meshObj.inputs.onPoint;
    rotationAroundX = meshObj.inputs.rotationAroundX;
    rotationAroundZ = meshObj.inputs.rotationAroundZ;
    tempFolder = meshObj.inputs.tempFolder;
    reconScript = meshObj.inputs.scripts.poissonRecon;
    reconParams = meshObj.inputs.scripts.params;

end


%%
                                                    
                                                    
%----------------------------------------------------------%


    if strcmp(template.Format, 'stlFile')
        [xunit, yunit, zunit, groundTruth] = GenerateTestDataFromSTLFile(template.File, nPoints);
        
        ImportData = Import1File(template.exampleDataFile);
        data = ImportData.Data;
        
        % From zero-noise points sampled on STL file, make some proxy SMLM data
        z = pointsToSMLMPointCloud(xunit, yunit, zunit, nPoints, nDetectsPerPoint,...
                                    data(:,7:8));

        % Filter points to remain in ROI bounds
        z = enforceROIBounds(z, sizeX, sizeY, sizeZ);
        
    elseif strcmp(template.Format, 'pointCloud')
        % Directly supplied an experimental point cloud
        % No need to generate from STL file
        
        z = Import1File(template.exampleDataFile);
        
    end



    %% Rotate input point cloud to avoid having one big plane parallel to the z axis
    % First step in 3d reconstruction analysis.  Required if object has a large
    % area that is parallel to one of the spatial planes.  Adding a tilt to the
    % object works around this problem. 

    zRot = rotateCloudAroundAxis(z, rotationAroundX, 'x');

    %% Loop over slices in x
    % Main fitting function along single axis

    ptRange = [max(zRot(:,(1:3))); min(zRot(:,1:3))];
    
    assignin('base', 'zRot', zRot);
    
    xPiecewise = fitPointSetLoop(zRot, ptRange(2,1), ptRange(1,1), xSlice, 2, subsampleStep, onPoint);
    
    
    %% Loop over slices in y
    % Repeat above section with same parameters but on perpendicular axis.
    yPiecewise = fitPointSetLoop(zRot, ptRange(2,2), ptRange(1,2), xSlice, 1, subsampleStep, onPoint);
    
    %% Loop over slices in z
    % Repeat above section with same parameters but on perpendicular axis.
    zPiecewise = fitPointSetLoop(zRot, ptRange(2,3), ptRange(1,3), xSlice, 3, subsampleStep, onPoint);


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
    [threshX, splitDistX] = rosinThreshold(xPiecewise, 'doPlot', false);
    [threshY, splitDistY] = rosinThreshold(yPiecewise, 'doPlot', false);
    [threshZ, splitDistZ] = rosinThreshold(zPiecewise, 'doPlot', false);

    [threshB, splitDistB] = rosinThreshold(bPiecewise, 'doPlot', false);
    [threshC, splitDistC] = rosinThreshold(cPiecewise, 'doPlot', false);
    [threshD, splitDistD] = rosinThreshold(dPiecewise, 'doPlot', false);

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

try    
    zP = [vertcat(xPiecewise{:,1}); vertcat(yPiecewise{:,1}); vertcat(zPiecewise{:,1})];
    zP2 = [vertcat(bPiecewise{:,1}); vertcat(cPiecewise{:,1}); vertcat(dPiecewise{:,1})];

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

%%
%%%%%%%%%%%%%%%%%%
% Probably best spot to try to do segementation on denoised data here
% Throw to DBSCAN
%%%%%%%%%%%%%%%%%%


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

    plot3(z(:,1), z(:,2), z(:,3), '.', 'color', [0.6, 0.6, 0.6], 'markersize', 1);
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

    dlmwrite(fullfile(tempFolder, 'ptsOut20200618.xyz'), zComb, 'delimiter', '\t');

    %%
    % In MeshLab
    % Load file
    % compute normals w/ defaults
    % reconstruct surface w/ Poisson + num pts 15
    % Discrete curvature w/ mean curvature
    % smooth: laplace vertex color x 3

    % To add : call out modifications to script as inputs to this function.
    % 
    scriptFile = modifyMeshLabScript(reconScript, 'saveInPlace', false, 'poissonDepth', reconParams.poissonDepth);

    [polygonReturn, ~, meshProps] = processWithMeshLab(fullfile(tempFolder, 'ptsOut20200618.xyz'), scriptFile); % Works!
    delete(scriptFile);
    
    if nargout == 6
        varargout{1} = z;
        varargout{1} = zComb;
        varargout{1} = polygonReturn;
        varargout{1} = pwOut;
        varargout{1} = groundTruth;
        varargout{1} = meshProps;
    elseif nargout == 1
        meshObj.results.z = z;
        meshObj.results.zComb = zComb;
        meshObj.results.polygonReturn = polygonReturn;
        meshObj.results.pwOut = pwOut;
        meshObj.results.meshProps = meshProps;
        
        meshObj.inputs.groundTruth = groundTruth;
        
        varargout{1} = meshObj;
    end

end

