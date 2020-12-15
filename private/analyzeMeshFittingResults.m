function analysis = analyzeMeshFittingResults(inputs, results, analysis)
% Analysis of outputs from meshFitting3DFcn or MeshFitting3D.m demo file.

    refVol = triangulation(inputs.groundTruth.faces, inputs.groundTruth.vertices);
    expVol = triangulation(results.polygonReturn.faces, results.polygonReturn.vertices);

    %%%%%%%
    % Volume of recon vs reference
    % Read out of file for reference
    scriptFile = modifyMeshLabScript(analysis.scripts.objectProperties, 'saveInPlace', false);
    [~, ~, meshProps] = processWithMeshLab(inputs.template.File, scriptFile, inputs.scripts.meshlabPath); % 
    delete(scriptFile);
    % Remember to scale
    analysis.volume.ref = meshProps.MeshVolume;

    % Recall from results.meshProps for expVol
    analysis.volume.exp = results.meshProps.MeshVolume;

    %%%%%%%
    % Per-recon'd-tile error
    % Align first to remove uninteresting errors
    % Decimate so this doesn't take forever
    [~, transMat] = ICP_finite(refVol.Points(1:10:end, :), expVol.Points(1:10:end, :));

    transPoints = pctransform(pointCloud(expVol.Points), affine3d(transMat'));
    analysis.alignedVolume.volume = triangulation(expVol.ConnectivityList, transPoints.Location);
    analysis.alignedVolume.alignedVolFile = fullfile(inputs.tempFolder, 'alignedOutputVol.stl');
    stlwrite(analysis.alignedVolume.volume, analysis.alignedVolume.alignedVolFile);

    %%%%%%
    % Displacement
    % Hausdorff distance
    % P Cignoni, C Rocchini, R Scopigno - Computer Graphics Forum, 1998 - Blackwell Synergy
    % Distance from calc'd + aligned volume -> Reference
    scriptFile = modifyMeshLabScript(analysis.scripts.hausdorff, 'saveInPlace', false);
    inputFiles = {inputs.template.File, analysis.alignedVolume.alignedVolFile};

    [analysis.displacement.alignedToReference, ~] = processWithMeshLab(inputFiles, scriptFile, inputs.scripts.meshlabPath, ...
        {}, fullfile(inputs.tempFolder, 'hausDorfAliToRef.ply')); % 

    % Distance from Reference -> calc'd and aligned volume
    inputFiles = {analysis.alignedVolume.alignedVolFile, inputs.template.File};

    [analysis.displacement.referenceToAligned, ~] = processWithMeshLab(inputFiles, scriptFile, inputs.scripts.meshlabPath, ...
        {}, fullfile(inputs.tempFolder, 'hausDorfRefToAli.ply')); % 
    delete(scriptFile);

    %%%%%%
    % Curvature of each tile <- get this in the processing script
    % Use APSS method to colorize faces, apply to vertices, export here

    scriptFile = modifyMeshLabScript(analysis.scripts.curvature, 'saveInPlace', false);

    [analysis.curvature.aligned, ~] = processWithMeshLab(analysis.alignedVolume.alignedVolFile, ...
                                        scriptFile, inputs.scripts.meshlabPath, ...
                                        {}, fullfile(inputs.tempFolder, 'curvatureAliVol.ply')); % 

    [analysis.curvature.reference, ~] = processWithMeshLab(inputs.template.File, ...
                                        scriptFile, inputs.scripts.meshlabPath, ...
                                        {}, fullfile(inputs.tempFolder, 'curvatureAliVol.ply')); % 

    delete(scriptFile);
end