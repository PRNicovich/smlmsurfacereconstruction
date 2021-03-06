classdef meshfitter < matlab.mixin.SetGet
    % meshfitter class for collecting parameters for meshFitting3DFcn.m,
    % analyzeMeshFittingResults.m
    
    properties
       
        inputs = struct('sizeX', 10000, ... % Size of volume bounding box (nm)
                        'sizeY', 10000, ...
                        'sizeZ', 10000, ...
                        'nPoints', 40000, ... % Number of emitters to generate
                        'nDetectsPerPoint', 12, ... % mean number of detection events per emitter, poisson distributed
                        'template', struct('exampleDataFile', fullfile(pwd, 'testData', 'CF647 KRas_8.txt'), ... % File containing real SMLM data.  Used to generate localization precision
                                           'Format', 'stlFile', ...
                                           'File', fullfile(pwd, 'testData', 'scaledStanfordBunny.stl')), ...
                        'tempFolder', fullfile(pwd, 'output'), ... % temp folder for intermediate files
                        'xSlice', 60, ... % Thickness of single slice, in nm
                        'subsampleStep', 10, ... % Subsample of fitted curve, in points
                        'onPoint', 5, ... % Start curve on this point.  Avoids deviation from center of point clound from using min span tree to find ridge
                        'rotationAroundX', pi/10, ...% Rotation around X axis to avoid bad analysis of area parallel to XY plane; in radians
                        'rotationAroundZ', pi/4, ... % Rotation around Z to re-slice second pass of point cloud fitting; in radians
                        'scripts', struct('poissonRecon', fullfile(pwd, 'meshlabScripts', 'poissonRecon2020.mlx'), ...
                                          'params', struct('poissonDepth', 9), ...
                                          'meshlabPath', 'C:\program files\vcg\meshlab\meshlabserver.exe'), ...
                        'groundTruth', [], ...
                        'dbScan', struct('minPts', 10, ...
                                         'eta', 3));
                    
        %-------------------------------------
        % Results from meshFitting3DFcn 
        %-------------------------------------
                    
        results = struct('z', [], ...
                         'zComb', [], ...
                         'polygonReturn', [], ...
                         'pwOut', [], ...
                         'meshProps', []);
                                          
                 

        %-------------------------------------
        % Analysis parameters
        %-------------------------------------

        analysis = struct('scripts', struct('objectProperties', fullfile(pwd, 'meshlabScripts', 'objectProperties.mlx'), ...
                                            'hausdorff', fullfile(pwd, 'meshlabScripts', 'hausDorff.mlx'), ...
                                            'curvature', fullfile(pwd, 'meshlabScripts', 'curvatureVertices.mlx')), ...
                          'volume', struct('ref', [], ...
                                           'exp', []), ...
                          'alignedVolume', struct('volume', [], ...
                                                  'alignedVolFile', ''), ...
                          'displacement', struct('alignedToReferenced', [], ...
                                                 'referenceToAligned', []), ...
                          'curvature', struct('aligned', [], ...
                                              'reference', []));
                                        
        
    end
    
end