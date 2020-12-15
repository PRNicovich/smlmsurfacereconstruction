function fileOut = modifyMeshLabScript(fileIn, varargin)
    %
    % Allow base meshlab script to be modified for functionalized execution and
    % testing.
    % Available parameters [defaults]:
    % 'saveInPlace' - boolean - to overwrite prototype script file w/ new. [[false]
    % 
    % Poisson reconstruction parameters:
    % 'normalsK' - integer - neighbors for calc'ing normals. [10]
    % 'normalsSmoothIter' - integer - number of smoothing iterations for normals [0]
    % 'poissonDepth' - integer - reconstruction depth [8]
    % 'poissonFullDepth' - integer - Adaptive octree depth [5]
    % 'poissonIters' - integer - Gauss-Seidel relaxations [8]
    % 'poissonPointWeight' - integer - Interpolation weight [4]
    % 'poissonSamplesPerNode' - numeric - Minimum number of samples per node [15]
    %
    % Curvature calculation parameters:
    % 'curvatureFilterScale' - integer - neighbors for calc'ing curvature. [10]
    % 'curvatureProjectionAccuracy' - float - threshold for curvatur projection. [0.0001]
    % 'curvatureProjectionMaxIters' - integer - max iterations for projection. [15]
    % 'curvatureSphereParam' - float - spherical parameter for calc; 1 is sphere, 0 is plane. [1]
    % 'curvatureType' - integer - Mean (0), Gauss (1), K1 (2), K2 (3), ApproxMean (4) - [0]

    p = inputParser;
    % Poisson reconsruction
    addRequired(p, 'fileIn', @ischar);
    addParameter(p, 'saveInPlace', false, @islogical);
    addParameter(p, 'normalsK', 10, @isnumeric);
    addParameter(p, 'normalsSmoothIter', 0, @isnumeric);
    addParameter(p, 'poissonDepth', 8, @isnumeric);
    addParameter(p, 'poissonFullDepth', 5, @isnumeric);
    addParameter(p, 'poissonIters', 8, @isnumeric);
    addParameter(p, 'poissonPointWeight', 4, @isnumeric);
    addParameter(p, 'poissonSamplesPerNode', 15, @isnumeric);
    
    % Curvature of vertices via APSS
    addParameter(p, 'curvatureFilterScale', 10, @isnumeric);
    addParameter(p, 'curvatureProjectionAccuracy', 0.0001, @isnumeric);
    addParameter(p, 'curvatureProjectionMaxIters', 15, @isnumeric);
    addParameter(p, 'curvatureSphereParam', 1, @isnumeric);
    addParameter(p, 'curvatureType', 0, @isnumeric);
    
    
    parse(p, fileIn, varargin{:});

    if p.Results.saveInPlace
        fileOut = fileIn;
    else
        fileOut = [tempname(pwd), '.mlx'];
    end

    % Read in specified file
    fid = fopen(fileIn, 'r');
    tline = fgetl(fid);
    t = 1;
    inFile{t} = tline;
    while ischar(tline)
       tline = fgetl(fid);
       inFile{t} = tline;
       t = t + 1;
    end
    inFile(end) = [];
    fclose(fid);

    % Change values identified in name, value pairs.  Focused on
    % poisson reconstruction workflow.

    for k = 1:length(p.Parameters)
        if ~ismember(p.Parameters{k}, p.UsingDefaults)

            switch p.Parameters{k}
                case 'normalsK'
                    line = find(~cellfun(@isempty, strfind(inFile, 'name="K"')));

                    [sP, eP] = startEndPositionFind(inFile{line});

                    inFile{line} = ...
                        strcat(inFile{line}(1:sP), num2str(p.Results.normalsK), inFile{line}(eP:end));

                case 'normalsSmoothInter'
                    line = find(~cellfun(@isempty, strfind(inFile, 'name="smoothIter"')));

                    [sP, eP] = startEndPositionFind(inFile{line});

                    inFile{line} = ...
                        strcat(inFile{line}(1:sP), num2str(p.Results.normalsSmoothInter), inFile{line}(eP:end));


                case 'poissonDepth'
                    line = find(~cellfun(@isempty, strfind(inFile, 'name="depth"')));
                    [sP, eP] = startEndPositionFind(inFile{line});
                    inFile{line} = ...
                            strcat(inFile{line}(1:sP), num2str(p.Results.poissonDepth), inFile{line}(eP:end));


                case 'poissonFullDepth'
                    line = find(~cellfun(@isempty, strfind(inFile, 'name="fullDepth"')));
                    [sP, eP] = startEndPositionFind(inFile{line});

                    inFile{line} = ...
                        strcat(inFile{line}(1:sP), num2str(p.Results.poissonFullDepth), inFile{line}(eP:end));


                case 'poissonIters'
                    line = find(~cellfun(@isempty, strfind(inFile, 'name="iters"')));
                    [sP, eP] = startEndPositionFind(inFile{line});

                    inFile{line} = ...
                        strcat(inFile{line}(1:sP), num2str(p.Results.poissonIters), inFile{line}(eP:end));


                case 'poissonPointWeight'
                    line = find(~cellfun(@isempty, strfind(inFile, 'name="pointWeight"')));
                    [sP, eP] = startEndPositionFind(inFile{line});

                    inFile{line} = ...
                        strcat(inFile{line}(1:sP), num2str(p.Results.poissonPointWeight), inFile{line}(eP:end));

                case 'poissonSamplesPerNode'
                    line = find(~cellfun(@isempty, strfind(inFile, 'name="samplesPerNode"')));
                    [sP, eP] = startEndPositionFind(inFile{line});

                    inFile{line} = ...
                        strcat(inFile{line}(1:sP), num2str(p.Results.poissonSamplesPerNode), inFile{line}(eP:end));

                case 'curvatureFilterScale'
                    line = find(~cellfun(@isempty, strfind(inFile, 'name="FilterScale"')));
                    [sP, eP] = startEndPositionFind(inFile{line});

                    inFile{line} = ...
                        strcat(inFile{line}(1:sP), num2str(p.Results.curvatureFilterScale), inFile{line}(eP:end));
                    
                case 'curvatureProjectionAccuracy'
                    line = find(~cellfun(@isempty, strfind(inFile, 'name="ProjectionAccuracy"')));
                    [sP, eP] = startEndPositionFind(inFile{line});

                    inFile{line} = ...
                        strcat(inFile{line}(1:sP), num2str(p.Results.curvatureProjectionAccuracy), inFile{line}(eP:end));
                    
                case 'curvatureProjectionMaxIters'
                    line = find(~cellfun(@isempty, strfind(inFile, 'name="MaxProjectionIters"')));
                    [sP, eP] = startEndPositionFind(inFile{line});

                    inFile{line} = ...
                        strcat(inFile{line}(1:sP), num2str(p.Results.curvatureProjectionMaxIters), inFile{line}(eP:end));
                    
                case 'curvatureSphereParam'
                    line = find(~cellfun(@isempty, strfind(inFile, 'name="SphericalParameter"')));
                    [sP, eP] = startEndPositionFind(inFile{line});

                    inFile{line} = ...
                        strcat(inFile{line}(1:sP), num2str(p.Results.curvatureSphereParam), inFile{line}(eP:end));
                    
                case 'curvatureType'
                    line = find(~cellfun(@isempty, strfind(inFile, 'name="CurvatureType"')));
                    [sP, eP] = startEndPositionFind(inFile{line});

                    inFile{line} = ...
                        strcat(inFile{line}(1:sP), num2str(p.Results.curvatureType), inFile{line}(eP:end));
                

                otherwise
                    % Skip since nothing to set

            end


        else
            % Do nothing - is default value
        end
    end

    fID = fopen(fileOut, 'w+');
    for k = 1:length(inFile)
        fprintf(fID, '%s\n', inFile{k});
    end
    fclose(fID);
end

function [sP, eP] = startEndPositionFind(line)
    % Find start and end positions to insert value in MeshLab filter script
    % file.
    sP = (strfind(line, 'value=')+numel('value='));
    eP = sP + min(strfind(line(sP:end), '"'))+1;
end
