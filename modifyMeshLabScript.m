function fileOut = modifyMeshLabScript(fileIn, varargin)

% Allow base meshlab script to be modified for functionalized execution and
% testing.
% Available parameters [defaults]:
% 'saveInPlace' - boolean - to overwrite prototype script file w/ new. [[false]
% 'normalsK' - integer - neighbors for calc'ing normals. [10]
% 'normalsSmoothIter' - integer - number of smoothing iterations for normals [0]
% 'poissonDepth' - integer - reconstruction depth [8]
% 'poissonFullDepth' - integer - Adaptive octree depth [5]
% 'poissonIters' - integer - Gauss-Seidel relaxations [8]
% 'poissonPointWeight' - integer - Interpolation weight [4]
% 'poissonSamplesPerNode' - numeric - Minimum number of samples per node [15]


p = inputParser;
addRequired(p, 'fileIn', @ischar);
addParameter(p, 'saveInPlace', false, @islogical);
addParameter(p, 'normalsK', 10, @isnumeric);
addParameter(p, 'normalsSmoothIter', 0, @isnumeric);
addParameter(p, 'poissonDepth', 8, @isnumeric);
addParameter(p, 'poissonFullDepth', 5, @isnumeric);
addParameter(p, 'poissonIters', 8, @isnumeric);
addParameter(p, 'poissonPointWeight', 4, @isnumeric);
addParameter(p, 'poissonSamplesPerNode', 15, @isnumeric);
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

% Change values identified in name, value pairs.  Obviously focused on
% poisson reconstruction workflow.

for k = 1:length(p.Parameters)
    if ~ismember(p.Parameters{k}, p.UsingDefaults)
        
        switch p.Parameters{k}
            case 'normalsK'
                line = find(~cellfun(@isempty, strfind(inFile, 'name="K"')));
                inFile{line} = ...
                        strcat(inFile{line}(1:(strfind(inFile{line}, 'value=')+numel('value='))), ...
                                num2str(p.Results.normalsK), '"/>');

            case 'normalsSmoothInter'
                line = find(~cellfun(@isempty, strfind(inFile, 'name="smoothIter"')));
                inFile{line} = ...
                        strcat(inFile{line}(1:(strfind(inFile{line}, 'value=')+numel('value='))), ...
                                num2str(p.Results.normalsSmoothInter), '"/>');


            case 'poissonDepth'
                line = find(~cellfun(@isempty, strfind(inFile, 'name="depth"')));
                inFile{line} = ...
                        strcat(inFile{line}(1:(strfind(inFile{line}, 'value=')+numel('value='))), ...
                                num2str(p.Results.poissonDepth), '"/>');

            case 'poissonFullDepth'
                line = find(~cellfun(@isempty, strfind(inFile, 'name="fullDepth"')));
                inFile{line} = ...
                        strcat(inFile{line}(1:(strfind(inFile{line}, 'value=')+numel('value='))), ...
                                num2str(p.Results.poissonFullDepth), '"/>');

            case 'poissonIters'
                line = find(~cellfun(@isempty, strfind(inFile, 'name="iters"')));
                inFile{line} = ...
                        strcat(inFile{line}(1:(strfind(inFile{line}, 'value=')+numel('value='))), ...
                                num2str(p.Results.poissonIters), '"/>');

            case 'poissonPointWeight'
                line = find(~cellfun(@isempty, strfind(inFile, 'name="pointWeight"')));
                inFile{line} = ...
                        strcat(inFile{line}(1:(strfind(inFile{line}, 'value=')+numel('value='))), ...
                                num2str(p.Results.poissonPointWeight), '"/>');

            case 'poissonSamplesPerNode'
                line = find(~cellfun(@isempty, strfind(inFile, 'name="samplesPerNode"')));
                inFile{line} = ...
                        strcat(inFile{line}(1:(strfind(inFile{line}, 'value=')+numel('value='))), ...
                                num2str(p.Results.poissonSamplesPerNode), '"/>');

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