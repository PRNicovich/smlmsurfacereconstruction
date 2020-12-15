
function polyOut = readSavedMeshReconstruction(inputFile, varargin)

% Read in ply file with extra parameters. Output with these properties
% included in output structure

    % Vertices and Connectivity in same file, just appended after vertex list.
    fID = fopen(inputFile);
    % Load file
    allData = textscan(fID,'%s','Delimiter','\n');
    fclose(fID);

    plyFile.headerLines = [find(strcmp(allData{1}, 'ply')), find(strcmp(allData{1}, 'end_header'))];
    plyFile.vertexElements = str2double(allData{1}{~cellfun(@isempty, strfind(allData{1}, 'element vertex'))}(2+numel('element vertex'):end));
    plyFile.faceElements = str2double(allData{1}{~cellfun(@isempty, strfind(allData{1}, 'element face'))}(2+numel('element face'):end));


    % Parse headers to get fields and order
    % Not perfect or general. 
    % Tested to pull out info when there's a vertex quality column
    plyFile.headerText = allData{1}(plyFile.headerLines(1):plyFile.headerLines(2));
    vertexStart = find(~cellfun(@isempty, strfind(plyFile.headerText, 'element vertex')))+1;
    vertexEnd = find(~cellfun(@isempty, strfind(plyFile.headerText, 'element face'))) - 1;
    for k = vertexStart:vertexEnd
        parseHeader = textscan(plyFile.headerText{k}, '%s %s %s');
        plyFile.vertexColLabels{k - vertexStart + 1} = parseHeader{end}{1};
    end
    
    faceStart = find(~cellfun(@isempty, strfind(plyFile.headerText, 'element face')))+1;
    faceEnd = find(~cellfun(@isempty, strfind(plyFile.headerText, 'end_header'))) - 1;
    for k = faceStart:faceEnd
        parseHeader = textscan(plyFile.headerText{k}, '%s %s %s');
        plyFile.faceColLabels{k - faceStart + 1} = parseHeader{end}{1};
    end
    
    
    plyFile.nVertexCols = numel(strfind(allData{1}{plyFile.headerLines(2) + 1}, sprintf(' ')));
    plyFile.nFaceCols = numel(strfind(allData{1}{plyFile.headerLines(2) + 1 + plyFile.vertexElements}, sprintf(' ')));
    
    % Read data from file
    plyFile.vertices = dlmread(inputFile, ' ', ...
        [plyFile.headerLines(2)+1, ...
        0, ...
        plyFile.vertexElements + plyFile.headerLines(2) - 1, ...
        plyFile.nVertexCols - 1]);

    plyFile.faces = dlmread(inputFile, ' ', ...
        [plyFile.headerLines(2)+plyFile.vertexElements, ...
        0, ...
        plyFile.vertexElements + plyFile.headerLines(2) + plyFile.faceElements - 1, ...
        plyFile.nFaceCols - 1]);


    % Vertices should only have 'x', 'y', 'z' in 'vertices' field.  
    % Move everything else to a 'vertexProperty' field.
    
    whichCols = zeros(1, plyFile.nVertexCols);
    for k = 1:plyFile.nVertexCols
       whichCols(k) = ismember(plyFile.vertexColLabels{k}, {'x', 'y', 'z'}); 
    end

    plyFile.vertexProperties = plyFile.vertices(:, ~whichCols);
    plyFile.vertices(:, ~whichCols) = [];
    
    polyOut = plyFile;


