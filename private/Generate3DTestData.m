function pts = GenerateTestDataFromSTLFile(repPALMData, STLFileName, NPointsToGenerate)

% Script to make test data for the 3D PALM analysis script

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pull photon distribution from a specified file
fileName = 'D:\MATLAB\JQ3DPALM\488.txt';

STLFileName = 'D:\Dropbox\Proposals\Discovery - 3D PALM\Figures\ResolutionTestSTL.stl';

NPointsToGenerate = round(5.512*10000);
CombineCoverslipAndObjects = 0; % Are coverslip and together (coverslip is base plane) or separate?

Objective_NA = 1.49; % NA of objective
SampleMediumRI = 1.518;
Wavelength = 680; % Emission wavelength (nm)
PixelSizeinnm = 100;
DarkCnts = 65;

RegionLimits = [-2.5e6 2.5e6]; % In pm

Save_path = 'D:\Dropbox\Proposals\Discovery - 3D PALM\Figures';

%%
%%%%%%%%%%%%%%%%%%%%%%
% SMLM Proxy Data import + generation

ImportData = Import1File(fileName);

data = ImportData.Data;

lNdata = data(:,10);
lnest = lognfit(lNdata);

[a, b] = hist(lNdata, 1000);

a = a/(numel(lNdata)*diff(b(1:2)));

% Make some values for the number of photons from each particle
makeData = lognrnd(lnest(1), lnest(2), NPointsToGenerate,1);
[aa, bb] = hist(makeData, 1000);

aa = aa/(numel(makeData)*diff(bb(1:2)));

figure(1)
plot(b, a);
hold on
plot(b, lognpdf(b, lnest(1), lnest(2)), 'r');
plot(bb, aa, 'g');
hold off
xlabel('N Photons'); ylabel('PDF');
legend({'Data', 'LogNormFit', 'GeneratedData'});

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import 3D surfaces made in CAD, exported as STL file

% Enter here if coverslip and objects modeled as separate STL files
if CombineCoverslipAndObjects == 1
    
    objFV = stlread('D:\MATLAB\3D PALM\JQ3DPALM\TestPatternObjects.stl');
    cslipFV = stlread('D:\MATLAB\3D PALM\JQ3DPALM\TestPatternCoverslip.stl');

    figure(2)
    patch(objFV,'FaceColor',       [0.6 .8 1.0], ...
        'EdgeColor',       'none',        ...
        'FaceLighting',    'gouraud',     ...
        'AmbientStrength', 0.15);
    hold on
    patch(cslipFV,'FaceColor',       [1 0.8 0.6], ...
        'EdgeColor',       'none',        ...
        'FaceLighting',    'gouraud',     ...
        'AmbientStrength', 0.15);
    hold off
    
elseif CombineCoverslipAndObjects == 0

    % With combined object + coverslip STL file
    objFV = stlread(STLFileName);

    figure(2)
    clf
    patch(objFV,'FaceColor',       [0.6 .8 1.0], ...
        'EdgeColor',       'none',        ...
        'FaceLighting',    'gouraud',     ...
        'AmbientStrength', 0.15);

end

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
view([12 12 10]);
lightHand = findobj('parent', gca, 'type', 'light');
for k = 1:numel(lightHand)
    lightangle(lightHand(k), 110, 30)
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate test STORM points
% Sample all faces with likelihood equal to their fraction of the total
% area.  


% Enter here if coverslip and objects modeled as separate STL files
if CombineCoverslipAndObjects == 1
    
    % Make top of coverslip object (y value of 0) as bottom of objects and
    % across whole image area. 
    % Any faces in object that were on the coverslip cannot be selected to
    % have a molecule simulated on them
    topPoints = find(cslipFV.vertices(:,2) == 0);
    cslipTop.faces = cslipFV.faces(all(ismember(cslipFV.faces, topPoints), 2), :);
    cslipTop.vertices = cslipFV.vertices;
   
    % Find all bottom-intersecting points in object
    bottPoints = find(objFV.vertices(:,2) <= 0);
    objBott.faces = objFV.faces(all(ismember(objFV.faces, bottPoints), 2), :);
    objBott.vertices = objFV.vertices;
    
%     figure(3)
%     patch(cslipTop, 'facecolor', [1 .8 .2]);
%     hold on
%     patch(objBott, 'facecolor', [.2 .2 1], 'edgecolor', 'k');
%     hold off
%     view([0 0])
    
    % Make triangulation of non-covered coverslip surface
    exposedCslip = [cslipTop.vertices(cslipTop.faces, :); objBott.vertices(objBott.faces, :)];
    exposedCslip = unique(exposedCslip, 'rows');
    
    
    figure(3)
    plot(exposedCslip(:,1), exposedCslip(:,3), 'o', 'markersize', 2);

    triRep = delaunayTriangulation(exposedCslip(:,1), exposedCslip(:,3));
    
    % Left here, unfinished, and going with simpler plan to just have a
    % single object for objects + coverslip at this point

else
    
    % With combined object, things below the top surface of the coverslip face
    % need to be removed
    lowPoints = find(objFV.vertices(:,2) < 0);
    cslipLower.faces = objFV.faces(all(ismember(objFV.faces, lowPoints), 2), :);
    cslipLower.vertices = objFV.vertices;
    % objFV.vertices(cslipLower.faces, :) = [];
    objFV.faces(all(ismember(objFV.faces, lowPoints), 2), :) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample faces to make 3D SMLM data

% http://www.mathworks.com/matlabcentral/answers/93023-is-there-a-matlab-function-that-can-compute-the-area-of-my-patch
vectA = objFV.vertices(objFV.faces(:, 2), :) - objFV.vertices(objFV.faces(:, 1), :);
vectB = objFV.vertices(objFV.faces(:, 3), :) - objFV.vertices(objFV.faces(:, 1), :);
crossVect = cross(vectA, vectB, 2);
area = 1/2 * (sqrt(sum(crossVect.^2, 2)));

% Randomly select faces to hold points
imagedFaces = randsample(1:numel(area), NPointsToGenerate, true, area);
[faceA, faceB] = histc(imagedFaces, 1:1:numel(area));
uniqueFaces = [sort(unique(imagedFaces))', faceA(faceA > 0)'];

% Within each face, choose a random position for the molecule to be
MolCoords = zeros(NPointsToGenerate, 3);
molCount = 1;

for k = 1:size(uniqueFaces, 1);
    Face = uniqueFaces(k,1);
    TriA = objFV.vertices(objFV.faces(Face, 1), :);
    TriB = objFV.vertices(objFV.faces(Face, 2), :);
    TriC = objFV.vertices(objFV.faces(Face, 3), :);
    randPoints = rand(uniqueFaces(k,2), 2);
    
    % Calculate random point within triangle patch defined by
    % objFV.vertices(objFV.faces(Face, N)) giving coordinates ABC
    % P=(1?r1???)A+(r1???(1?r2))B+(r2r1???)C
    
    samplePoints = (1 - sqrt(randPoints(:,1)))*TriA + (sqrt(randPoints(:,1)).*(1 - randPoints(:,2)))*TriB + (randPoints(:,2).*sqrt(randPoints(:,1)))*TriC;
    
    MolCoords(molCount:molCount+uniqueFaces(k,2)-1, :) = samplePoints;
    
    molCount = molCount + uniqueFaces(k,2);
    
end

% Change dimensions of MolCoords into nm, get dimension order right
MolCoords = MolCoords(:,[1 3 2])*1000;
% MolCoords(:,1:2) = MolCoords(:,1:2) + 12500;

%%
% And blur it out by the uncertainty in x, y, and z

Background_holder = poissrnd(DarkCnts*ones(NPointsToGenerate, 1));
Bkgd_variance = 1e-2*round(var(Background_holder(randi(NPointsToGenerate, NPointsToGenerate, 50)), [], 2)/1e-2);



PrecisionXY = (Thompson((1.22*Wavelength/(Objective_NA*2))*ones(NPointsToGenerate, 1), PixelSizeinnm, makeData(1:NPointsToGenerate), Bkgd_variance));
% PrecisionXY = 0.5*ones(NPointsToGenerate, 1);
% All other parameters match Zeiss data fairly well except for Z precision.
%  This number fudged to come up with something at least reasonable, though
%  data shows values ~20 nm for z precision.  
PrecisionZ = (Thompson(Wavelength*ones(NPointsToGenerate, 1), PixelSizeinnm, makeData(1:NPointsToGenerate), Bkgd_variance));
% PrecisionZ = 0.5*ones(NPointsToGenerate, 1);

% Fudge for vesicle figure - Guess we can do 2x better than done before
% with MicAO
PrecisionXY = PrecisionXY/4;
PrecisionZ = PrecisionZ/4;

for k = 1:size(MolCoords, 1)
    
    MolCoords(k, 1:2) = sqrt(PrecisionXY(k))*randn(1, 2) + MolCoords(k, 1:2);
    MolCoords(k, 3) = sqrt(PrecisionZ(k))*randn(1, 1) + MolCoords(k, 3);
    
end

%%

% Delete any molecular coordinates that go out of XY working area
PrecisionXY(MolCoords(:,1) < RegionLimits(1) | MolCoords(:,2) < RegionLimits(1), :) = [];
PrecisionZ(MolCoords(:,1) < RegionLimits(1) | MolCoords(:,2) < RegionLimits(1), :) = [];
makeData(MolCoords(:,1) < RegionLimits(1) | MolCoords(:,2) < RegionLimits(1), :) = [];
Bkgd_variance(MolCoords(:,1) < RegionLimits(1) | MolCoords(:,2) < RegionLimits(1), :) = [];
MolCoords(MolCoords(:,1) < RegionLimits(1) | MolCoords(:,2) < RegionLimits(1), :) = [];

PrecisionXY(MolCoords(:,1) > RegionLimits(2) | MolCoords(:,2) > RegionLimits(2), :) = [];
PrecisionZ(MolCoords(:,1) > RegionLimits(2) | MolCoords(:,2) > RegionLimits(2), :) = [];
makeData(MolCoords(:,1) > RegionLimits(2) | MolCoords(:,2) > RegionLimits(2), :) = [];
Bkgd_variance(MolCoords(:,1) > RegionLimits(2) | MolCoords(:,2) > RegionLimits(2), :) = [];
MolCoords(MolCoords(:,1) > RegionLimits(2) | MolCoords(:,2) > RegionLimits(2), :) = [];



%%%% Plot it as a big ol' point cloud
figure(4)
plot3(MolCoords(:,1), MolCoords(:,2), MolCoords(:,3), 'bo', 'markersize', 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate output file

% Assemble into 1.txt file

DataMatrix = ones(size(MolCoords, 1), 14);

DataMatrix(:,1) = 1:size(DataMatrix, 1);
% DataMatrix(:,2) = 
% DataMatrix(:,3) = 
% DataMatrix(:,4) = 
DataMatrix(:,5) = MolCoords(:,1); % Position X (nm)
DataMatrix(:,6) = MolCoords(:,2); % Position Y (nm)
DataMatrix(:,7) = MolCoords(:,3) - min(MolCoords(:,3)); % Position Z (nm)
DataMatrix(:,8) = PrecisionXY; % Precision (nm)
DataMatrix(:,9) = PrecisionZ; % PrecisionZ (nm)
DataMatrix(:,10) = makeData'; % Number photons
DataMatrix(:,11) = Bkgd_variance; % Background variance, synthesized above
DataMatrix(:,12) = rand(size(DataMatrix, 1), 1); % Dummy values, Chi square
% DataMatrix(:,13) = 
DataMatrix(:,14) = ones(size(DataMatrix, 1), 1); % Channel
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output to format of Zeiss 1.txt file.

Header_string{1} = 'Index';
Header_string{2} = 'First Frame';
Header_string{3} = 'Number Frames';
Header_string{4} = 'Frames Missing';
Header_string{5} = 'Position X [nm]';
Header_string{6} = 'Position Y [nm]';
Header_string{7} = 'Position Z [nm]';
Header_string{8} = 'Precision [nm]';
Header_string{9} = 'Precision [nm]';
Header_string{10} = 'Number Photons';
Header_string{11} = 'Background variance';
Header_string{12} = 'Chi square';
Header_string{13} = 'PSF width [nm]';
Header_string{14} = 'Channel';

Footer_string{1}{1} = 'VoxelSizeX';
Footer_string{1}{2} = 'VoxelSizeY';
Footer_string{1}{3} = 'ResolutionX';
Footer_string{1}{4} = 'ResolutionY';
Footer_string{1}{5} = 'SizeX';
Footer_string{1}{6} = 'SizeY';

% Always use 25.6 x 25.6 mm (converted to micron) size working area in the CAD program
SizeX = 2560;
SizeY = 2560;

VoxelSizeX = 0.1;
VoxelSizeY = 0.1;
ResolutionX = 0.1;
ResolutionY = 0.1;

Footer_matrix = [VoxelSizeX; VoxelSizeY; ResolutionX; ResolutionY; SizeX; SizeY];

[Footer_string{3}{1}, Footer_string{3}{2}] = deal('um');
[Footer_string{3}{3}, Footer_string{3}{4}, Footer_string{3}{5}, Footer_string{3}{6}] = deal([]);


% Pull save name for 1.txt file out of original STL file name
[~, sNameHold] = fileparts(STLFileName);
save_name = sprintf('%s_%gpoints_MicAO.txt', sNameHold, NPointsToGenerate);


fID = fopen(fullfile(Save_path, save_name), 'w');

% Header
for k = 1:(length(Header_string))
    
    fprintf(fID, '%s\t', Header_string{k});
    
end

fprintf(fID, '\r\n');


% Body
% Remake this with fprintf
fprintf(fID, '%d\t%d\t%d\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.0f\t%.2f\t%.2f\t%.1f\t%d \r\n', DataMatrix');
%dlmwrite(fullfile(Save_path, save_name), Output_matrix, 'delimiter', '\t', 'roffset', 0, 'newline', 'pc', '-append');
fseek(fID, 0, 1);
fprintf(fID, '\r\n');

% Footer
for k = 1:length(Footer_matrix)
    
    fprintf(fID, '%s : %f %s\r\n', Footer_string{1}{k}, Footer_matrix(k,:), Footer_string{3}{k});
    
end

% Add in metadata for generation of this 1.txt file at the end of the
% footer for repetition's sake
fprintf(fID, '%s : %s %s\r\n', 'NPhotonsFile', fileName, []);
fprintf(fID, '%s : %s %s\r\n', 'ModelFile', STLFileName, []);
fprintf(fID, '%s : %g %s\r\n', 'NPoints', NPointsToGenerate, []);
fprintf(fID, '%s : %.3f %s\r\n', 'Objective_NA', Objective_NA, []);
fprintf(fID, '%s : %.3f %s\r\n', 'SampleMediumRI', SampleMediumRI, []);
fprintf(fID, '%s : %.0f %s\r\n', 'Wavelength', Wavelength, []);
fprintf(fID, '%s : %.0f %s\r\n', 'PixelSizeinnm', PixelSizeinnm, []);
fprintf(fID, '%s : %.0f %s\r\n', 'DarkCnts', DarkCnts, []);

fclose(fID);






