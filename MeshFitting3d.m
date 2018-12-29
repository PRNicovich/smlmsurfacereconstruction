
%% Simulate points around ground truth

% Double-check units here
sizeX = 10000; % nm
sizeY = 10000; % nm
sizeZ = 10000; % nm

nPoints = 40000;
nDetectsPerPoint = 12; % mean, poisson distributed

exampleDataFile = 'C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\testData\CF647 KRas_8.txt';

ImportData = Import1File(exampleDataFile);
data = ImportData.Data;

% Make a sphere
circRadius = 700; % nm
circCenter = [5000, 5000, 5000]; % nm
% [xunit, yunit, zunit] = sampleSpherePoints(nPoints, circRadius, circCenter);

% Make a teapot
% [xunit, yunit, zunit] = GenerateTestDataFromSTLFile('C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\classic-teapot\TEAPOT.stl', ...
%     nPoints, circCenter, .07);

[xunit, yunit, zunit] = GenerateTestDataFromSTLFile('C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\High_Resolution_Stanford_Bunny\FLATFOOT_StanfordBunny_jmil_HIGH_RES_Smoothed.stl', ...
    nPoints, circCenter, .08);


nPtsActual = repmat(nDetectsPerPoint, [nPoints, 1]);
nPtsActual = poissrnd(nPtsActual);

pts = zeros(sum(nPtsActual(:)), 5); %[x, y, z, nPhotons, locPrecision]
ptsNow = 1;
for k = 1:size(nPtsActual, 1)
    
    rSamp = randsample(numel(data(:,7)), nPtsActual(k), 'true');
    precHere = data(rSamp, 8).^2;
    
    pts(ptsNow:(ptsNow + nPtsActual(k)-1),:) = ...
        [xunit(k) + precHere.*randn(nPtsActual(k), 1), ...
         yunit(k) + precHere.*randn(nPtsActual(k), 1), ...
         zunit(k) + precHere.*randn(nPtsActual(k), 1), ...
         data(rSamp, 8), ...
         precHere] ;
     
     ptsNow = ptsNow + nPtsActual(k);
    
end

pts(pts(:,1) > sizeX | pts(:,2) > sizeY | pts(:,3) > sizeZ, :) =[ ];
pts(pts(:,1) < 0 | pts(:,2) < 0 | pts(:,3) < 0, :) =[ ];

%% Rotate input point cloud to avoid having one big plane parallel to the z axis

rotationAroundX = pi/10; % radians

X = pts(:,1);
Y = pts(:,2)*cos(rotationAroundX) - pts(:,3)*sin(rotationAroundX);
Z = pts(:,2)*sin(rotationAroundX) + pts(:,3)*cos(rotationAroundX);

pts = [X, Y, Z, pts(:,4:5)];

%% Loop over slices in x

xSlice = 60; % Thickness of single slice, in nm
subsampleStep = 10; % Subsample of fitted curve, in points
onPoint = 5; % Start curve on this point.  Avoids deviation from center of point clound from using min span tree to find ridge

ptRange = [max(pts(:,(1:3))); min(pts(:,1:3))];

xPiecewise = fitPointSetLoop(pts, ptRange(2,1), ptRange(1,1), xSlice, 2, subsampleStep, onPoint);
    
%% Loop over slices in y
yPiecewise = fitPointSetLoop(pts, ptRange(2,2), ptRange(1,2), xSlice, 1, subsampleStep, onPoint);

%% Loop over slices between X and Y axes

rotationAroundZ = pi/4; % radians

X = pts(:,1)*cos(rotationAroundZ) - pts(:,2)*sin(rotationAroundZ);
Y = pts(:,1)*sin(rotationAroundZ) + pts(:,2)*cos(rotationAroundZ);
Z = pts(:,3);

rotPts = [X, Y, Z, pts(:,4:5)];

ptRange = [max(rotPts(:,(1:3))); min(rotPts(:,1:3))];
bPiecewise = fitPointSetLoop(rotPts, ptRange(2,1), ptRange(1,1), xSlice, 2, subsampleStep, onPoint);
cPiecewise = fitPointSetLoop(rotPts, ptRange(2,2), ptRange(1,2), xSlice, 1, subsampleStep, onPoint);



%% Double-check for erroneously long links and split if needed
% Determine threshold for splitting
% Can feed eithe xPiecewise or yPiecewise values
% Do both and then take mean as best choice
[threshX, splitDistX] = rosinThreshold(xPiecewise, 'doPlot', true);
[threshY, splitDistY] = rosinThreshold(yPiecewise, 'doPlot', true);

[threshB, splitDistB] = rosinThreshold(bPiecewise, 'doPlot', true);
[threshC, splitDistC] = rosinThreshold(cPiecewise, 'doPlot', true);

threshVal = (threshX + threshY + threshB + threshC)/4;


%%

% Split up cells at bad links
% Generates extra points that may or may not matter. 
% Def an error, but don't care enough to find out right now.


xPiecewise = splitByLargeDistance(xPiecewise, splitDistX, threshVal);
yPiecewise = splitByLargeDistance(yPiecewise, splitDistY, threshVal);
bPiecewise = splitByLargeDistance(bPiecewise, splitDistB, threshVal);
cPiecewise = splitByLargeDistance(cPiecewise, splitDistC, threshVal);



%% Z slice fitted curves
% Generate regularly-spaced points that hopefully don't have too much trash
% in the set

% Make all connected point pairs into a single variable
% No way this isn't done better in a graph class, but don't know how to
% implement that....

zSlice = xSlice/2;

zVect = ptRange(2,3) : zSlice : ptRange(1,3);

zPiecewise = cell(numel(zVect), 2);
zP = [];
zP2 = [];

for z = 1:numel(zVect)

    
    for x = 1:length(xPiecewise)
        
       if ~isempty(xPiecewise{x, 1})
        
           transPts = find(diff(xPiecewise{x, 1}(:,2) > zVect(z)) ~= 0);

           if ~isempty(transPts)

               for m = 1:numel(transPts)

                   oneSide = xPiecewise{x, 1}(transPts(m), [3, 1, 2]);
                   otherSide = xPiecewise{x, 1}(transPts(m)+1, [3, 1, 2]);

                   rat = 1 - ((oneSide(3) - zVect(z))/(oneSide(3) - otherSide(3)));

                   zP = [zP; rat*(oneSide - otherSide) + otherSide];

               end

           end
           
       end
        
    end
    
    
    for y = 1:length(yPiecewise)

        if ~isempty(yPiecewise{y, 1})
            
            transPts = find(diff(yPiecewise{y, 1}(:,2) > zVect(z)) ~= 0);

            if ~isempty(transPts)

               for m = 1:numel(transPts) 

                   oneSide = yPiecewise{y, 1}(transPts(m), [1, 3, 2]);
                   otherSide = yPiecewise{y, 1}(transPts(m)+1, [1, 3, 2]);

                   rat = 1 - ((oneSide(3) - zVect(z))/(oneSide(3) - otherSide(3)));

                   zP = [zP; rat*(oneSide - otherSide) + otherSide];

               end

            end
            
        end

    end
    
    for b = 1:length(bPiecewise)
        
       if ~isempty(bPiecewise{b, 1})
        
           transPts = find(diff(bPiecewise{b, 1}(:,2) > zVect(z)) ~= 0);

           if ~isempty(transPts)

               for m = 1:numel(transPts)

                   oneSide = bPiecewise{b, 1}(transPts(m), [3, 1, 2]);
                   otherSide = bPiecewise{b, 1}(transPts(m)+1, [3, 1, 2]);

                   rat = 1 - ((oneSide(3) - zVect(z))/(oneSide(3) - otherSide(3)));

                   zP2 = [zP2; rat*(oneSide - otherSide) + otherSide];

               end

           end
           
       end
        
    end
    
	for c = 1:length(cPiecewise)

        if ~isempty(cPiecewise{c, 1})
            
            transPts = find(diff(cPiecewise{c, 1}(:,2) > zVect(z)) ~= 0);

            if ~isempty(transPts)

               for m = 1:numel(transPts) 

                   oneSide = cPiecewise{c, 1}(transPts(m), [1, 3, 2]);
                   otherSide = cPiecewise{c, 1}(transPts(m)+1, [1, 3, 2]);

                   rat = 1 - ((oneSide(3) - zVect(z))/(oneSide(3) - otherSide(3)));

                   zP2 = [zP2; rat*(oneSide - otherSide) + otherSide];

               end

            end
            
        end

    end
    
    
    
end

%% Rotate back zP2

X = zP2(:,1)*cos(-rotationAroundZ) - zP2(:,2)*sin(-rotationAroundZ);
Y = zP2(:,1)*sin(-rotationAroundZ) + zP2(:,2)*cos(-rotationAroundZ);
Z = zP2(:,3);

zP2 = [X, Y, Z];


%% Triangulate fitted points
% xlist = vertcat(xPiecewise{:,1});
% ylist = vertcat(yPiecewise{:,1}); 
% 
% fitList = unique([xlist(:,[3, 1, 2]); ylist(:,[1, 3, 2])], 'rows');
% 
% dt = delaunayTriangulation(fitList(:,1), fitList(:,2), fitList(:,3));    


% shp = alphaShape(zP, 80);

%% Un-rotate point cloud

X = zP(:,1);
Y = zP(:,2)*cos(-rotationAroundX) - zP(:,3)*sin(-rotationAroundX);
Z = zP(:,2)*sin(-rotationAroundX) + zP(:,3)*cos(-rotationAroundX);

zP = [X, Y, Z];

X = zP2(:,1);
Y = zP2(:,2)*cos(-rotationAroundX) - zP2(:,3)*sin(-rotationAroundX);
Z = zP2(:,2)*sin(-rotationAroundX) + zP2(:,3)*cos(-rotationAroundX);

zP2 = [X, Y, Z];

zP = [zP; zP2];

%% Unskew original point cloud

X = pts(:,1);
Y = pts(:,2)*cos(-rotationAroundX) - pts(:,3)*sin(-rotationAroundX);
Z = pts(:,2)*sin(-rotationAroundX) + pts(:,3)*cos(-rotationAroundX);

pts = [X, Y, Z, pts(:,4:5)];

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

plot3(zP(:,1), zP(:,2), zP(:,3), 'k.')

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

%% Generate a fitted mesh to the extracted points
% Can pass off to MeshLab to do the crunching here

dlmwrite('ptsOut20181228_02.xyz', zP, 'delimiter', '\t');

% In MeshLab
% Load file
% compute normals w/ defaults
% reconstruct surface w/ Poisson + num pts 15
% Discrete curvature w/ mean curvature
% smooth: laplace vertex color x 3

