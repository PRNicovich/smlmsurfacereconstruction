% FilamentGraphFitting
% Sandbox for testing fitting of filaments in same graph-based fitting as
% 3D section-wise method

FieldSize = 1000; % in nm, for simulations
RandFieldPoints = 5e2;
PeakSigma = [40 60]; % nm
PtsPerPeak = [0 100]; 
NFilaments = 40;
FilamentCoverage = 0.9;
OversampleImage = 4;
PtsPerFilSpot = [2 6];
FilSpotSigma = [10 20]; % nm
SMLMImageSize = 5000; % in nm
PSFWidth = 100; % Widefield PSF width
widefieldImgSize = 75; % Pixels

%%
% Generate test data


% Make filament image
% Draw a bunch of polynomial curves, then rotate a few of them slightly
% Overlay should look something like microtubules

baseFig = false(FieldSize*OversampleImage, FieldSize*OversampleImage);

figure(1)
clf(1)
hold on

polyDom = linspace(0, pi, FieldSize*OversampleImage) - pi/2;

for k = 1:NFilaments
   
    poly = rand(1,1)*sin((pi + rand(1,1))+5*polyDom*(rand(1,1))) + polyDom + pi*rand(1,1);
    
    rotTheta = pi*rand(1,1);
    rotMat = [cos(rotTheta) -sin(rotTheta); sin(rotTheta), cos(rotTheta)];
    
    poly = rotMat*[poly' polyDom']';
    poly = poly';
    
    domPix = OversampleImage*FieldSize*(poly(:,1) - min(poly(:,1)))/(max(poly(:,1)) - min(poly(:,1)));
    domPoly = OversampleImage*FieldSize*(poly(:,2) - min(poly(:,2)))/(max(poly(:,2)) - min(poly(:,2)));
    domPix(1) = [];
    domPoly(1) = [];
    
    domPoly(domPix < 1) = [];
    domPix(domPix < 1) = [];

    domPix(domPoly < 1) = [];
    domPoly(domPoly < 1) = [];
    
    for m = 1:numel(domPix)
        baseFig(round(domPix(m)), round(domPoly(m))) = baseFig(round(domPix(m)), round(domPoly(m))) + 1;
    end
    
end

baseFig = baseFig(0.5*size(baseFig,1):0.6*size(baseFig,1), 0.5*size(baseFig,2):0.6*size(baseFig,2));


% % Bin down to OversampleImage size
% % Downsample baseImg to dsImg
% dsImg = zeros(FieldSize, FieldSize);
% for p = 2:FieldSize
%     for m = 2:FieldSize
%         
%         kernHere = baseFig(((p-1)*OversampleImage+1):((p*OversampleImage)), ...
%             ((m-1)*OversampleImage+1):((m*OversampleImage)));
%         dsImg(p, m) = sum(kernHere(:));
%         
%     end
% end
snapnow;





%%
skelImg = bwmorph(bwmorph(baseFig, 'dilate', 2), 'erode', 1);
skelImg = bwmorph(skelImg, 'thin', Inf);

crossPoints = bwmorph(skelImg, 'branchpoints');
[cpY, cpX] = find(crossPoints);
% 
imagesc(baseFig);
hold on
for k = 1:numel(cpX)
    plot(cpX, cpY, 'rx');
end
hold off

% From idealized image above, make a SMLM image + widefield image

[c1x, c1y] = find(baseFig);
% Stretch to SMLMImageSize 
c1x = c1x*SMLMImageSize/size(baseFig,1) + rand(1,1);
c1y = c1y*SMLMImageSize/size(baseFig,2) + rand(1,1);

cwhich = randsample(1:numel(c1x), round(numel(c1x)*FilamentCoverage));

chan1 = [c1x(cwhich), c1y(cwhich)];

cpY = cpY*SMLMImageSize/size(baseFig,1) + rand(1,1);
cpX = cpX*SMLMImageSize/size(baseFig,2) + rand(1,1);
chan2 = [cpY, cpX];

% Widefield image from chan1 and chan2 data points
wfImgDom = linspace(0, SMLMImageSize, widefieldImgSize);
wfImgChan1 = hist2(chan1(:,1), chan1(:,2), wfImgDom, wfImgDom);
wfImgChan2 = hist2(chan2(:,1), chan2(:,2), wfImgDom, wfImgDom);

PSFInPixels = PSFWidth/(SMLMImageSize/widefieldImgSize);
PSFMatrix = fspecial('gaussian', 71, PSFInPixels); % Ensure second value here is odd

wfImgChan1 = imfilter(wfImgChan1, PSFMatrix,  'symmetric', 'same');
wfImgChan2 = imfilter(wfImgChan2, PSFMatrix,  'symmetric', 'same');

wfImgChan1 = 4095*(wfImgChan1 - min(wfImgChan1(:)))/(max(wfImgChan1(:)) - min(wfImgChan1(:)));
wfImgChan1 = wfImgChan1 + 200;
wfImgChan1 = poissrnd(wfImgChan1);

wfImgChan2 = 4095*(wfImgChan2 - min(wfImgChan2(:)))/(max(wfImgChan2(:)) - min(wfImgChan2(:)));
wfImgChan2 = wfImgChan2 + 200;
wfImgChan2 = poissrnd(wfImgChan2);


figure(2)
wfImg = imfuseSpecColors((wfImgChan1/max(wfImgChan1(:))), (wfImgChan2/max(wfImgChan2(:))), rgb(39, 174, 96), rgb(142, 68, 173), 'white');
image(wfImg);
axis image
set(gca, 'xtick', [], 'ytick', []);
set(gcf, 'color', [1 1 1]);

% SMLM image from same data set

chan1Pts = [];
for k = 1:size(chan1, 1)
   
    pS = diff(FilSpotSigma)*rand(1)+min(FilSpotSigma);
    pN = round(diff(PtsPerFilSpot)*rand(1)+min(PtsPerFilSpot));
    chan1Pts = [chan1Pts; pS*randn(pN, 2)+repmat(chan1(k,:), pN, 1)];
    
    
end

randPts = SMLMImageSize*rand(RandFieldPoints, 2);
chan1Pts = [chan1Pts; randPts];

chan2Pts = [];
for k = 1:size(chan2, 1)
   
    pS = diff(PeakSigma)*rand(1)+min(PeakSigma);
    pN = round(diff(PtsPerPeak)*rand(1)+min(PtsPerPeak));
    chan2Pts = [chan2Pts; pS*randn(pN, 2)+repmat(chan2(k,:), pN, 1)];
    
    
end

randPts = SMLMImageSize*rand(RandFieldPoints/2, 2);
chan2Pts = [chan2Pts; randPts];

SRdom = linspace(1, SMLMImageSize, SMLMImageSize);
chan1SRhist = hist2(chan1Pts(:,1), chan1Pts(:,2), SRdom, SRdom);
chan2SRhist = hist2(chan2Pts(:,1), chan2Pts(:,2), SRdom, SRdom);

PSFInPixels = mean(FilSpotSigma);
PSFMatrix = fspecial('gaussian', 71, PSFInPixels); % Ensure second value here is odd
chan1SRImg = imfilter(chan1SRhist, PSFMatrix,  'symmetric', 'same');
chan2SRImg = imfilter(chan2SRhist, PSFMatrix,  'symmetric', 'same');

%%

figure(3)
srImg = imfuseSpecColors((chan1SRImg/max(chan1SRImg(:))), (chan2SRImg/max(chan2SRImg(:))), rgb(39, 174, 96), rgb(142, 68, 173), 'white');
image(srImg);
axis image
set(gca, 'xtick', [], 'ytick', []);
set(gcf, 'color', [1 1 1]);

figure(31)
% plot(chan1Pts(:,2), chan1Pts(:,1), 'o', 'markerfacecolor', rgb(39, 174, 96), 'markeredgecolor', 'none', 'markersize', 2);
scatter(chan1Pts(:,2), chan1Pts(:,1),10,'filled', ...
       'MarkerFaceAlpha',4/8,'MarkerFaceColor',rgb(39, 174, 96)) 
hold on
% plot(chan2Pts(:,2), chan2Pts(:,1), 'o', 'markerfacecolor', rgb(142, 68, 173), 'markeredgecolor', 'none', 'markersize', 2);
scatter(chan2Pts(:,2), chan2Pts(:,1),10,'filled', ...
       'MarkerFaceAlpha',4/8,'MarkerFaceColor',rgb(142, 68, 173)) 
hold off
axis image
set(gca, 'ydir', 'reverse');
set(gca, 'xtick', [], 'ytick', []);
set(gcf, 'color', [1 1 1]);
box on

%%
% Pass to 2D graph tracing to find shortest connecting path, backbone
% tracing
% piecewisePoints = meshfitAndSmooth2d(chan1Pts, 10, 5);

slicePts = chan1Pts;

    dt = delaunayTriangulation(slicePts(:,1:2));

    connList = zeros(size(dt.ConnectivityList, 1), 2);
    
    if isempty(connList)
        smoothedFit = [];
        return
    end
    
    i = 1;
    for k = 1:size(dt.ConnectivityList, 1)

        connList(i:i+1,:) = [dt.ConnectivityList(k,1:2); ... 
                             dt.ConnectivityList(k,[1,3])];

        i = i+2;
    end

    [~, idxConn] = unique(connList, 'rows');

    startPts = dt.Points(connList(:,1),:);
    endPts = dt.Points(connList(:,2), :);

    connLength = sqrt((startPts(:,1) - endPts(:,1)).^2 + (startPts(:,2) - endPts(:,2)).^2);

    DG = sparse(connList(idxConn,1), connList(idxConn,2), connLength(idxConn), ...
        max(max(connList(idxConn, :))), max(max(connList(idxConn, :))));
    
    
%     DG = sparse(connList(:,1), connList(:,2), connLength);
    UG = tril(DG + DG');
    [ST,~] = graphminspantree(UG);

    
    %% Find min traverse of spanning tree

    [minA, minB, minWts] = find(ST);
    minGraph = graph(minA, minB, minWts);
    pairDist = distances(minGraph);
    [maxDistA, maxDistB] = find(pairDist == max(pairDist(:)));

    % Choose first point returned for max distance above to remove possible
    % ambiguity
    shortSpan = shortestpath(minGraph, maxDistA(1), maxDistB(1));
    
    minSet = dt.Points(shortSpan, :);
    
    %%
    nSamples = 1e4;
    sampleVector = reshape(randsample(size(pairDist, 1), 1e4, 'true'), [], 2);
    mA = sampleVector(:,1);
    mB = sampleVector(:,2);
%     
%     [mA, mB] = find(pairDist > prctile(pairDist(:), 99.99));
    
    testMinSet = cell(numel(mA), 1);
    for k = 1:numel(mA)
        
       testMinSet{k} = shortestpath(minGraph, mA(k), mB(k)); 
        
    end
    
    testMinSpots = horzcat(testMinSet{:})';

        
%% Plot minspantree
    
    tS = histc(testMinSpots, 1:max(testMinSpots)); 
    ff = find(tS > 10);
    
    ms = dt.Points(ff, :);
   

figure(1)
plot(chan1Pts(:,1), chan1Pts(:,2), 'o', 'markersize', 3, 'markeredgecolor', [0.8, 0.8, 0.8]);
hold on
plot(ms(:,1), ms(:,2), 'b.', 'linewidth', 2);
hold off
set(gca, 'xlim', [0, SMLMImageSize], 'ylim', [0, SMLMImageSize]);

%%

% Break excessively large jumps



% Does this deal with branches at all?



% Make an attempt to trace individual actins w/ enforcing straightness