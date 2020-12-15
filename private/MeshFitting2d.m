%% Generate test pattern ground truth
% Drawing a circle

circRadius = 1000; % nm
circCenter = [1200,100]; % nm

%% Simulate points around ground truth

sizeX = 2560; % nm
sizeY = 2560; % nm

nPoints = 2000;
nDetectsPerPoint = 2; % mean, poisson distributed

nRandPoints = 0;

precMod = 5; % Multiplier to precision here

freqMultList = [1.0, 1.5, 2.0];


%%
exampleDataFile = fullfile(pwd, 'testData', 'CF647 KRas_8.txt');

ImportData = Import1File(exampleDataFile);
data = ImportData.Data;

    
for fM = 1:length(freqMultList)
    
    freqMult = freqMultList(fM);
    
    % th = rand(nPoints, 1)*6*pi;
    % xunit = circRadius * cos(th) + circCenter(1);
    % yunit = circRadius * sin(th) + circCenter(2);
    th = rand(nPoints, 1)*4*pi;
    yunit = circRadius * cos(th*freqMult) + circCenter(1);
    xunit = 200*th+ circCenter(2);

    calcPeriod = 200*pi*2/freqMult;

    thBase = 0:0.05:4*pi;
    xBase = circRadius * cos(thBase*freqMult) + circCenter(1);
    yBase = 200*thBase+ circCenter(2);

    xunit = [xunit; rand(nRandPoints,1)*sizeX];
    yunit = [yunit; rand(nRandPoints,1)*sizeY];

    nPtsActual = repmat(nDetectsPerPoint, [nPoints+nRandPoints, 1]);
    nPtsActual = poissrnd(nPtsActual);

    pts = zeros(sum(nPtsActual(:)), 4); %[x, y, nPhotons, locPrecision]
    ptsNow = 1;
    for k = 1:size(nPtsActual, 1)

        rSamp = randsample(numel(data(:,7)), nPtsActual(k), 'true');
        precHere = data(rSamp, 8)*precMod;

        pts(ptsNow:(ptsNow + nPtsActual(k)-1),:) = ...
            [xunit(k) + precHere.*randn(nPtsActual(k), 1), ...
             yunit(k) + precHere.*randn(nPtsActual(k), 1), ...
             data(rSamp, 8), ...
             precHere] ;

         ptsNow = ptsNow + nPtsActual(k);

    end

    pts(pts(:,1) > sizeX | pts(:,2) > sizeY, :) =[ ];
    pts(pts(:,1) < 0 | pts(:,2) < 0, :) =[ ];

    %% Render points into image
    renderPixSize = 100;

    % rendImg = zeros(numel(0:renderPixSize:sizeX), ...
    %                 numel(0:renderPixSize:sizeY));
    %             
    % for k = 1:size(pts, 1)
    %    
    %     rendImg = rendImg + gauss_2D([pts(k, 3), ...
    %                                   pts(k,1:2), ...
    %                                   pts(k,4), ...
    %                                   pts(k,4), ...
    %                                   0], ...
    %                                   0:renderPixSize:sizeX);  
    %     
    % end

    histImg = hist3(pts(:,1:2), 'Edges', {0:renderPixSize:sizeX, 0:renderPixSize:sizeY});
    % [cont, chand] = contour(histImg, 10);
    % 
    % slmsettings =slmset('weights', pts(:,3));
    % [slm, xp, yp] = slmengine(pts(:,1), pts(:,2), slmsettings);

    %% Minimum spanning tree of points
    dt = delaunayTriangulation(pts(:,1:2));

    connList = zeros(size(dt.ConnectivityList, 1), 2);
    i = 1;
    for k = 1:size(dt.ConnectivityList, 1)

        connList(i:i+1,:) = [dt.ConnectivityList(k,1:2); ... 
                             dt.ConnectivityList(k,[1,3])];

        i = i+2;
    end

    [unConn, idxConn] = unique(connList, 'rows');

    startPts = dt.Points(connList(:,1),:);
    endPts = dt.Points(connList(:,2), :);

    connLength = sqrt((startPts(:,1) - endPts(:,1)).^2 + (startPts(:,2) - endPts(:,2)).^2);

    DG = sparse(connList(idxConn,1), connList(idxConn,2), connLength(idxConn), ...
        max(max(connList(idxConn, :))), max(max(connList(idxConn, :))));

    % DG = sparse(connList(idxConn,1), connList(idxConn,2), connLength(idxConn));
    UG = tril(DG + DG');
    [ST,pred] = graphminspantree(UG);

    %% Find min traverse of spanning tree

    [minA, minB, minWts] = find(ST);
    minGraph = graph(minA, minB, minWts);
    pairDist = distances(minGraph);
    [maxDistA, maxDistB] = find(pairDist == max(pairDist(:)));

    shortSpan = shortestpath(minGraph, maxDistA(1), maxDistB(1));

    minSet = dt.Points(shortSpan, :);


    %% Fit to piecewise segments at a reasonable spacing

    % subsampleStep = 20;
    % onPoint = 5;
    % 
    % nPointsForCycle = floor((size(minSet, 1) - onPoint)/subsampleStep);
    % piecewisePoints = zeros(nPointsForCycle, 2);
    % 
    % options = optimset('display', 'off');
    % 
    % m = 1;
    % piecewisePoints(m, :) = minSet(onPoint, :);
    % startPoint = minSet(onPoint, :);
    % 
    % for k = 2:size(piecewisePoints, 1)
    % 
    %     pointsHere = minSet(onPoint:(onPoint+subsampleStep), :);
    % 
    %     testSlope = (pointsHere(end, 2) - startPoint(2))/(pointsHere(end,1) - startPoint(1));
    % 
    %     fitSlope = lsqcurvefit(@fitSpan, testSlope, ...
    %         [startPoint(1); startPoint(2); pointsHere(:,1)], pointsHere(:,2), [], [], options); 
    % 
    %     nextPoint = fitSpan(fitSlope, [startPoint(:); pointsHere(end,1)]);
    %   
    %     m = m + 1;
    %     onPoint = onPoint + subsampleStep;
    %     piecewisePoints(m, :) = [pointsHere(end, 1), nextPoint];
    %     startPoint = [pointsHere(end,1), nextPoint];
    %     
    % end

    %% Smooth points by built-in functions
    % subsampleStep = 20;


    % Sub-sample at some reasonable spacing
    subsampleStep = 20;
    onPoint = 5;

    xx = smooth(1:length(minSet(:,1)), minSet(:,1), subsampleStep, 'moving');
    yy = smooth(1:length(minSet(:,2)), minSet(:,2), subsampleStep, 'moving');

    nPointsForCycle = floor((size(minSet, 1) - onPoint)/subsampleStep);

    subsampleVector = onPoint:subsampleStep:size(minSet, 1);
    piecewisePoints = [xx(subsampleVector), yy(subsampleVector)];

    % piecewisePoints = meshfitAndSmooth2d(pts, 20, 5);

    %% Display
    figure(1)
    set(gcf, 'color', [1 1 1]);
    plot(pts(:,1), pts(:,2), 'o', 'markersize', 3, 'markeredgecolor', rgb(189, 195, 199));
    hold on
    plot(minGraph,'XData', dt.Points(:,1),'YData',dt.Points(:,2), 'edgecolor', rgb(149, 165, 166), 'marker', 'none');
    plot(minSet(:,1), minSet(:,2), '-x', 'color', rgb(155, 89, 182), 'markersize', 8);
    plot(piecewisePoints(:,1), piecewisePoints(:,2), '-x', 'linewidth', 2, 'color', rgb(231, 76, 60));
    hold off
    set(gca, 'xlim', [0, sizeX], 'ylim', [0, sizeY]);
    set(gca, 'FontSize', 12);
    xlabel('Position (nm)', 'FontSize', 12);
    ylabel('Position (nm)', 'FontSize', 12);
    patch([100, 200, 200, 100, 100], [100, 100, 120, 120, 100], 'k');
     text(150, 170, '100 nm', 'HorizontalAlignment', 'center');
    axis square
    set(gcf, 'position', [680   285   805   693])
    title(sprintf('Period = %.2f nm', calcPeriod))

    %% 
    % set(gca, 'xlim', [0, 400], 'ylim', [1850, 2250]);
    % patch([10, 110, 110, 10, 10], [1900, 1900, 1905, 1905, 1900], 'k');
    % text(55, 1915, '100 nm', 'HorizontalAlignment', 'center');

    %% Display
    figure(2)
    set(gcf, 'color', [1 1 1]);
    plot(yBase, xBase, '-', 'markersize', 3, 'color', rgb(44, 62, 80));
    hold on
    plot(pts(:,1), pts(:,2), 'o', 'markersize', 3, 'markeredgecolor', rgb(189, 195, 199));
    plot(xunit, yunit, 'o', 'markersize', 4, 'markeredgecolor', rgb(39, 174, 96));
    hold off
    set(gca, 'xlim', [0, sizeX], 'ylim', [0, sizeY]);
    set(gca, 'FontSize', 12);
    xlabel('Position (nm)', 'FontSize', 12);
    ylabel('Position (nm)', 'FontSize', 12);
    patch([100, 200, 200, 100, 100], [100, 100, 120, 120, 100], 'k');
     text(150, 170, '100 nm', 'HorizontalAlignment', 'center');
    axis square
    set(gcf, 'position', [680   285   805   693])
    title(sprintf('Period = %.2f nm', calcPeriod))

    %% Histogram
    [a, b] = hist(pts(:,end), 50);
    a = a/sum(a(:));
    figure(3);
    set(gcf, 'color', [1 1 1]);
    plot(b, a, 'linewidth', 2, 'color', rgb(231, 76, 60))
    hold on
    plot(mean(pts(:,end))*ones(2, 1), [0, max(a)*1.15], '--', 'color', rgb(52, 73, 94));
    hold off
    set(gca, 'xlim', [0, max(b)*1.1], 'ylim', [0, max(a)*1.1]);
    set(gca, 'FontSize', 12);
    xlabel('Precision (\sigma)', 'FontSize', 12);
    ylabel('PDF', 'FontSize', 12);
    set(gcf, 'position', [680   285   805   693])
     text(mean(pts(:,end))*1.1, max(a)*0.8, sprintf('\\sigma = %.1f nm',mean(pts(:,end))), 'fontsize', 12);

%% Save figures to disk
    
    print(1,fullfile(pwd, 'output', ...
                    sprintf('FitPoints2d_pts-%.f_period-%.2f.png', nPoints, calcPeriod)), '-dpng');
    print(2,fullfile(pwd, 'output', ...
                    sprintf('GroundTruth2d_pts-%.f_period-%.2f.png', nPoints, calcPeriod)), '-dpng');
    print(3,fullfile(pwd, 'output', ...
                    sprintf('Histogram_pts-%.f_period-%.2f.png', nPoints, calcPeriod)), '-dpng');
    
end

