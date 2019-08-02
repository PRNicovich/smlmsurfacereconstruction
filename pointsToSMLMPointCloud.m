% Generate proxy SMLM point cloud given base point cloud and example data. 
% Example data usually pulled from experimental dataset
% Base point cloud represents true molecular positions and returned point
% cloud is what would be detected in experiment. 

function pts = pointsToSMLMPointCloud(xunit, yunit, zunit, ...
                                        nPoints, nDetectsPerPoint,photonsAndPrecision)

nPtsActual = repmat(nDetectsPerPoint, [nPoints, 1]);
nPtsActual = poissrnd(nPtsActual);

pts = zeros(sum(nPtsActual(:)), 5); %[x, y, z, nPhotons, locPrecision]
ptsNow = 1;
for k = 1:size(nPtsActual, 1)
    
    rSamp = randsample(numel(photonsAndPrecision(:,1)), nPtsActual(k), 'true');
    precHere = photonsAndPrecision(rSamp, 2).^2;
    
    pts(ptsNow:(ptsNow + nPtsActual(k)-1),:) = ...
        [xunit(k) + precHere.*randn(nPtsActual(k), 1), ...
         yunit(k) + precHere.*randn(nPtsActual(k), 1), ...
         zunit(k) + precHere.*randn(nPtsActual(k), 1), ...
         photonsAndPrecision(rSamp, 2), ...
         precHere] ;
     
     ptsNow = ptsNow + nPtsActual(k);
    
end
