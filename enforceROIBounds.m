% Remove all points that fall outsize of [sizeX, sizeY, sizeZ] bounds

function pts = enforceROIBounds(pts, sizeX, sizeY, sizeZ)

pts(pts(:,1) > sizeX | pts(:,2) > sizeY | pts(:,3) > sizeZ, :) =[ ];
pts(pts(:,1) < 0 | pts(:,2) < 0 | pts(:,3) < 0, :) =[ ];