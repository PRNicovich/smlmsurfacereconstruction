% Generate test surface for mesh fitting parameter space search

sizeX = 5000; % nm
sizeY = 5000; % nm
sizeZ = 5000; % nm

downSampleRate = 50;

amp = 1000; % amplitude in microns

k1 = 3/sizeX; % cycles over range [-1, 1], converted to cycles over [-sizeX, sizeX] range
% Would be good to get this into period, at least eventually.
k2 = 0; 

[X, Y] = meshgrid(linspace(-sizeX, sizeX, sizeX/downSampleRate));

Z = amp*cos(k1*X*pi).*cos(k2*Y*pi);

surf(X, Y, Z, 'edgecolor', 'none');

% Write to STL file
DT = delaunay(X(:), Y(:));
tr = triangulation(DT,X(:),Y(:),Z(:));

stlwrite(tr, 'parametricSurface.stl')
