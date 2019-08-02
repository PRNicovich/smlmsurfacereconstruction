% strOut = '"C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\AdaptiveSolvers.x64\PoissonRecon.exe" --in "C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\ptCloud20190527.ply" --out out.ply --depth 10 --samplesPerNode 15 --bType 2';

strOut = '"C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\AdaptiveSolvers.x64\SSDRecon.exe" --in "C:\Users\Rusty Nicovich\Documents\MATLAB\meshFitting\ptCloud20190527.ply" --out out.ply --depth 10 --samplesPerNode 15';


[~, ~] = dos(strOut, '-echo');