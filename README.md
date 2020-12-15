## SMLM Surface Reconstruction

Reconstruct a surface mesh from single molecule localization microscopy (SMLM) data.  Take an arbitrary structure (as STL file), generate a simulated 3D SMLM dataset around it, then denoise and reconstruct a surface.  

![smlmReconBunny](/output/smlm3dreconBunny.gif)

Code does following steps:
*  Take in reference structure as STL file
*  Read example SMLM data file (nb - Nikon, Zeiss, and some versions of ThunderSTORM output supported)
*  Simulate 3D SMLM dataset based on STL file and precision information of example SMLM file
*  Denoise data by model-free fitting of regularly-spaced meridians along orthogonal axes
*  Use a noise-tolerant surface reconstruction to generate surface
*  Analyze surface features

./manuscript/smlm3dmeshfitting.pdf has more details on execution of code and math behind it.

./MeshFitting3d.m is demo function.  Run to generate example output.

## Requires:

MeshLab
https://www.meshlab.net/

ICP_finite.m
https://www.mathworks.com/matlabcentral/fileexchange/24301-finite-iterative-closest-point

Tested on MATLAB 2019a on PC (Intel i5; 8 GB RAM)

## Distribution:
Code is open for those who wish to explore this approach for 3D sMLM reconstructions.  It is currently not published elsewhere.  If you intend to publish, please contact regarding attribution. Released with GNU license.

