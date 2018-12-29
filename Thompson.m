function sigma2_xy = Thompson(s, a, N, b)

% Modified Thompson equation for estimating fluorophore localization precision
% From:
% K.I. Mortensen, L.S. Churchman, J.A. Spudich, and Flyvbjerg, 
% “Optimized localization analysis for singlemolecule tracking and 
% super-resolution microscopy,” Nature Methods, vol. 7, no. 5, 
% pp. 377–384, 2010.
% As presented in:
% Bernd Rieger and Sjoerd Stallinga
% THE EFFECT OF BACKGROUND ON LOCALIZATION UNCERTAINTY IN 
% SINGLE EMITTER IMAGING
% ISB IEEE (2012) , p. 988-991.
% Actual Thompson equation appears as sig_u value below.
% From Thompson, Larson, and Webb, Biophysical (82) 2002 2775-2783.

%
% Inputs:
% s - standard deviation of PSF (nm)
% a - pixel dimension (nm)
% N - number of fluorophore photons
% b - background variance
% Output - sigma_xy (2D uncertainty) in fluorophore localization, units nm^2.  In Zeiss
% 1.txt file, this value is reported as sqrt(sigma_xy) in nm.  

sig_a = s.^2 + ((a.^2)/12); % sig_a^2 in 

sigma2_xy = ((sig_a)./N).*(16/9 + ((8*pi*(sig_a).*(b))./(N.*(a.^2))));

% MLE approximation from IEEE paper
% tau = (2*pi*(sig_a).*b)./(N.*(a.^2));
% 
% sig_MLE = ((sig_a)./N).*(1 + 4.*tau + sqrt((2*tau)./(1 + 4*tau)));

% From Zeiss website:  
% sig_u = ((s.^2)./N) + (((a.^2)/12)./N) + ((8*pi*(s.^4).*(b))./((a.^2).*(N.^2)));
% sigma2_xy = sig_u;
% None match result from 1.txt file.  All underreport precision.