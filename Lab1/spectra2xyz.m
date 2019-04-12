function XYZ=spectra2xyz(reflectance,illumination)
% XYZ=spectral2xyz(reflectance,illumination)
%
% reflectance = Reflektansspektra
% illumination = Belysningens spektralf�rdelning.
% XYZ = CIEXYZ-v�rden.
%
% Reflektansspektra multipliceras med belysningens spektralf�rdelning
% och �vers�tts sedan med hj�lp av CIE:s standardobservat�r till
% CIEXYZ-v�rden. F�rgmatchningsfunktionerna �terfinns i filen
% "spectra.mat". Gl�m inte att ber�kna k-v�rdet f�r belysningen och skala
% CIEXYZ-v�rdena med erh�llet k-v�rde! (Formler f�r detta �terfinns i
% laborationsh�ftet, %Ekvation 1.4 och 1.5)
% Om CIEXYZ-v�rden ska ber�knas f�r en belysning (dvs en belysnings
% vitpunkt) s� kan reflektanspektrat ers�ttas med ettor (full reflektans i
% alla v�gl�ngder). Ex: reflectance = ones(81,1);
% Anv�nd funktionen load f�r att h�mta f�rgmatchningsfunktioner.
% Din kod h�r:

load('spectra.mat');
load('xyz.mat');

E = reflectance.*illumination;

K = 100/sum(illumination.*xyz(:,2));
XYZ(1,1) = K * sum(E.*xyz(:, 1)); 
XYZ(1,2) =K * sum(E.*xyz(:, 2));  
XYZ(1,3) = K * sum(E.*xyz(:, 3));  

return;