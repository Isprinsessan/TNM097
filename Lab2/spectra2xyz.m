function XYZ=spectra2xyz(reflectance,illumination)
% XYZ=spectral2xyz(reflectance,illumination)
%
% reflectance = Reflektansspektra
% illumination = Belysningens spektralfördelning.
% XYZ = CIEXYZ-värden.
%
% Reflektansspektra multipliceras med belysningens spektralfördelning
% och översätts sedan med hjälp av CIE:s standardobservatör till
% CIEXYZ-värden. Färgmatchningsfunktionerna återfinns i filen
% "spectra.mat". Glöm inte att beräkna k-värdet för belysningen och skala
% CIEXYZ-värdena med erhållet k-värde! (Formler för detta Återfinns i
% laborationshäftet, %Ekvation 1.4 och 1.5)
% Om CIEXYZ-värden ska beräknas för en belysning (dvs en belysnings
% vitpunkt) så kan reflektanspektrat ersättas med ettor (full reflektans i
% alla våglängder). Ex: reflectance = ones(81,1);
% Använd funktionen load för att hämta färgmatchningsfunktioner.
% Din kod här:

load('spectra.mat');
load('xyz.mat');

E = reflectance.*illumination;

K = 100/sum(illumination.*xyz(:,2));
XYZ(1,1) = K * sum(E.*xyz(:, 1)); 
XYZ(1,2) =K * sum(E.*xyz(:, 2));  
XYZ(1,3) = K * sum(E.*xyz(:, 3));  

return;