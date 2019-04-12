
%%
%upg1.1
load('Ad.mat');
load('Ad2.mat');

plot(Ad);
figure
plot(Ad2);
%%
%upg1.2
load('chips20.mat');
load('illum.mat');
d1 = Ad'*(chips20(:,:).*CIED65)';
d2 = Ad2'*(chips20(:,:).*CIED65)';

showRGB(d1');

showRGB(d2');

%%
%upg 2.1
e = ones(1,61);
norm1 = Ad'*e';
norm2 = Ad2'*e';


plot(norm1);
figure
plot(norm2);
%%
%upg 2.2

RGB_cal_D65_Ad1 = (Ad'*(chips20(:,:).*CIED65)')./norm1;
RGB_cal_D65_Ad2 = (Ad2'*(chips20(:,:).*CIED65)')./norm2;


showRGB(RGB_cal_D65_Ad1');
showRGB(RGB_cal_D65_Ad2');

%%
%upg 2.3


plot(CIED65);
figure
plot(CIEA);

%%
%upg 2.4

RGB_A = (Ad'*(chips20(:,:).*CIEA)')./norm1;
RGB_D65 = (Ad'*(chips20(:,:).*CIED65)')./norm1;


showRGB(RGB_A');
showRGB(RGB_D65');


%%
%upg 2.5
R= ones(1,61);
cal_D65 = (Ad'*(R.*CIED65)');
cal_A = (Ad'*(R.*CIEA)');


RGB_cal_A = (Ad'*(chips20(:,:).*CIEA)')./cal_A;
RGB_cal_D65 = (Ad'*(chips20(:,:).*CIED65)')./cal_D65;


showRGB(RGB_cal_A');
showRGB(RGB_cal_D65');


%%
%upg 3.1
clc
for i = 1:1:20
    XYZ_D65_ref(i,:) = spectra2xyz(chips20(i,:)',CIED65');
end



%%
%upg 3.2
clc


MXYZ = inv(M_XYZ2RGB) * RGB_cal_D65_Ad1;

%Matrix
for i = 1:1:20
    [L, a, b] = xyz2lab(MXYZ(1,i), MXYZ(2,i), MXYZ(3,i));
    MLAB(i,:) = [L, a, b];
end

%reference
for i = 1:1:20
    [L, a, b] = xyz2lab(XYZ_D65_ref(i,1), XYZ_D65_ref(i,2), XYZ_D65_ref(i,3));
    MLAB_ref(i,:) = [L, a, b];
end

%Difference in color
for i = 1:1:20
    DeltaE(i) = sqrt((MLAB(i,1) - MLAB_ref(i,1))^2 + (MLAB(i,2) - MLAB_ref(i,2))^2 + (MLAB(i,3) - MLAB_ref(i,3))^2); 
end
    
DeltaEMean = mean(DeltaE);
DeltaEMax = max(DeltaE);
    
    
%%
%upg 3.3

plot(Ad);
figure
plot(xyz)



%%
% upg 3.4


A = pinv(RGB_cal_D65_Ad1')*XYZ_D65_ref;
    
conversionMatrix =  RGB_cal_D65_Ad1'*A;


%Matrix
for i = 1:1:20
    [L, a, b] = xyz2lab(conversionMatrix(i,1), conversionMatrix(i,2), conversionMatrix(i,3));
    ALAB(i,:) = [L, a, b];
end



%Difference in color
for i = 1:1:20
    Adelta(i) = sqrt((ALAB(i,1) - MLAB_ref(i,1))^2 + (ALAB(i,2) - MLAB_ref(i,2))^2 + (ALAB(i,3) - MLAB_ref(i,3))^2); 
end


Amean = mean(Adelta);
Amax = max(Adelta);
    
    
%%    
%upg 3.5
clc
AOpt = Optimize_poly(RGB_cal_D65_Ad1,XYZ_D65_ref');

estOpt = Polynomial_regression(RGB_cal_D65_Ad1,AOpt);


%Matrix
for i = 1:1:20
    [L, a, b] = xyz2lab(estOpt(1,i), estOpt(2,i), estOpt(3,i));
    OptLAB(i,:) = [L, a, b];
end



%Difference in color
for i = 1:1:20
    OptDelta(i) = sqrt((OptLAB(i,1) - MLAB_ref(i,1))^2 + (OptLAB(i,2) - MLAB_ref(i,2))^2 + (OptLAB(i,3) - MLAB_ref(i,3))^2); 
end


Optmean = mean(OptDelta);
Optmax = max(OptDelta);

%%
    
    