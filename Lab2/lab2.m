%%
% upg 1.1
load('TRC_display.mat');
plot(TRCb, 'b');
hold on
plot(TRCg,'g');
hold on
plot(TRCr, 'r');
linear = 0:0.01:1;
hold on
plot(linear,'y');

%%
%upg 1.2.1
load('ramp_display.mat');
load('ramp_linear.mat');
imshow(Ramp_display);
figure
imshow(Ramp_linear);

%%
%upg 1.2.2
TRCMatrix(:,:,1) = TRCr;
TRCMatrix(:,:,2) = TRCg;
TRCMatrix(:,:,3) = TRCb;

linearized_img = Linearize(TRCMatrix,Ramp_display);

imshow(linearized_img);
figure
imshow(Ramp_linear);

%%
%upg 1.3

GammaR = 2.1;
GammaG = 2.4;
GammaB = 1.8;

Dmax = 1;

D(:,:,1) = Dmax.*(Ramp_display(:,:,1).^(1/GammaR));
D(:,:,2) = Dmax.*(Ramp_display(:,:,2).^(1/GammaG));
D(:,:,3) = Dmax.*(Ramp_display(:,:,3).^(1/GammaB));

imshow(D);

%%
%upg 2.1
load('DLP.mat');

plot(DLP);

%% 
%upg 2.2
load('RGB_raw');
load('illum.mat');

for i = 1:1:20
    SRGB(:,i) = RGB_raw(1, i)*DLP(:,1) + RGB_raw(2, i)*DLP(:,2) + RGB_raw(3, i)*DLP(:,3);
end

for i = 1:1:20
    SXYZ(i,:) = spectra2xyz(SRGB(:,i),CIED65');
end

[DeltaMean, DeltaMax] = DeltaE(XYZ_ref', SXYZ');

%%
%upg 2.3
load('RGB_cal');
load('illum.mat');

for i = 1:1:20
    SRGB(:,i) = RGB_cal(1, i)*DLP(:,1) + RGB_cal(2, i)*DLP(:,2) + RGB_cal(3, i)*DLP(:,3);
end

for i = 1:1:20
    SXYZ(i,:) = spectra2xyz(SRGB(:,i),CIED65');
end

[DeltaMean, DeltaMax] = DeltaE(XYZ_ref', SXYZ');

%%
%upg 3.1
load('xyz.mat');

k = 100/sum(CIED65*xyz(:,2));

Acrt = k*xyz'*DLP;


%%
%upg 3.2
load('Ad.mat');
%3.4

e = ones(1,61);
norm1 = Ad'*e';

for i = 1:1:20
    XYZ_D65_ref(i,:) = spectra2xyz(chips20(i,:)',CIED65');
end


RGB_cal_D65_Ad1 = (Ad'*(chips20(:,:).*CIED65)')./norm1;

A = pinv(RGB_cal_D65_Ad1')*XYZ_D65_ref;
    
conversionMatrix =  RGB_cal_D65_Ad1'*A;

%3.5
AOpt = Optimize_poly(RGB_cal_D65_Ad1,XYZ_D65_ref');

estOpt = Polynomial_regression(RGB_cal_D65_Ad1,AOpt);

%

DPrim = Acrt\estOpt;

S = DLP * DPrim;

k = 100/sum(CIED65*xyz(:,2));
XYZ =  k*S'*xyz;

[DeltaMean, DeltaMax] = DeltaE(XYZ_ref', XYZ');
%%
%upg 3.4


DPrimRange = min(max(DPrim,0),1);

S = DLP * DPrimRange;

k = 100/sum(CIED65*xyz(:,2));
XYZ =  k*S'*xyz;

[DeltaMean, DeltaMax] = DeltaE(XYZ_ref', XYZ');
%%
%3.5

plot_chrom_sRGB(Acrt);

%%
%3.6

plot(S(:,1),'y');
hold on
ch = (chips20(1,:).*CIED65);
plot(ch, 'b');



