%%
%Upg 1.1
load('Dell.mat');
load('Inkjet.mat')
plot_chrom(XYZdell, 'b');
hold on 
plot_chrom(XYZinkjet, 'g');

%%
%Upg 2.1.1
img = imread('peppers_gray.tif');

imgNe = imresize(imresize(img,0.25,'nearest'),4,'nearest');
imgLi = imresize(imresize(img,0.25,'bilinear'),4,'bilinear');
imgBi = imresize(imresize(img,0.25,'bicubic'),4,'bicubic');

SnrImgNe = mysnr(img, img-imgNe);
MseImgNe = immse(img, imgNe);

SnrImgLi = mysnr(img, img-imgLi);
MseImgLi = immse(img, imgLi);

SnrImgBi = mysnr(img, img-imgBi);
MseImgBi = immse(img, imgBi);

imshow(img);
figure
imshow(imgNe);
figure
imshow(imgLi);
figure
imshow(imgBi);

%%
%upg 2.1.2
imgDouble = im2double(img);
imgThres = im2double(imgDouble >= 0.5);
imgDith = im2double(dither(imgDouble));imgBi

SnrImgThres = mysnr(imgDouble, imgDouble-imgThres); 
MseImgThres = immse(imgDouble, imgThres);

SnrImgDith = mysnr(imgDouble, imgDouble-imgDith);
MseImgDith = immse(imgDouble, imgDith);

imshow(imgDouble);
figure
imshow(imgThres);
figure
imshow(imgDith);



%%
%upg 2.2
imrgb = im2double(imread('peppers_color.tif'));

imLab = rgb2lab(imrgb);

imTresh(:,:,1) = (imrgb(:,:,1) >= 0.5);
imTresh(:,:,2) = (imrgb(:,:,2) >= 0.5);
imTresh(:,:,3) = (imrgb(:,:,3) >= 0.5);

imDith(:,:,1) = dither(imrgb(:,:,1));
imDith(:,:,2) = dither(imrgb(:,:,2));
imDith(:,:,3) = dither(imrgb(:,:,3));


treshRGB = double(imTresh);
ditchRGB = double(imDith);

treshLab = rgb2lab(treshRGB);
ditchLab = rgb2lab(ditchRGB);

deltaETresh= 0;

for i = 1:1:512
    for j = 1:1:512
    
    deltaETresh = deltaETresh + sqrt((treshLab(i,j,1)-imLab(i,j,1))^2 +  (treshLab(i,j,2)-imLab(i,j,2))^2+ (treshLab(i,j,3)-imLab(i,j,3))^2);
    
    end

end

deltaETresh = deltaETresh/(512*512);

deltaEDith= 0;

for i = 1:1:512
    for j = 1:1:512
    
    deltaEDith = deltaEDith + sqrt((ditchLab(i,j,1)-imLab(i,j,1))^2 +  (ditchLab(i,j,2)-imLab(i,j,2))^2 + (ditchLab(i,j,3)-imLab(i,j,3))^2);
    
    end

end

deltaEDith = deltaEDith/(512*512);

imshow(imrgb);
figure
imshow(treshRGB);
figure
imshow(ditchRGB);

%%
%upg 3
img = imread('peppers_gray.tif');
imgDouble = im2double(img);
imgThres = im2double(imgDouble >= 0.5);
imgDith = im2double(dither(imgDouble));


snrDith = snr_filter(imgDouble, imgDouble-imgDith);
snrTresh = snr_filter(imgDouble, imgDouble-imgThres);


imshow(imgDouble)
figure
imshow(imgThres)
figure
imshow(imgDith)


%%
%upg 3.2

imrgb = im2double(imread('peppers_color.tif'));

imLab = rgb2lab(imrgb);

imTresh(:,:,1) = (imrgb(:,:,1) >= 0.5);
imTresh(:,:,2) = (imrgb(:,:,2) >= 0.5);
imTresh(:,:,3) = (imrgb(:,:,3) >= 0.5);

imDith(:,:,1) = dither(imrgb(:,:,1));
imDith(:,:,2) = dither(imrgb(:,:,2));
imDith(:,:,3) = dither(imrgb(:,:,3));

treshRGB = double(imTresh);
ditchRGB = double(imDith);


f = MFTsp(15,0.0847,500);

%ref
ref(:,:,1) = conv2(imrgb(:,:,1),f,'same');
ref(:,:,2) = conv2(imrgb(:,:,2),f,'same');
ref(:,:,3) = conv2(imrgb(:,:,3),f,'same');

%Thresh
Thresh(:,:,1) = conv2(treshRGB(:,:,1),f,'same');
 Thresh(:,:,1)=(Thresh(:,:,1)>0).*Thresh(:,:,1);
Thresh(:,:,2) = conv2(treshRGB(:,:,2),f,'same');
 Thresh(:,:,2)=(Thresh(:,:,2)>0).*Thresh(:,:,2);
Thresh(:,:,3) = conv2(treshRGB(:,:,3),f,'same');
 Thresh(:,:,3)=(Thresh(:,:,3)>0).*Thresh(:,:,3);

%Thresh
Dith(:,:,1) = conv2(imDith(:,:,1),f,'same');
 Dith(:,:,1)=(Dith(:,:,1)>0).*Dith(:,:,1);
Dith(:,:,2) = conv2(imDith(:,:,2),f,'same');
 Dith(:,:,2)=(Dith(:,:,2)>0).*Dith(:,:,2);
Dith(:,:,3) = conv2(imDith(:,:,3),f,'same');
 Dith(:,:,3)=(Dith(:,:,3)>0).*Dith(:,:,3);
refLab = rgb2lab(ref);
treshLab = rgb2lab(Thresh);
ditchLab = rgb2lab(Dith);

deltaETresh= 0;

for i = 1:1:512
    for j = 1:1:512
    
    deltaETresh = deltaETresh + sqrt((treshLab(i,j,1)-refLab(i,j,1))^2 +  (treshLab(i,j,2)-refLab(i,j,2))^2+ (treshLab(i,j,3)-refLab(i,j,3))^2);
    
    end

end

deltaETresh = deltaETresh/(512*512);

deltaEDith= 0;

for i = 1:1:512
    for j = 1:1:512
    
    deltaEDith = deltaEDith + sqrt((ditchLab(i,j,1)-refLab(i,j,1))^2 +  (ditchLab(i,j,2)-refLab(i,j,2))^2 + (ditchLab(i,j,3)-refLab(i,j,3))^2);
    
    end

end

deltaEDith = deltaEDith/(512*512);



imshow(ref)
figure
imshow(Thresh)
figure
imshow(Dith)


%%
%upg 4.1

img = imread('peppers_color.tif');


imgNe = imresize(imresize(img,0.25,'nearest'),4,'nearest');
imgLi = imresize(imresize(img,0.25,'bilinear'),4,'bilinear');
imgBi = imresize(imresize(img,0.25,'bicubic'),4,'bicubic');

imgXYZ = rgb2xyz(img);
nearestXYZ = rgb2xyz(imgNe);
linearXYZ = rgb2xyz(imgLi);
bicubicXYZ = rgb2xyz(imgBi);

sampDegree = visualAngle(-1,200/2.5,300,1);
whitePoint = [95.05,100,108.9];
near = scielab(sampDegree,imgXYZ,nearestXYZ,whitePoint,'xyz');
linear = scielab(sampDegree,imgXYZ,linearXYZ,whitePoint,'xyz');
bicubic = scielab(sampDegree,imgXYZ,bicubicXYZ,whitePoint,'xyz');

nearMean = mean(mean(near));
linearMean = mean(mean(linear));
bicubicMean = mean(mean(bicubic));

imshow(img);
figure
imshow(imgNe);
figure
imshow(imgLi);
figure
imshow(imgBi);




%%
%upg 4.2
load('colorhalftones.mat');

sampDegree = visualAngle(-1,500/2.5,300,1);

c1xyz = rgb2xyz(c1);
c2xyz = rgb2xyz(c2);


c1Lab = scielab(sampDegree,c1xyz);
c2Lab = scielab(sampDegree,c2xyz);

c1Sum = sum(std(c1Lab(:,:,1))) + sum(std(c1Lab(:,:,2)))+sum(std(c1Lab(:,:,3)));
c2Sum = sum(std(c2Lab(:,:,1)))+sum(std(c2Lab(:,:,2)))+sum(std(c2Lab(:,:,3)));


imshow(c1);
figure;
imshow(c2);

%%
%upg 4.2.2

load('colorhalftones.mat');

sampDegree = visualAngle(-1,500/2.5,300,1);

c3xyz = rgb2xyz(c3);
c4xyz = rgb2xyz(c4);
c5xyz = rgb2xyz(c5);


c3Lab = scielab(sampDegree,c3xyz);
c4Lab = scielab(sampDegree,c4xyz);
c5Lab = scielab(sampDegree,c5xyz);

c3Sum = sum(std(c3Lab(:,:,1))) + sum(std(c3Lab(:,:,2)))+sum(std(c3Lab(:,:,3)));
c4Sum = sum(std(c4Lab(:,:,1)))+sum(std(c4Lab(:,:,2)))+sum(std(c4Lab(:,:,3)));
c5Sum = sum(std(c5Lab(:,:,1)))+sum(std(c5Lab(:,:,2)))+sum(std(c5Lab(:,:,3)));


imshow(c3);
figure;
imshow(c4);
figure;
imshow(c5);




%%
%upg5.1
img = imread('peppers_gray.tif');
 img = im2double(img);
imgA1 = img;
imgA2 = img;

for i= 1:1:512
   for j = 1:1:512
    if( mod(i,2) == 0)
        imgA1(i,j) = img(i,j) - 0.1;
    else
        imgA1(i,j) = img(i,j) + 0.1;
    end
   end
end

for i= 1:1:512
   for j = 1:1:512
    if( i< 256)
        imgA2(i,j) = img(i,j) + 0.1;
    else
        imgA2(i,j) = img(i,j) - 0.1;
    end
   end
end

SnrA1 = mysnr(img, img-imgA1); 
SnrA2 = mysnr(img, img-imgA2); 


[A1, SSIMA1] = ssim(img, imgA1);
[A2,SSIMA2] = ssim(img, imgA2);
imshow(SSIMA1);
figure
imshow(SSIMA2);

%%
%upg 5.2
img = imread('peppers_gray.tif');
img = im2double(img);

dist1=img+0.2*(rand(size(img))-0.5);

f=fspecial('gauss',21,10);
dist2 = imfilter(img,f);


SnrB1 = mysnr(img, img-dist1); 
SnrB2 = mysnr(img, img-dist2); 


[B1, SSIMB1] = ssim(img, dist1);
[B2,SSIMB2] = ssim(img, dist2);

imshow(img);
figure
imshow(SSIMB1);
figure
imshow(SSIMB2);


