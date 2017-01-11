%% Shadow Removal for Strongly Cast shadow
% CE264 - ZUO XIA, HAOYU WU - WINTER2015
close all
clear all;
Im = imread('3.JPG'); % Read image into matlab, RGB.
Threshold = 0.3; % Threshold for shadow detection
Factor = 18; % Lightness correction Factor
bluef = 0.82; % Blue correction Factor
greenf = 0.95; % Green correction Factor
k = 1; % number of times to run border removing filter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shadow_map = zeros(size(Im));
subplot(2,3,1);
imshow(Im); title('Original');
HSI = double(Im);
[x, y, c] = size(Im); % x = columns, y = rows
M = [1/3 1/3 1/3; -sqrt(6)/6 -sqrt(6)/6 -sqrt(6)/6; -1/sqrt(6) -2/sqrt(6) 0];
Ps = 0.95; % THRESHOLD
count = 0;
red = 0;
%% Shadow Detection
for i = 1:x
for j = 1:y
R = Im(i,j,1);
G = Im(i,j,2);
B = Im(i,j,3);
% Use matrix M given in research
IV1V2 = M*double([R; G; B]);
HSI(i,j,3) = IV1V2(1);
HSI(i,j,2) = sqrt(IV1V2(2)^2+IV1V2(3)^2);
HSI(i,j,1) = (atan(IV1V2(2)/IV1V2(1)+pi))*255/(2*pi);
RMh(i,j) = HSI(i,j,1)/(HSI(i,j,3)+1);
if Im(i,j,3) >= .75*Im(i,j,1) && (Im(i,j,1) >=15 && Im(i,j,2) >=15 && Im(i,j,3) >=15)
% This also removes any singularities, namely when the pixel is
% black or very close to black
if RMh(i,j) > Threshold
shadow_map(i,j) = 1;
count = count+1;
end
end
red = red + double(R);
end
end
if red/(x*y) >= 180;
% Somtimes the background is predominantly red hued, in which case the
% blue detection won't work and the shadow map will be null.
for i = 1:x
for j = 1:y
R = Im(i,j,1);
G = Im(i,j,2);
B = Im(i,j,3);
% Use matrix M given in research
IV1V2 = M*double([R; G; B]);
HSI(i,j,3) = IV1V2(1);
HSI(i,j,2) = sqrt(IV1V2(2)^2+IV1V2(3)^2);
HSI(i,j,1) = (atan(IV1V2(2)/IV1V2(1)+pi))*255/(2*pi);
RMh(i,j) = HSI(i,j,1)/(HSI(i,j,3)+1);
if (Im(i,j,1) >=15 && Im(i,j,2) >=15 && Im(i,j,3) >=15)
if RMh(i,j) > Threshold
shadow_map(i,j) = 1;
count = count+1;
end
end
end
end
end
shadow_map = bwareaopen(shadow_map,500); % Remove small sections of false shadow points
s = round(x/100);
dil = ones(s,s);
if count/(x*y) >= 0.04
shadow_map = imdilate(shadow_map,dil); % Dilate shadow region to close any holes
% But only for large images with large number of shadow pixels. Smaller
% images are left alone
shadow_map = imerode(shadow_map,dil);
% erode shadow map back to original size, now without any holes.
end
subplot(2,3,2);
imagesc(shadow_map); axis off; title('Shadow Map');
%% Shadow Removal
LAB = RGB2Lab(Im); % Convert image to CIELAB color space
sumc = 0;
counter = 0;
% Find average lightness of non-shadow pixels to create an upper limit for
% lightness correction
for i = 1:x
for j = 1:y
if shadow_map(i,j) == 0
sumc = sumc + LAB(i,j,1);
counter = counter+1;
end
end
end
avg = sumc/counter; % Average L of non-shadow pixels, used as a guide.
%% Lightness correction
for i = 1:x
for j = 1:y
if shadow_map(i,j) ~= 0
LAB(i,j,1) = LAB(i,j,1) + Factor*avg/LAB(i,j,1);
end
end
end
ImB = Lab2RGB(LAB);
ImC = ImB;
subplot(2,3,4);
imshow(ImB); title('Lightness Corrected');
%% Color correction
for i = 1:x
for j = 1:y
if shadow_map(i,j) ~=0
ImC(i,j,3) = ImC(i,j,3)*bluef; % Correct Blue
ImC(i,j,2) = ImC(i,j,2)*greenf; % Correct Green
end
end
end
subplot(2,3,5);
imshow(ImC); title('Color Corrected');
%% Border Removal
BR = zeros(size(shadow_map));
for i = 2:x-1
for j = 2:y-1
summ = shadow_map(i,j)+shadow_map(i+1,j)+shadow_map(i-1,j)+shadow_map(i,j+1)+shadow_map(i,j-1);
% If any of a pixel's surrounding pixels are different from itself,
% it is considered a border pixel. Use 4-connectivity.
if summ ~= 0 && summ ~= 5
BR(i,j) = 1;
end
end
end
BR = bwareaopen(BR,100);
ImBR = ImC;
dil2 = ones(5);
BRm = imdilate(BR,dil2); % thicken border since the affected border area on the image
% is approximately 5 pixels wide.
subplot(2,3,3); imagesc(BRm);title('Shadow Edge'); axis off;
filter = ones(3)/9;
for kk = 1:k
for i = 5:x-4
for j = 5:y-4
if BRm(i,j) == 1
% apply averaging filter,covering a 9x9 grid, excluding the
% affected pixel itself since we know it is incorrect.
RR = sum(sum(double(ImBR((i-4:i+4),(j-4:j+4),1))))/80;
GG = sum(sum(double(ImBR((i-4:i+4),(j-4:j+4),2))))/80;
BB = sum(sum(double(ImBR((i-4:i+4),(j-4:j+4),3))))/80;
ImBR(i,j,1) = uint8(round((RR)-double(ImC(i,j,1))/80));
ImBR(i,j,2) = uint8(round((GG)-double(ImC(i,j,2))/80));
ImBR(i,j,3) = uint8(round((BB)-double(ImC(i,j,3))/80));
end
end
end
end
subplot(2,3,6); imshow(ImBR);
title('Final Image');