%addpath('D:\KTH\흆 4\Image Analysis and Computer Vision (DD2423)\MATLABPATH\DD2423_Lab_Files\Functions','D:\KTH\흆 4\Image Analysis and Computer Vision (DD2423)\MATLABPATH\DD2423_Lab_Files\Images','D:\KTH\흆 4\Image Analysis and Computer Vision (DD2423)\MATLABPATH\DD2423_Lab_Files\Images-m','D:\KTH\흆 4\Image Analysis and Computer Vision (DD2423)\MATLABPATH\DD2423_Lab_Files\Images-mat')
%addpath('E:\KTH\흆 4\Image Analysis and Computer Vision (DD2423)\MATLABPATH\DD2423_Lab_Files\Functions','E:\KTH\흆 4\Image Analysis and Computer Vision (DD2423)\MATLABPATH\DD2423_Lab_Files\Images','E:\KTH\흆 4\Image Analysis and Computer Vision (DD2423)\MATLABPATH\DD2423_Lab_Files\Images-m','E:\KTH\흆 4\Image Analysis and Computer Vision (DD2423)\MATLABPATH\DD2423_Lab_Files\Images-mat')
%addpath('F:\KTH\흆 4\Image Analysis and Computer Vision (DD2423)\MATLABPATH\DD2423_Lab_Files\Functions','F:\KTH\흆 4\Image Analysis and Computer Vision (DD2423)\MATLABPATH\DD2423_Lab_Files\Images','F:\KTH\흆 4\Image Analysis and Computer Vision (DD2423)\MATLABPATH\DD2423_Lab_Files\Images-m','F:\KTH\흆 4\Image Analysis and Computer Vision (DD2423)\MATLABPATH\DD2423_Lab_Files\Images-mat')

% x increasing downwards, y increasing to the right
%% Question 1
clear all, close all, clc

tools = few256;

deltax = [0     -0.5    0
          0     0       0
          0     0.5     0];
      
deltay = deltax';

dxtools = conv2(tools, deltax, 'valid');    % Convolute the kernel to the image
dytools = conv2(tools, deltay, 'valid');

figure
subplot(1, 3, 1)
showgrey(tools)
title('original image')
axis on
xlabel('y')
ylabel('x')
subplot(1, 3, 2)
showgrey(dxtools)
title('dx')
subplot(1, 3, 3)
showgrey(dytools)
title('dy')

% We loose all the edges because the kernel cannot operate on the outer most pixels
sizeTools = size(tools);
sizeDxTools = size(dxtools);

%% Question 2-3
clear all, close all, clc

tools = few256;
blurredtools = discgaussfft(tools, 2);

dx = [0     -0.5    0
      0     0       0
      0     0.5     0];
      
dy = dx';

gradmagntools = Lv(tools, 'same');
% dxtoolsconv = conv2(tools, dx, 'valid');
% dytoolsconv = conv2(tools, dy, 'valid');
% gradmagntools = sqrt(dxtoolsconv .^2 + dytoolsconv .^2);

blurredgradmagntools = Lv(blurredtools, 'same');
% blurreddxtoolsconv = conv2(blurredtools, dx, 'valid');
% blurreddytoolsconv = conv2(blurredtools, dy, 'valid');
% blurredgradmagntools = sqrt(blurreddxtoolsconv .^2 + blurreddytoolsconv .^2);

figure
subplot(2, 4, 1)
showgrey(tools)
title('original')
axis on
xlabel('y')
ylabel('x')
subplot(2, 4, 2)
showgrey(gradmagntools)
title('gradmagntools')
subplot(2, 4, 3)
showgrey(blurredgradmagntools)
title('smoothed gradmagntools')
subplot(2, 4, 4)
histogram(gradmagntools)
title('smoothed gradmagntools')
axis([0 40 0 16000])
xlabel('magnitude')

threshold_tools = [5 10 15 20];

i = 0;
for i = 1:length(threshold_tools)
    subplot(2,4,4+i)
    showgrey((blurredgradmagntools - threshold_tools(i)) > 0)
    title(['threshold = ', num2str(threshold_tools(i))])
end


house = godthem256;
blurredhouse = discgaussfft(house, 2);

gradmagnhouse = Lv(house, 'same');
% dxhouseconv = conv2(house, dx, 'valid');
% dyhouseconv = conv2(house, dy, 'valid');
% gradmagnhouse = sqrt(dxhouseconv .^2 + dyhouseconv .^2);

blurredgradmagnhouse = Lv(blurredhouse, 'same');
% blurredxhouseconv = conv2(blurredhouse, dx, 'valid');
% blurredyhouseconv = conv2(blurredhouse, dy, 'valid');
% blurredgradmagnhouse = sqrt(blurredxhouseconv .^2 + blurredyhouseconv .^2);

figure
subplot(2, 4, 1)
showgrey(house)
title('original')
axis on
xlabel('y')
ylabel('x')
subplot(2, 4, 2)
showgrey(gradmagnhouse)
title('gradmagnhouse')
subplot(2, 4, 3)
showgrey(blurredgradmagnhouse)
title('smoothed gradmagnhouse')
subplot(2, 4, 4)
histogram(gradmagnhouse)
title('smoothed gradmagnhouse')
axis([0 60 0 16000])
xlabel('magnitude')

threshold_house = [5 10 15 20];

i = 0;
for i = 1:length(threshold_house)
    subplot(2,4,4+i)
    showgrey((blurredgradmagnhouse - threshold_house(i)) > 0)
    title(['threshold = ', num2str(threshold_house(i))])
end


smoothing = [0 1 2 3 6 10];
threshold = 10;

figure
for i = 1:length(smoothing)
    subplot(2,3,i)
    smoothedhouse = discgaussfft(house, smoothing(i));
    smoothedgradmagnhouse = Lv(smoothedhouse, 'same');
    showgrey((smoothedgradmagnhouse - threshold) > 0)
    title(['Smoothing = ', num2str(smoothing(i))])
end

%% Question 4-6
clear all, close all, clc

house = godthem256;
tools = few256;

scale = [0.0001, 1.0, 4.0, 16.0, 64.0];

figure
i = 0;
for i = 1:length(scale)
    subplot(3,2,i)
    contour(Lvvtilde(discgaussfft(house, scale(i)), 'same'), [0 0]);
    axis('image')
    axis('ij')
    title(['scale = ', num2str(scale(i))])
end

figure
i = 0;
for i = 1:length(scale)
    subplot(3,2,i)
    showgrey(Lvvvtilde(discgaussfft(tools, scale(i)), 'same') < 0);
    axis('image')
    axis('ij')
    title(['scale = ', num2str(scale(i))])
end

i = 0;
figure
for i = 1:length(scale)
    subplot(3,2,i)
    showgrey(Lvvvtilde(discgaussfft(tools, scale(i)), 'same') < -10);
    axis('image')
    axis('ij')
    title(['scale = ', num2str(scale(i))])
end

%% Derivative check
dy = [0     0       0       0       0;
    0     0       -0.5    0       0;
    0     0       0       0       0;
    0     0       0.5     0       0;
    0     0       0       0       0];
dx = dy';

dyy = [0     0       0       0       0;
    0     0       1       0       0;
    0     0       -2      0       0;
    0     0       1       0       0;
    0     0       0       0       0];
dxx = dyy';

dxy = conv2(dx, dy, 'same');
dxxx = conv2(dx, dxx, 'same');
dyyy = conv2(dy, dyy, 'same');
dxxy = conv2(dx, dxy, 'same');

[x, y] = meshgrid(-5:5, -5:5);

checkxxx = filter2(dxxx, x .^3, 'valid');
checkxx = filter2(dxx, x .^3, 'valid');
checkxxy = filter2(dxxy, x .^2 .* y, 'valid');

%% Question 7
clear all, close all, clc

tools = few256;
house = godthem256;
threshold = [1 3 5 7 9 11 13 15];
scale = 4;

figure
rows = 3;
columns = 3;
subplot(rows, columns, 1)
showgrey(house)
title('Original')
axis on
xlabel('y')
ylabel('x')

for i = 1:length(threshold)
    edgecurves = extractedge(house, scale, threshold(i), 'same');
    
    subplot(rows, columns, i+1)
    overlaycurves(house, edgecurves)
    title(['threshold = ', num2str(threshold(i))])
end
% sgtitle(['scale = ', num2str(scale)])

%% Question 8
clear all
% close all
clc

testimage1 = triangle128;
smalltest1 = binsubsample(testimage1);

testimage2 = houghtest256;
smalltest2 = binsubsample(binsubsample(testimage2())); 

tools = few256;
phone = phonecalc256;
house = godthem256;

pic = house;
scale = 3;
gradmagnthreshold = 10;
nrho = 800;
ntheta = 100;
nlines = 15;
verbose = 1;

[linepar, acc] = houghedgeline(pic, scale, gradmagnthreshold, nrho, ntheta, nlines, verbose);

cellsChange = nrho * ntheta;
cells = [10000, 21000, 30000, 43750, 60000];
time = [0.79, 0.89, 1.01, 1.21, 1.35];

% subplot(1,2,2)
% showgrey(acc)
% title('Hough')

% figure
% plot(time, cells)
% xlabel('Computational time')
% ylabel('Number of cells in the accumulator')
