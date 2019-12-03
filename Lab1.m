%% Question 1
clear all, close all, clc

figure
p = 5;
q = 9;
sz = 128;
fftwavetemplate(p, q, sz)

figure
p = 9;
q = 5;
fftwavetemplate(p, q, sz)

figure
p = 17;
q = 9;
fftwavetemplate(p, q, sz)

figure
p = 17;
q = 121;
fftwavetemplate(p, q, sz)

figure
p = 5;
q = 1;
fftwavetemplate(p, q, sz)

figure
p = 125;
q = 1;
fftwavetemplate(p, q, sz)

%% Question 2
clear all, close all, clc

figure
p = 1;
q = 1;
sz = 128;
fftwavetemplate(p, q, sz)

figure
p = 3;
q = 1;
sz = 128;
fftwavetemplate(p, q, sz)

figure
p = 10;
q = 1;
sz = 128;
fftwavetemplate(p, q, sz)

%% Question 3
clear all, close all, clc

u = 5;
v = 9;
sz = 128;

Fhat = zeros(sz);           % zero matrix with sz x sz
Fhat(u, v) = 1;             % set (u, v)-coordinate to white
F = ifft2(Fhat);            % inverse discrete Fourier transform
Fabsmax = max(abs(F(:)));   % max-value of the absolute of F


subplot(3, 2, 1);
showgrey(Fhat);
title(sprintf('Fhat: (u, v) = (%d, %d)', u, v))
% What is done by these instructions?   % centering
if (u <= sz/2)  % if point is in the upper half
    uc = u - 1;
else            % if point is in bottom part
    uc = u - 1 - sz;
end
if (v <= sz/2)  % if point is in left half
    vc = v - 1;
else            % if point is in right half
    vc = v - 1 - sz;
end

wvl = sz/(sqrt(uc^2 + vc^2));       % wavelength

wavelength = wvl;       % Replace by correct expression
amplitude = 1/sz^2;     % Replace by correct expression

subplot(3, 2, 2);
showgrey(fftshift(Fhat));
title(sprintf('centered Fhat: (uc, vc) = (%d, %d)', uc, vc))

subplot(3, 2, 3);
showgrey(real(F), 64, -Fabsmax, Fabsmax);
title('real(F)')

subplot(3, 2, 4);
showgrey(imag(F), 64, -Fabsmax, Fabsmax);
title('imag(F)')

subplot(3, 2, 5);
showgrey(abs(F), 64, -Fabsmax, Fabsmax);
title(sprintf('abs(F) (amplitude %f)', amplitude))

subplot(3, 2, 6);
showgrey(angle(F), 64, -pi, pi);
title(sprintf('angle(F) (wavelength %f)', wavelength))

%% Question 5
clear all, close all, clc

uVek = [5 60 80];
j = 1;

for i = 1:3
    u = uVek(i);
    v = 9;
    sz = 128;
    
    Fhat = zeros(sz);
    Fhat(u, v) = 1;
    
    subplot(3, 2, j);
    j = j + 1;
    showgrey(Fhat);
    title(sprintf('Fhat: (u, v) = (%d, %d)', u, v))
    
    if (u <= sz/2)
        uc = u - 1;
    else
        uc = u - 1 - sz;
    end
    if (v <= sz/2)
        vc = v - 1;
    else
        vc = v - 1 - sz;
    end
    
    subplot(3, 2, j);
    j = j + 1;
    showgrey(fftshift(Fhat));
    title(sprintf('centered Fhat: (uc, vc) = (%d, %d)', uc, vc))
end


vVek = [9 60 80];
j = 1;
figure

for i = 1:3
    u = 5;
    v = vVek(i);
    sz = 128;
    
    Fhat = zeros(sz);
    Fhat(u, v) = 1;
    
    subplot(3, 2, j);
    j = j + 1;
    showgrey(Fhat);
    title(sprintf('Fhat: (u, v) = (%d, %d)', u, v))
    
    if (u <= sz/2)
        uc = u - 1;
    else
        uc = u - 1 - sz;
    end
    if (v <= sz/2)
        vc = v - 1;
    else
        vc = v - 1 - sz;
    end
    
    subplot(3, 2, j);
    j = j + 1;
    showgrey(fftshift(Fhat));
    title(sprintf('centered Fhat: (uc, vc) = (%d, %d)', uc, vc))
end

%% Question 7-9
clear all, close all, clc

F = [ zeros(56, 128); ones(16, 128); zeros(56, 128)];
%F = [ zeros(36, 128); ones(56, 128); zeros(36, 128)];
G = F';
H = F + 2 * G;

subplot(3,3,1)
showgrey(F)
title('F')
subplot(3,3,2)
showgrey(G)
title('G')
subplot(3,3,3)
showgrey(H)
title('H')

Fhat = fft2(F);
Ghat = fft2(G);
Hhat = fft2(H);

subplot(3,3,4)
showgrey(log(1 + abs(Fhat)));
title('log(1 + abs(Fhat))')
subplot(3,3,5)
showgrey(log(1 + abs(Ghat)));
title('log(1 + abs(Ghat))')
subplot(3,3,6)
showgrey(log(1 + abs(Hhat)));
title('log(1 + abs(Hhat))')

subplot(3,3,7)
showgrey(log(1 + abs(fftshift(Fhat))));
title('log(1 + abs(fftshift(Fhat)))')
subplot(3,3,8)
showgrey(log(1 + abs(fftshift(Ghat))));
title('log(1 + abs(fftshift(Ghat)))')
subplot(3,3,9)
showgrey(log(1 + abs(fftshift(Hhat))));
title('log(1 + abs(fftshift(Hhat)))')

Hhat2 = fft2(F) + 2*fft2(G);
subplot(1,2,1)
showgrey(log(1 + abs(fftshift(Hhat))));
title('fft2(F+2*G)')
subplot(1,2,2)
showgrey(log(1 + abs(fftshift(Hhat2))));
title('fft2(F) + 2*fft2(G)')

% figure
% surf(log(1 + abs(Fhat)))

% Square wave in spatial domain equals sinc function in frequency domain

%% Question 10
clear all, close all, clc
sz = 128;

F = [ zeros(56, 128); ones(16, 128); zeros(56, 128)];
G = F';

figure
subplot(1,3,1)
showgrey(F .* G);
title('F .* G');
subplot(1,3,2)
showfs(fft2(F .* G));
title('fft2(F .* G)');

Fhat = fft2(F);
Ghat = fft2(G);

subplot(1,3,3)
showgrey(log(1 + abs(conv2(fftshift(Fhat), fftshift(Ghat), 'same')/sz^2)));
title('conv2(Fhat, Ghat)')

% apply fftshift (matlab problem)
% showfs - auto centering. showgrey - have to have log(1 + abs(...))

%% Question 11, 12
% clear all, close all, clc

alpha = 30;

F = [zeros(60, 128); ones(8, 128); zeros(60, 128)] .* ...
    [zeros(128, 48) ones(128, 32) zeros(128, 48)];

figure
% subplot(2,3,1)
% showgrey(F);
% title('F')

% subplot(2,3,4)
% Fhat = fft2(F);
% showfs(Fhat);
% title('Fhat')

G = rot(F, alpha);

subplot(3,3,1)
showgrey(G)
title('G for alpha = 30')

subplot(3,3,4)
Ghat = fft2(G);
showfs(Ghat)
title('Ghat for alpha = 30')

subplot(3,3,7)
Hhat = rot(fftshift(Ghat), -alpha );
showgrey(log(1 + abs(Hhat)))
title('Hhat for alpha = 30')

alpha = 45;
G = rot(F, alpha);

subplot(3,3,2)
showgrey(G)
title('G for alpha = 45')

subplot(3,3,5)
Ghat = fft2(G);
showfs(Ghat)
title('Ghat for alpha = 45')

subplot(3,3,8)
Hhat = rot(fftshift(Ghat), -alpha );
showgrey(log(1 + abs(Hhat)))
title('Hhat for alpha = 45')

alpha = 60;
G = rot(F, alpha);

subplot(3,3,3)
showgrey(G)
title('G for alpha = 60')

subplot(3,3,6)
Ghat = fft2(G);
showfs(Ghat)
title('Ghat for alpha = 60')

subplot(3,3,9)
Hhat = rot(fftshift(Ghat), -alpha );
showgrey(log(1 + abs(Hhat)))
title('Hhat for alpha = 60')

% txt = ['alpha = ',num2str(alpha)];
% sgtitle(txt)

% figure
% subplot(1,2,1)
% showgrey(F);
% title('F')
% subplot(1,2,2)
% Fhat = fft2(F);
% showfs(Fhat);
% title('Fhat')

%% Question 13
clear all, close all, clc

a = 10e-10;
phone = phonecalc128;
few = few128;
nallo = nallo128;

% Raw images
figure
subplot(3,3,1)
showgrey(phone)
title('Raw images')
subplot(3,3,2)
showgrey(few)
subplot(3,3,3)
showgrey(nallo)


%  pow2image performs a transformation in the Fourier domain such
%  that the phase information is preserved, whereas the magnitude
%  is REPLACED BY a power spectrum of the form

subplot(3,3,4)
showgrey(pow2image(phone, a));
title('Phase preserved, magnitude changed')
subplot(3,3,5)
showgrey(pow2image(few, a));
subplot(3,3,6)
showgrey(pow2image(nallo, a));


% randphaseimage keeps the magnitude but randomises the distribution
% of the phase.
subplot(3,3,7)
showgrey(randphaseimage(phone))
title('Magnitude preserved, phase randomized')
subplot(3,3,8)
showgrey(randphaseimage(few))
subplot(3,3,9)
showgrey(randphaseimage(nallo))


%% Question 14
clear all, close all, clc

t = [0.1 0.3 1.0 10.0 100.0];                                          % Variance
pic = few128;

i = 0;
figure
for i = 1:length(t)
    psf = gaussfft(deltafcn(128, 128), t(i));          % impulse signal. Impulse response = after the filter
    varPsf = variance(psf);
    subplot(3,2,i)
    showgrey(psf)
    title(['t = ', num2str(t(i))])
end
sgtitle('Impulse response, 2D')

i = 0;
figure
for i = 1:length(t)
    psf = gaussfft(deltafcn(128, 128), t(i));          % impulse signal. Impulse response = after the filter
    subplot(3,2,i)
    surf(psf)
    title(['t = ', num2str(t(i))])
end
sgtitle('Impulse response, 3D')


t = [1.0 4.0 16.0 64.0 256.0];
i = 0;
figure
for i = 1:length(t)
    C = gaussfft(pic, t(i));          % impulse signal. Impulse response = after the filter
    subplot(3,2,i)
    showgrey(C)
    title(['t = ', num2str(t(i))])
end
sgtitle('Gaussian filter on images')


% C = gaussfft(pic, t);

%var = variance(C);

% figure
% subplot(2, 1, 1)
% showgrey(pic)
% subplot(2, 1, 2)
% showgrey(C)

% when it is small it will not change as much.
% Convolve the impulse signal as an image. image =

%% Question 17-18
clear all, close all, clc

office = office256;
add = gaussnoise(office, 16);           % Adds white Gaussian noise with
% with standard deviation SDEV to
% INPIC. (gaussnoise(INPIC, SDEV)

sap = sapnoise(office, 0.1, 255);       % Adds salt-and-peppar noise to an
% image by resetting a fraction
% FRAC/2 to ZMIN and a similar
% fraction to ZMAX in a
% pixel-to-pixel independant
% manner. With no ZMIN and ZMAX
% values they are set to the true
% minimum and maximum values of
% inpic. (sapnoise(inpic, FRAC, ZMIN, ZMAX)

figure
subplot(1, 2, 1)
showgrey(add)
title('Gauss noise')
subplot(1, 2, 2)
showgrey(sap)
title('salt-and-peppar noise')
sgtitle('Noisy images')

t = [0.1, 0.3, 1.0, 10.0, 100.0];

figure
for i = 1:length(t)
    GSadd = gaussfft(add, t(i));           % Gaussian smoothing of add
    subplot(3, 2, i)
    showgrey(GSadd)
    title(['t = ', num2str(t(i))])
end
sgtitle('Gaussian smoothing of Gauss noise')

figure
for i = 1:length(t)
    GSsap = gaussfft(sap, t(i));           % Gaussian smoothing of sap
    subplot(3, 2, i)
    showgrey(GSsap)
    title(['t = ', num2str(t(i))])
end
sgtitle('Gaussian smoothing of salt-and-peppar noise')

width = [1:3:13];

figure
for i = 1:length(width)
    medadd = medfilt(add, width(i));           % Median fitlering of add
    subplot(3, 2, i)
    showgrey(medadd)
    title(['Width = ', num2str(width(i))])
end
sgtitle('Median fitlering of Gauss noise')

figure
for i = 1:length(width)
    medsap = medfilt(sap, width(i));           % Median filtering of sap
    subplot(3, 2, i)
    showgrey(medsap)
    title(['Width = ', num2str(width(i))])
end
sgtitle('Median fitlering of salt-and-peppar noise')

COfreq = [0.01, 0.05, 0.1, 0.5, 1];
figure
for i = 1:length(COfreq)
    LPadd = ideal(add, COfreq(i));           % Low-pass filter of add
    subplot(3, 2, i)
    showgrey(LPadd)
    title(['cut-off frequency = ', num2str(COfreq(i))])
end
sgtitle('Low-pass filter of Gauss noise')

figure
for i = 1:length(COfreq)
    LPsap = ideal(sap, COfreq(i));           % Low-pass filter of sap
    subplot(3, 2, i)
    showgrey(LPsap)
    title(['cut-off frequency = ', num2str(COfreq(i))])
end
sgtitle('Low-pass filter of salt-and-peppar noise')

%% Question 19, 20
%clear all, close all, clc

figure
t = 2;
freq = 0.3;
img = phonecalc256;
GFsmoothimg = img;
LPsmoothimg = img;
N=5;
for i=1:N
    if i>1 % generate subsampled versions
        img = rawsubsample(img);
        GFsmoothimg = gaussfft(GFsmoothimg, t);
        GFsmoothimg = rawsubsample(GFsmoothimg);
        LPsmoothimg = ideal(LPsmoothimg, freq);
        LPsmoothimg = rawsubsample(LPsmoothimg);
    end
    subplot(3, N, i)
    showgrey(img)
    title([num2str(length(img)),'^{2} pixels'])
    subplot(3, N, i+N)
    showgrey(GFsmoothimg)
    title([num2str(length(img)),'^{2} pixels'])
    subplot(3, N, i+2*N)
    showgrey(LPsmoothimg)
    title([num2str(length(img)),'^{2} pixels'])
end
sgtitle(['t = ', num2str(t), ', freq = ', num2str(freq)])