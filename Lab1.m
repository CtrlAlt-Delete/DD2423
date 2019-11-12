% addpath('C:\MATLABPATH\DD2423_Lab_Files\Functions','C:\MATLABPATH\DD2423_Lab_Files\Images','C:\MATLABPATH\DD2423_Lab_Files\Images-m','C:\MATLABPATH\DD2423_Lab_Files\Images-mat')

%% Question 1
clear all, close all, clc

figure
p = 5;
q = 9;
sz = 128;
fftwavetemplate(p, q, sz)

% figure
% p = 9;
% q = 5;
% fftwavetemplate(p, q, sz)
%
% figure
% p = 17;
% q = 9;
% fftwavetemplate(p, q, sz)
%
% figure
% p = 17;
% q = 121;
% fftwavetemplate(p, q, sz)
%
% figure
% p = 5;
% q = 1;
% fftwavetemplate(p, q, sz)
%
% figure
% p = 125;
% q = 1;
% fftwavetemplate(p, q, sz)

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

%% Question 3 - 4
%function fftwave(u, v, sz)
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

%% Question 5 - 6
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

%% Question 7 - 9
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

% figure
% surf(log(1 + abs(Fhat)))

% Square wave in spatial domain equals sinc function in frequency domain

%% Question 10
clear all, close all, clc
F = [ zeros(56, 128); ones(16, 128); zeros(56, 128)];
G = F';

figure
subplot(3,1,1)
showgrey(F .* G);
title('F .* G');
subplot(3,1,2)
showfs(fft2(F .* G));
title('fft2(F .* G)');

Fhat = fft2(F);
Ghat = fft2(G);

subplot(3,1,3)
showfs(conv2(Fhat, Ghat));
title('conv2(Fhat, Ghat)')

%% Question 11 - 12
clear all, close all, clc

alpha = 30;

F = [zeros(60, 128); ones(8, 128); zeros(60, 128)] .* ...
    [zeros(128, 48) ones(128, 32) zeros(128, 48)];

figure
subplot(2,3,1)
showgrey(F);
title('F')

subplot(2,3,4)
Fhat = fft2(F);
showfs(Fhat);
title('Fhat')

G = rot(F, alpha);

subplot(2,3,2)
showgrey(G)
title('G')
axis on

subplot(2,3,5)
Ghat = fft2(G);
showfs(Ghat)
title('Ghat')

subplot(2,3,6)
Hhat = rot(fftshift(Ghat), -alpha );
showgrey(log(1 + abs(Hhat)))
title('Hhat')

txt = ['alpha = ',num2str(alpha)];
sgtitle(txt)

%% Question 13
clear all, close all, clc

a = 10e-10;
phone = phonecalc128;
few = few128;
nallo = nallo128;

subplot(2,3,1)
image(phone)
subplot(2,3,2)
image(few)
subplot(2,3,3)
image(nallo)
title('Raw images')

figure
subplot(2,3,4)
showgrey(pow2image(phone, a));
subplot(2,3,5)
showgrey(pow2image(few, a));
subplot(2,3,6)
showgrey(pow2image(nallo, a));
title('Phase preserved, magnitude changed')

figure
subplot(2,3,1)
image(randphaseimage(phone))
subplot(2,3,2)
image(randphaseimage(few))
subplot(2,3,3)
image(randphaseimage(nallo))
title('Magnitude preserved, phase randomized')

%% Question 14-16
clear all, close all, clc

t = 10;                                          % Variance
pic = few128;

C = gaussfft(pic, t);

psf = gaussfft(deltafcn(128, 128), t);          % impulse signal. Impulse response = after the filter 
figure
showgrey(psf)
title('Impulse response')

var = variance(C);

figure 
subplot(2, 1, 1)
showgrey(pic)
subplot(2, 1, 2)
showgrey(C)

% when it is small it will not change as much.
% Convolve the impulse signal as an image. image = 




