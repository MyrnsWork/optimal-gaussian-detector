%  Parseval theorem [1]:
%  -------------------
%
%  In this tutorial we show that the total energy of waveform X(t) (1D,2D and 3D) 
%  computed in time domain is equal to the total energy of the waveform's Fourier Transform
%  F(X(t))=x(f) computed in the frequency domain . 
% Reference :
%[1] :  Parseval des Chênes, Marc-Antoine "Mémoire sur les séries et sur
%       l'intégration complète d'une équation aux différences partielles
%       linéaire du second ordre,à coefficients constants" Paris 1799.
%
% (c)  KHMOU Youssef,  Applied Mathematics 

close all;

% 1D case :  Sinusoidal function 
Fs=40;f=4;Ts=1/Fs;T=2;t=0:Ts:T-Ts;N=length(t);
x=2*cos(2*pi*f*t);
fx=fft(x)/N;
figure,
subplot(1,2,1), area(t,abs(x.^2)),title(' Time Domain');
subplot(1,2,2),area(abs(fx)), title(' Frequency Domain');
E1_timedomain  = sum(abs(x.^2))/N
E1_frequdomain = sum(abs(fx.^2))

%  2D case : Sinusoidal function
[X,Y]=meshgrid(0:Ts:T);R=sqrt(X.^2+Y.^2);
y=sin(R*f);
fy=fftshift(fft2(y));
figure,
subplot(1,2,1),surf(abs(y.^2)),shading interp, title(' Time Domain');
subplot(1,2,2), surf(abs(fy)), shading interp, title(' Frequency Domain');
E2_timedomain=sum(sum(abs(y.^2)))
E2_frequdomain=sum(sum(abs(fy.^2)))/N^2

% 3D case :  Gaussian function
z=randn(N,N,N);
fz=fftn(z);
E3_timedomain=sum(sum(sum(abs(z.^2))))
E3_frequdomain=sum(sum(sum(abs(fz.^2))))/N^3