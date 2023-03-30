%Demonstrating Weiner Filter using Weiner-Hopf equations
% Code to remove a noisy portion of a corrupted signal, back to the
% reference signal

clear
close all
clc
fs = 6325; %frequency of samples
T = 3; %total recording tme 
L = T.*fs; %length of vector
tt = (0:L-1)/fs; %time vector
ff = (0:L-1).*fs / L;
y = 3*sin(2*pi*90.*tt); y = y(:); %sinusoidal reference signal
x = 0.70*randn(L,1) + y; x = x(:); % reference signal with noise added
N = 500;
[xest, b, MSE] = WeinerFilt(x,y,N);

% we then plots the results as below
figure
subplot(311);

plot(tt,x,'k'), hold on,  plot (tt,y,'r');
title('Weiner Filter Demo');

legend('Noisy signal', 'reference Signal');

subplot(312);
plot(tt(N+1:end), xest, 'k');
legend('Estimated signal');
subplot(313);

plot(tt(N+1:end), (x(N+1:end) - xest), 'k');
legend('Residue signal after removal of noise');
xlabel('time (s)')


%Define the WeinFilter function

function [xest,B,MSE] = WeinerFilt(x,y,N)
X = 1/N .* fft(x(1:N));
Y = 1/N .* fft(y(1:N));
X = X(:);
Y = Y(:);
Rxx = N .* real(ifft(X .* conj(X))); % Autocorrelation function
Rxy = N .* real(ifft(X .* conj(Y))); % Crosscorrelation function
Rxx = toeplitz(Rxx);
Rxy = Rxy';
B = Rxy / Rxx; B = B(:); % Wiener-Hopf eq. B = inv(Rxx) Rxy
xest = fftfilt(B,x);
xest = xest(N+1:end); % cut first N samples due to distorsion during filtering operation
MSE = mean(y(N+1:end) - xest) .^2; % mean squared error


end




