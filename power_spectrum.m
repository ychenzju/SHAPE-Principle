function [pows,freq] = power_spectrum(ts,fs,mode)

N = length(ts);

% take FFT and shift it for symmetry
amp = fftshift(fft(ts));

% make frequency range
fN = N - mod(N,2);
k = -fN/2 : fN/2+1;
T = N / fs;
freq = k/T;

% select the positive domain FFT and range
one_idx = fN/2 + 2;
amp = amp(one_idx:end);
freq = freq(one_idx:end);

%return power spectrum
pows = abs(amp).^2;
end
