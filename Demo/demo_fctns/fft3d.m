% Fourier transform function
% (transform along 3 spatial dimensions, at each time)

function [A] = fft3d(A)

A = fftshift3d(fft(fft(fft(ifftshift3d(A),[],1),[],2),[],3));

end