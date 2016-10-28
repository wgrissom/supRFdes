function y = ft2(x)
    y = fftshift(fft2(fftshift(x)));
end