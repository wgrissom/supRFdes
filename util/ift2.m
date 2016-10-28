function y = ift2(x)
    y = ifftshift(ifft2(ifftshift(x)));
end