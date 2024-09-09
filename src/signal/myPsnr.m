function [psnr] = myPsnr(x)

psnr = 10*log10(max(x).^2./var(x));
% [mx, ix] = max(x);
% x(ix)=[];
% psnr = 10*log10(mx^2/mean(abs(x).^2));

end