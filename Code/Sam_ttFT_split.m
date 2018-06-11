function x=Sam_ttFT_split(x,ftIFT)
% Performs FT/iFT on one physical dimension without quantics (so only 1
% core)
r=x.r;
n=x.n;
co=x.core;
co=reshape(co,[r(1) n(1) r(2)]);
co=permute(co,[2 1 3]);
if ftIFT==1
    co=fft(co);
else
    co=ifft(co);
end
co=permute(co,[2 1 3]);
x.core=co(:);
end