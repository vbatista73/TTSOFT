function [res] = Plotspec(xi,tau,d)
% Plot Spectrum
dtt=tau*2.418E-2; % femtosec
nsteps=size(xi,2);
TC=dtt*(nsteps-1);
if d==4
    % converting eV to rad/fs
    Es=(3.94+0.238663)*8065.0*3.0E-5*2.0*pi; % Expt. E(S1)
else
    Es=(4.06-1.7015+.11)*8065.0*3.0E-5*2.0*pi; %+0.0758
end
n0=nsteps;
nsteps=2^11;  % Apply phenomenological broadening
if d==24
    broadTau=150;
else
    broadTau=30;
end
tt=((1:n0)-1)*dtt;
xin=xi.*exp(-1i*Es*tt).*cos(pi*tt/2/TC).*exp(-tt/broadTau);
xin(n0+1:nsteps)=0;  % Extend xi with zeroes
dw=2*pi/(nsteps*tau);   % Build frequency grid
for j=1:nsteps
    if(j<(nsteps/2+1))
        w(j)=(j-1)*dw;
    else
        w(j)=(j-1-nsteps)*dw;
    end
end
sigma=real(ifft(xin));  % Compute spectrum as FT of xi
%*tau/sqrt(2*pi)*nsteps; % Normalization factor
E_ev=circshift(w,nsteps/2)*27.2; % Convert frequencies to eV
sigc=circshift(sigma,nsteps/2);
for j=1:nsteps
    sigc(j)=sigc(j)*E_ev(j)*nsteps*dtt;
end

res=[E_ev;sigc];

end

