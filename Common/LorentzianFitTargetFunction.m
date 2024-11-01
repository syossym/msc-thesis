function y = LorentzianFitTargetFunction(params, E)

global E_0 Amp;

N = length(params(1:end))/3;
E_0_vec = params(1:N);
amp_vec = params(N+1:2*N);
fwhm_vec = params(2*N+1:end);

%E_0_vec = params;

y = 0;
for (ii=1:length(E_0_vec))
    y = y + Lorentz(E,E_0_vec(ii),amp_vec(ii),fwhm_vec(ii));
    %y = y + Lorentz(E,E_0(ii),Amp(ii),params);
end

%figure(1); box on; hold on;
%plot(E, y, 'r');

function out = Lorentz(x,x0,amplitude,fwhm)

out = amplitude*(fwhm/(2*pi))*((x-x0).^2+(fwhm/2)^2).^(-1);