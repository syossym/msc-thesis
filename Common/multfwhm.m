function width_vec = multfwhm(x,y)

% function width = fwhm(x,y)
%
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
%
% Rev 1.2, April 2006 (Patrick Egan)

x_interp = min(x):(max(x)-min(x))/10000:max(x);
y = interp1(x, y, x_interp, 'pchip');
x = x_interp;
y = y / max(y);
N = length(y);
lev50 = 0.5;

[maxtab, mintab]=peakdet(y, 0.05);
width_vec = [];

for (ii=1:length(maxtab(:,1)))
    centerindex = maxtab(ii,1);
    
    i = 2;
    while sign(y(i)-lev50) == sign(y(i-1)-lev50)
        i = i+1;
    end                                   %first crossing is between v(i-1) & v(i)
    interp = (lev50-y(i-1)) / (y(i)-y(i-1));
    tlead = x(i-1) + interp*(x(i)-x(i-1));
    i = centerindex+1;                    %start search for next crossing at center
    while ((sign(y(i)-lev50) == sign(y(i-1)-lev50)) & (i <= N-1))
        i = i+1;
    end
    if i ~= N
        Ptype = 1;
        %disp('Pulse is Impulse or Rectangular with 2 edges')
        interp = (lev50-y(i-1)) / (y(i)-y(i-1));
        ttrail = x(i-1) + interp*(x(i)-x(i-1));
        width = ttrail - tlead;
    else
        Ptype = 2;
        %disp('Step-Like Pulse, no second edge')
        ttrail = NaN;
        width = NaN;
    end
    width_vec = [width_vec, width];
end
