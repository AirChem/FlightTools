function [ylagged,lag2use] = lagit(x,y,nlags,plotlags,method,setlag)
% function [ylagged,lagused] = lagit(x,y,nlags,plotlags,method)
% Applies a lag to a vector to align it with another vector.
% Useful for time-aligning flight data before scattering.
% The lag is determined by optimizing the correlation coefficient (r)
% between the two vectors.
%
% INPUTS:
% x: fixed vector.
% y: vector to be lagged.
% nlags: optional scalar specifying how many lags to look over.
%        Default value is 60.
% plotlags: flag for plotting lag-correlation plots
% method: optional flag for method to determine best lag. 0 = max
% correlation (DEFAULT), 1 = 2nd derivative.
% setlag: option to override lag2use with pre-defined lag, setlag.
%
% OUTPUTS:
% ylagged: y-vector lagged to optimize correlation with x.
%          Same size as y.
%          missing values are filled in as NaN.
% lag2use: lag applied.
%
% 20130610 GMW
% 20130717 GMW  Added flag for plotting lag-correlation and "lag" output
% 20130723 GMW  Modified to use 2nd derivative of r to identify best lag
% 20130806 GMW  Modified to only use 2nd derivative if r changes sign.
% 20171026 GMW  Added method input to force 2nd derivative or max method.
% 20171031 RAH  Added optional input setlag to override lag2use.

%defaults
y = y(:);

if nargin<6, setlag=[]; end
if nargin<5, method=0; end
if nargin<4, plotlags = 0; end
if nargin<3, nlags = 60; end

% calculate optimized lag shift
[r,lags] = lagcorr(x,y,nlags,0);
s = sign(r);
if length(unique(s))>1 || method==1 %if it changes sign, use gradient method
    del2r = gradient(gradient(r)); %second derivative
    [~,imax] = max(abs(del2r));
elseif length(unique(s))==1 || method==0 
    [~,imax] = max(abs(r)); %if not, use maximum
end

if ~isempty(setlag)
    lag2use = setlag;
else
    lag2use = lags(imax);
end

%build lagged vector
n = nan(abs(lag2use),1); %filler
if lag2use<0
    ylagged = [y(-lag2use+1:end); n]; %shift backward
elseif lag2use>0
    ylagged = [n; y(1:end-lag2use)]; %shift forward
else
    ylagged = y;
end

%plot if desired
if plotlags
    rnew = lagcorr(x,ylagged,nlags,0);
    figure
    subplot(2,2,[1 2])
    plot(lags,r,'b.-',lags,rnew,'r.-')
    xlabel('Lags')
    ylabel('r_x_y')
    legend('Original','Lagged')
    title(['Optimal lag = ' num2str(lag2use) ' points'])
    grid on
    
    subplot(2,2,3)
    plot(x,y,'.')
    xlabel('x')
    ylabel('y')
    title('Original')
    
    subplot(2,2,4)
    plot(x,ylagged,'.')
    xlabel('x')
    ylabel('y')
    title('Lagged')
end


