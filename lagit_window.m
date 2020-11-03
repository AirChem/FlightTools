function [ylagged,ylagged_avg,lag2use] = lagit_window(x,tx,y,ty,nlags,plotlags,method,setlag)
% function [ylagged,lagused] = lagit_window(x,tx,y,ty,nlags,plotlags,method)
% Applies a lag to a vector to align it with another vector.
% Useful for time-aligning flight data before scattering.
% The lag is determined by optimizing the correlation coefficient (r)
% between the two vectors.
%
% This version differs from lagit.m in that it is designed to work with two variables on different
% timebases, one "fast" and the other "slow." The latter is typically an averaged variable, like WAS
% or TOGA data.
%
% INPUTS:
% x: fixed vector. This must be the slower of the two datasets.
% tx: time for x. Must be a 2-column matrix of start and stop times.
% y: vector to be lagged.
% ty: time for y. Must be a 1-column matrix of start times.
% nlags: optional scalar specifying how many lags to look over.
%        Default value is 60.
% plotlags: flag for plotting lag-correlation plots
% method: optional flag for method to determine best lag. 0 = max
% correlation (DEFAULT), 1 = 2nd derivative.
% setlag: option to override lag2use with pre-defined lag, setlag.
%
% OUTPUTS:
% ylagged: y-vector lagged to optimize correlation with x.
%          Same size as y. Ends padded with NaNs.
% ylagged_avg: lagged y, averaged to input tx.
% lag2use: lag applied.
%
% 20130610 GMW
% 20130717 GMW  Added flag for plotting lag-correlation and "lag" output
% 20130723 GMW  Modified to use 2nd derivative of r to identify best lag
% 20130806 GMW  Modified to only use 2nd derivative if r changes sign.
% 20171026 GMW  Added method input to force 2nd derivative or max method.
% 20171031 RAH  Added optional input setlag to override lag2use.
% 20200318 GMW  Modified from lagit.m


%% Check inputs
assert(size(tx,2)==2,'Input tx must have two columns for start and stop times.')
assert(size(ty,2)==1,'Input ty must have one column of start times.')

%% Do stuff
%defaults
y = y(:);

if nargin<8, setlag=[]; end
if nargin<7, method=[]; end
if nargin<6, plotlags = 0; end
if nargin<5, nlags = 60; end

% calculate optimized lag shift
lags_fast = -nlags:nlags;
r_fast = nan(size(lags_fast));
for i = 1:length(lags_fast)
    
    y_shifted = lagVar(y,lags_fast(i)); % shift fast
    y_avg = BinAvg(ty,y_shifted,tx); % average fast to slow
    [r,~] = lagcorr(x,y_avg,nlags,0); % do lag correlation
    
    % pick optimum shift
    s = sign(r);
    if method==0 || length(unique(s))==1
        [~,imax] = max(abs(r)); %default to maximum if well-behaved
    else % length(unique(s))>1 || method==1 %if it changes sign, use gradient method
        del2r = gradient(gradient(r)); %second derivative
        [~,imax] = max(abs(del2r));
    end
    
    r_fast(i) = r(imax);
end

% pick optimum shift for all fast shifts
r = r_fast;
lags = lags_fast;
s = sign(r);
if method==0 || length(unique(s))==1
    [~,imax] = max(abs(r)); %default to maximum if well-behaved
else % length(unique(s))>1 || method==1 %if it changes sign, use gradient method
    del2r = gradient(gradient(r)); %second derivative
    [~,imax] = max(abs(del2r));
end

% choose lag
if ~isempty(setlag)
    lag2use = setlag;
else
    lag2use = lags(imax);
end

%build lagged y
n = nan(abs(lag2use),1); %filler
if lag2use<0
    ylagged = [y(-lag2use+1:end); n]; %shift backward
elseif lag2use>0
    ylagged = [n; y(1:end-lag2use)]; %shift forward
else
    ylagged = y;
end
ylagged_avg = BinAvg(ty,ylagged,tx); % average fast to slow

%plot if desired
if plotlags
    figure
    subplot(2,2,[1 2])
    plot(lags,r,'b.-')
    xlabel('Lags')
    ylabel('r_x_y')
    title(['Optimal lag = ' num2str(lag2use) ' points'])
    grid on
    
    subplot(2,2,3)
    yold_avg = BinAvg(ty,y,tx);
    plot(x,yold_avg,'.')
    xlabel('x')
    ylabel('y')
    title('Original')
    
    subplot(2,2,4)
    plot(x,ylagged_avg,'.')
    xlabel('x')
    ylabel('y')
    title('Lagged')
end


