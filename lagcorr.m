% function [r,lags,dt] = lagcorr(x,y,maxlag,plotme,t)
% Calculated a lagged correlation coefficients between two vectors.
% The function lags y relative to x from -maxlag:maxlag.
% Also makes a lag-correlation plot.
% INPUTS:
% x,y: vectors to be compared. Should be column vectors.
% maxlag: maximum number of points to lag y relative to x.
% plotme: flag for making lag-correlation plot. default is 1 (yes).
% t: time stamps for x and y, used to determine time-step of each lag (optional).
%    t MUST be evenly-spaced for this to work properly.
% OUTPUTS:
% r: correlation coefficient of x and y at each lag.
% lags: # of lag points for each r. Positve values correspond to y being shifted forward.
% dt: corresponding delta-t of each lag.
%
% 090717 GMW
%
% Improved to deal better with NaN-containing vectors and calculate time-steps.
% 120617 GMW

function [r,lags,dt] = lagcorr(x,y,maxlag,plotme,t)

Ly = length(y);

%%%%%DEFAULTS%%%%%
if nargin<4
    plotme = 1;
end

if nargin<5
    t = nan(Ly,1);
end

%%%%%CALCULATE LAGS%%%%%
lags = -maxlag:maxlag;
r = nan(length(lags),1);
dt = nan(length(lags),1);
for i=1:length(lags)
    l = lags(i);
    if l>0
        yl = [nan(l,1); y(1:Ly-l)];
        dt(i) = t(l+1) - t(1);
    elseif l<0
        l=-l;
        yl = [y(l+1:end); nan(l,1)];
        dt(i) = t(end-l) - t(end);
    else
        yl = y;
        dt(i) = 0;
    end
    
    rnow = corrcoef(x,yl,'rows','pairwise');
    r(i) = rnow(1,2);
end

if plotme
    figure;
    if nargin==5
        plot(dt,r,'.-')
        xlabel('\Delta t_y')
    else
        plot(lags,r,'.-')
        xlabel('y lag (# of points)')
    end
    ylabel('r_x_y')
    grid on
end


% beginnings of a faster method, if you had enough memory
% ypad = [nan(maxlag,1); x(:); nan(maxlag,1)];
% imat = repmat((1:length(ypad))',[1 length(ypad)]);
% index = spdiags(imat,-maxlag:maxlag);


