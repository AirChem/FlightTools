function age = LGage(llstart,llstop,wspd)
% function age = LGage(llstart,llstop,wspd)
% Calculates lagrangian age given start and stop coordinates and wind speed.
% INPUTS:
% llstart: 2-column matrix of starting lat/lon (e.g, source point).
% llstop: 2-column matrix of stop lat/lons.
% wspd: vector of wind speeds, m/s.
%
% OUTPUTS:
% age: lagrangian age (distance/wind speed), hours.
%
% 20190802 GMW

% extend start if needed
nrows = size(llstop,1);
if size(llstart,1)==1
    llstart = repmat(llstart,nrows,1);
end

dist = lldistkm(llstart,llstop); %get distances in km
age = dist*1000./wspd/3600;