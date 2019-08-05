function P = plumeAverage(D,legTimes,vars2fit,vars2avg,timeVarName,fitVarName)
% function P = plumeAverage(D,legTimes,vars2fit,vars2avg,timeVarName,fitVarName)
% Performs a few simple calculations on defined plumes.
%
% INPUTS:
% D: structure of data.
% legTimes: 2-column matrix of start and stop times defining each plume.
% vars2fit: cell array of names of variables for calculating enhancment ratios (slope relative to fitVar).
%   Note, variables are lagged relative to fitVar for each plume prior to fitting.
% vars2avg: cell array of names of variables to average over whole plume.
% timeVarName: name of time variable (to use with legTimes for indexing plumes).
% fitVarName: name of variable to fit against (e.g., CO).
%
% OUTPUTS:
% P: structure containing info for each plume.
%     P.ER: enhanement ratios (for vars2fit)
%     P.ER.r: correlation coefficients for fits
%     P.AV: averages (for vars2avg)
%
% 20190802 GMW

% get some variables
time = D.(timeVarName);
xfit = D.(fitVarName);

% initialize outputs
n = nan(length(legTimes),1);
P = struct;
for i = 1:length(vars2fit)
    P.ER.(vars2fit{i}) = n; %slopes
    P.ER.r.(vars2fit{i}) = n; %correlation coefficients
end

for i = 1:length(vars2avg)
    P.AV.(vars2avg{i}) = n;
end

% loop through legs
for i = 1:length(legTimes)
    j = time>=legTimes(i,1) & time<=legTimes(i,2);
    
    for k = 1:length(vars2fit)
        x = xfit(j);
        y = lagit(x,D.(vars2fit{k})(j),5,0); %assumes timing good to w/in 5 seconds
        g = ~isnan(x+y);
    	[P.ER.(vars2fit{k})(i),~,P.ER.r.(vars2fit{k})(i)] = lsqfitgm(x(g),y(g));
    end
    
    for k = 1:length(vars2avg)
        P.AV.(vars2avg{k})(i) = nanmean(D.(vars2avg{k})(j));
    end
    
end


