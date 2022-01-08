function [yout,youtstd,n,xout] = BinAvg(xin,yin,x2avg,num_pts,num_std,meanmed,som)
% function [yout,youtstd,n,xout] = BinAvg(xin,yin,x2avg,num_pts,num_std,meanmed,som)
%
% Bins and averages data.
% INPUTS:
% xin: input independent variable (vector).
% yin: input data. If a matrix, each column is assumed to be a variable, and rows correspond to xin.
% x2avg: Value specifying how to average data:
%        1) If a scalar, this specifies the spacing of the averaging vector xout. 
%           In this case, xout = min(xin):x2avg:max(xin).
%        2) If a vector, this gives the x-values to which yin will be averaged.
%           The bin width is assumed equal to the median spacing of adjacent x2avg values.
%        3) If a 2-D array, this is assumed to be a 2-column matrix with start and stop times.
%           In this case, bin edges are x2avg(i,1) to x2avg(i,2), and xout is the mid-point of each window.
%
% OPTIONAL INPUTS:
% num_pts: minimum number of points in a bin. Bins with less than this number of points will return
%           NaNs as averages and std.
% num_std: number of standard deviations for outlier filter. Points lying outside of 
%           a window defined by (mean +/- std) in a bin will be excluded from the average.
% meanmed: flag for averaging. 0=means, 1=medians, 2=max, 3 = min;
% som:     flag for standard deviations. 0=normal std, 1=std of mean.
%
% OUTPUTS:
% yout: bin-averaged data.
% youtstd: standard deviations for each bin
% n: number of valid points in each bin
% xout: averaging vector.
%
% 20120627 GMW
% 20120803 GMW    Fixed issue with indexing for repeated bins (e.g. with diel averaging).
% 20131002 GMW    Addeded cabability to input time windows for xavg.
% 20200318 GMW    Added option to do min or max of bin (meanmed flag)

%%%%%DEFAULT INPUTS%%%%%
if nargin<4, num_pts=0; end
if nargin<5, num_std=0; end
if nargin<6, meanmed=0; end
if nargin<7, som=0; end

%%%%%GENERATE X VECTORS%%%%%
if isscalar(x2avg)
    inc = x2avg;
    xedge = (min(xin):inc:max(xin))'; %bin edges
    xout = xedge + inc/2; %bin centers
elseif isvector(x2avg)
    inc = median(diff(x2avg),'omitnan'); %use bin spacing to determine increment
    xedge = x2avg - inc/2;
    xout = x2avg;
else %assume start and stop windows given
    nwin = size(x2avg,1); %number of windows
    xedge = nan(2*nwin,1);
    xedge(1:2:2*nwin-1) = x2avg(:,1);%interleave start and stop values
    xedge(2:2:2*nwin) = x2avg(:,2);
    xout = mean(x2avg,2,'omitnan'); %use bin centers
    inc = max(diff(xedge)); %prevents bins from being thrown away for wrong size
end

nbins = length(xedge);
xedge(end+1) = xedge(end)+inc; %need to define outer edge of last bin

%%%%%GENERATE BIN INDICES%%%%%
%ptsperbin gives number points in each bin; bincol is an index for which bin each point goes into
[ptsperbin,bincol] = histc(xin,xedge);
ptsperbin(nbins+1,:) = []; %remove extra "bin" generated by specifying outer edge of last bin
bincol(bincol==nbins+1) = 0;

%remove in-between bins if window x2avg was given as a 2-column window
if exist('nwin','var')
    ptsperbin(2:2:nbins) = [];
    i = rem(bincol,2)==0; %even indices
    bincol(i) = 0;
    bincol(~i) = (bincol(~i)+1)/2; %go from 1,3,5,7 to 1,2,3,4
end
    

%deal with gaps in xedge if it was input as a vector
dx = diff(xedge);
i = find(dx>inc*1.01); %eliminate bins bigger than 1% of "median" bin size
j = ismember(bincol,i);
bincol(j) = 0;
ptsperbin(i)=0;

%remove out-of-range points
i = bincol==0;
bincol(i)=[];
yin(i,:)=[]; 

nbin = length(ptsperbin); %number of bins
maxbin = max(ptsperbin); %max number of points in a bin

%complimentary index for mapping to averaging matrix
binrow = nan(size(bincol));
if any(diff(bincol)<0) %slower but safer for non-monotonic series
    for i=1:nbin
        if ptsperbin(i)==0, continue; end
        binrow(bincol==i) = 1:ptsperbin(i);
    end
else
    j=1;
    for i=1:nbin
        if ptsperbin(i)==0, continue; end
        k = j:j+ptsperbin(i)-1;
        binrow(k) = 1:ptsperbin(i);     %faster but only works for monotonic series
        j = k(end)+1;
    end
end

i2map = sub2ind([maxbin nbin],binrow,bincol); %linear index for mapping onto averaging matrix

%%%%%DO AVERAGING%%%%%
ncol = size(yin,2); %number of columns
yout = nan(nbin,ncol);
youtstd = nan(nbin,ncol);
n = nan(nbin,ncol);
for i=1:ncol
    ynow = yin(:,i);
    y2avg = nan(maxbin,nbin);
    y2avg(i2map) = ynow;

    %%%%%CALCULATE STATISTICS%%%%%
    switch meanmed
        case 0
            yavg = mean(y2avg,1,'omitnan');
        case 1
            yavg = median(y2avg,1,'omitnan');
        case 2
            yavg = max(y2avg,[],1,'omitnan');
        case 3
            yavg = min(y2avg,[],1,'omitnan');
    end
    
    if som, ystd = stdom(y2avg);
    else    ystd = std(y2avg,1,'omitnan');
    end

    %%%%%FILTER OUTLIERS%%%%%
    if num_std>0
        lowerlimit = ones(maxbin,1)*(yavg - num_std*ystd);
        upperlimit = ones(maxbin,1)*(yavg + num_std*ystd);
        out = y2avg>upperlimit | y2avg<lowerlimit;
        y2avg(out) = nan;
        
    switch meanmed
        case 0
            yavg = nanmean(y2avg,1);
        case 1
            yavg = nanmedian(y2avg,1);
        case 2
            yavg = nanmax(y2avg,[],1);
        case 3
            yavg = nanmin(y2avg,[],1);
    end

        if som, ystd = stdom(y2avg,1);
        else    ystd = nanstd(y2avg,1);
        end
    end

    %%%%%FILTER FOR MINIMUM POINTS%%%%%
    nvalid = sum(~isnan(y2avg),1);
    yavg(nvalid<num_pts) = nan;
    ystd(nvalid<num_pts) = nan;
    
    %%%%%ASSIGN OUTPUTS%%%%%
    yout(:,i) = yavg;
    youtstd(:,i) = ystd;
    n(:,i) = nvalid;
end


