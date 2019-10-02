function llend = llwalk(llstart,latdist,londist)
% function llend = llwalk(llstart,latdist,londist)
% Given a point on the globe, calculates geographics coordinates for another point at a given x and
% y distance away from that point.
% Calculation is based on iterative minimization of a guessed coordinate using lldistkm.
%
% INPUTS:
% llstart: starting lat/lon coordinate.
% latdist: distance along latitude coordinate, km.  Pos = N, neg = S.
% londist: distance along longitude coordinate, km. Pos = E, neg = W.
%
% OUTPUTS:
% llend: lat/lon coordinate for ending point.
%
% 20190708 GMW

x_coord = walk(llstart,abs(londist),[0 sign(londist)]); % longitude
y_coord = walk(llstart,abs(latdist),[sign(latdist) 0]); % latitude

llend = [y_coord(1) x_coord(2)];

% function to calculate lat-lon pair some distance from a starting point
% llstart is the lat-lon pair
% dist is the distance to go in km
% direction says which way to go: N=[1 0], E=[0 1], S=[-1 0], W=[0 -1]
function llend = walk(llstart,dist,direction)
    ddeg = 0.01;
    dx = lldistkm(llstart,llstart + ddeg.*direction); % estimate distance per degree
    dxError = abs(dx - dist);
    Tol = 0.001; %position error tolerance, km
    rep = 0;
    while dxError>Tol
        ddeg = ddeg./dx.*dist;
        dx = lldistkm(llstart,llstart + ddeg.*direction);
        dxError = abs(dx - dist);
            
        rep = rep+1;
        if rep>1000
            disp('llbox: max repetitions exceeded before convergence.')
            break
        end
    end
    llend = llstart + ddeg.*direction;

