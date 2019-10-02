function [boxLat,boxLon] = llbox(center,width,height)
% function [boxLat,boxLon] = llbox(center,width,height)
% calculates the bounds of a geographic box centered on some point.
% INPUTS:
% center: latitude-longitude pair specifying centerpoint of box.
% width: full width (E-W distance) of box, km.
% height: full height (N-S distance) of box, km.
%
% OUTPUTS:
% boxLat: S and N edges of box
% boxLon: W and E edges of box
%
% 20170110 GMW

llE = walk(center,width/2,[0 1]);
llW = walk(center,width/2,[0 -1]);
llN = walk(center,height/2,[1 0]);
llS = walk(center,height/2,[-1 0]);

boxLat = [llS(1) llN(1)];
boxLon = [llW(2) llE(2)];


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


