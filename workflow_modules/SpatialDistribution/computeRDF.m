%% Oliver Grimm and Lucy Hutchinson 2021

function [Lnorm, RD, D] = computeRDF(tile_curr, radius_max)
%% Calculate the RDF with periodic BCs for one tile and output Lnorm, Lstd,
% Get area of square ring
sq_a = 3:2:radius_max*2 +1; sq_b = 1:2:(radius_max*2)-1;
A = sq_a.^2-sq_b.^2;
num_CD8 = sum(sum(tile_curr));

%% Break if too few CD8s
if num_CD8 < 10 || num_CD8 > 1800
    Lnorm = nan(1,radius_max);
    RD = nan(1,radius_max);
    D = nan;
    return
end


% Replicate tile so it is a 3*3 repetition of the same tile
% (for periodicity)
tileDim = size(tile_curr,1);
tile_ext = repmat(tile_curr,3,3);
% Get the positions of the cells in the centre tile
[CD8_x, CD8_y] = find(tile_curr); CD8 = [CD8_x CD8_y];
% Get the positions of the cells in the extended tile
[CD8ext_x, CD8ext_y] = find(tile_ext); CD8ext = [CD8ext_x CD8ext_y];
% Want to get the RDF for the cells in the centre tile
CD8_centre = CD8+tileDim;

[~,D] = knnsearch(CD8ext,CD8_centre,'K',4*radius_max*radius_max,'Distance','chebychev');
D(D>radius_max)=[]; 
D(D==0)=[]; 
RD = histcounts(D,[1:radius_max+1]);

cell_num_density = 2*4*[1:radius_max]*num_CD8*(num_CD8-1)/(tileDim^2-1);% Factor of 2 from the fact every distance is measured twice
Lnorm = RD./cell_num_density;
