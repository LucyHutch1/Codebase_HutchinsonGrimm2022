%% Edited by Roisin Stephens, Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, inspired by Jan Poleszczuk

function m = getAdjacent_2D_LH_OG(L,MYcells,nh)

% *** Periodic boundary conditions ***


dims = size(L);
N = dims(1);
M = dims(2);
m.S=[];
mS_index=[];

%find the extreme cases
MYcells_toprow = MYcells(MYcells>1 & MYcells<N);
MYcells_bottomrow = MYcells(MYcells>(N*(M-1)+1) & MYcells<N*M); %
MYcells_leftcol = MYcells(mod(MYcells,N) == 1 & MYcells>1 & MYcells<(N*(M-1)));
MYcells_rightcol = MYcells(mod(MYcells,N)==0 & MYcells>N & MYcells<(N*M));

[MYcells_topLcorner,iTL] = ismember(1,MYcells); %iTL - keeps track of ind in original array MYcells
[MYcells_topRcorner,iTR] = ismember(N,MYcells);
[MYcells_botLcorner,iBL] = ismember((N*(M-1))+1,MYcells);
[MYcells_botRcorner,iBR] = ismember(N*M,MYcells);

MYcells_middle = setdiff(MYcells, [MYcells_toprow,MYcells_bottomrow,MYcells_leftcol,MYcells_rightcol,[1,N,(N*(M-1))+1,N*M]],'stable');

% instead of concatenating the lists of cells for each domain area, we need
% to make the whole list at the end, and maintain the original order of the
% cells as in MYcells.

if numel(MYcells_middle)>0
    m.S = bsxfun(@plus,MYcells_middle,nh.aux(nh.Pms(:,randi(nh.nP,1,length(MYcells_middle)))));
    [vals,ia] = intersect(MYcells,MYcells_middle,'stable');
    ia=ia';
    mS_index=ia;
end

if MYcells_topLcorner
    nh_temp = int32([(N*M)-1 N*(M-1) (N*(M-1))+1 N-1 1 (2*N)-1 N N+1])';
    temp_var = bsxfun(@plus,1,nh_temp(nh.Pms(:,randi(nh.nP,1,1))));
    m.S = cat(2, m.S, temp_var);
    mS_index = cat(2,mS_index,iTL);
end

if MYcells_topRcorner
    nh_temp = int32([(N*M)-(N+1) (N*M)-N (N*(M-2))+1 -1 -(N-1) N-1 N 1])';
    temp_var = bsxfun(@plus,N,nh_temp(nh.Pms(:,randi(nh.nP,1,1))));
    m.S = cat(2, m.S, temp_var);
    mS_index = cat(2,mS_index,iTR);
end

if MYcells_botLcorner
    nh_temp = int32([-1 -N -N+1 N-1 1 -(N*(M-2))-1 -(N*(M-1)) -(N*(M-1))+1])';
    temp_var = bsxfun(@plus,(N*(M-1))+1,nh_temp(nh.Pms(:,randi(nh.nP,1,1))));
    m.S = cat(2, m.S, temp_var);
    mS_index = cat(2,mS_index,iBL);
end

if MYcells_botRcorner
    nh_temp = int32([-N-1 -N -(2*N)+1 -1 -N+1 -(N*(M-1))-1 -(N*(M-1)) -(N*M)+1])';
    temp_var = bsxfun(@plus,N*M,nh_temp(nh.Pms(:,randi(nh.nP,1,1))));
    m.S = cat(2, m.S, temp_var);
    mS_index = cat(2,mS_index,iBR);    
end

if numel(MYcells_toprow)>0
    nh_temp = int32([(N*(M-1))-1 N*(M-1) (N*(M-1))+1 -1 1 N-1 N N+1])';
    temp_var = bsxfun(@plus,MYcells_toprow,nh_temp(nh.Pms(:,randi(nh.nP,1,length(MYcells_toprow)))));
    m.S = cat(2, m.S, temp_var);
    mS_index = horzcat(mS_index,find(MYcells>1 & MYcells<N));
end

if numel(MYcells_bottomrow)>0
    nh_temp = int32([-N-1 -N -N+1 -1 1 (-N*(M-1))-1 -N*(M-1) (-N*(M-1))+1])';
    temp_var = bsxfun(@plus,MYcells_bottomrow,nh_temp(nh.Pms(:,randi(nh.nP,1,length(MYcells_bottomrow)))));
    m.S = cat(2, m.S, temp_var);
    mS_index = cat(2,mS_index,find(MYcells>(N*(M-1)+1) & MYcells<N*M));   
end

if numel(MYcells_leftcol)>0
    nh_temp = int32([-1 -N -N+1 N-1 1 (2*N)-1 N N+1])';
    temp_var = bsxfun(@plus,MYcells_leftcol,nh_temp(nh.Pms(:,randi(nh.nP,1,length(MYcells_leftcol)))));
    m.S = cat(2, m.S, temp_var);
    mS_index = cat(2,mS_index,find(mod((MYcells+N),N) == 1 & MYcells>1 & MYcells<(N*(M-1))));   
end

if numel(MYcells_rightcol)>0
    nh_temp = int32([-N-1 -N (-2*N)+1 -1 -N+1 N-1 N 1])';
    temp_var = bsxfun(@plus,MYcells_rightcol,nh_temp(nh.Pms(:,randi(nh.nP,1,length(MYcells_rightcol)))));
    m.S = cat(2, m.S, temp_var);
    mS_index = cat(2,mS_index,find(mod(MYcells,N)==0 & MYcells>N & MYcells<(N*M)));   
end
%}
%last column: if MYcells is multiple of N
%first column: if (MYcells+N) is multiple of N+1
%top row: if MYcells<N+1
%bottom row: if MYcells>N*(M-1)
%top left corner: if MYcells=1
%top right corner: if MYcells=N
%bottom left corner: if MYcells = N*(M-1) + 1
%bottom right corner: if MYcells = M*N

% Need to sort m.S according to the indicies in mS_index
[mS_index,indicies2sort_mS] = sort(mS_index,'ascend');
m.S = m.S(:,indicies2sort_mS);
m.S(L(m.S)) = 0; 			% setting occupied grid cells to false
m.indxF = find(any(m.S)); 	% selecting agents with at least one free spot
m.nC = length(m.indxF); 	% number of agents with free spot
m.randI = rand(1,m.nC); 	% initialize random number vector

end