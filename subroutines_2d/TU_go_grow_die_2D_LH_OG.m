%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, inspired by Jan Poleszczuk

function [L, TUcells, TUprop] = TU_go_grow_die_2D_LH_OG( L, nh, TUcells, ...
    TUprop, TUpprol, TUpmig, TUpdeath, TUps, TUintmax)


try
    m = getAdjacent_2D_LH_OG(L,TUcells,nh); % create masks for adjacent position
catch
    warning('severe error: could not get neighborhood.');
    disp(TUcells)
    whos, disp(char(10));
    save SEVERE_ERROR; pause(60);
    warning('error log saved. will continue...');
end

% P, D and Mi are mutually exclusive; Ps and De are dependent on P
[P,D,Mi] = CellWhichAction(m.randI,TUpprol,TUpdeath,TUpmig);
Ps = P & rand(1,m.nC) <= TUps & TUprop.isStem(m.indxF) & TUprop.engaged(m.indxF) ==0; % symmetric division, can't proliferate if engaged
De = P & (TUprop.Pcap(m.indxF) == 0); % proliferation capacity exhaution -> Die
del = D | De; % find dead / dying cells
act = find((P | Mi) & ~del & TUprop.engaged(m.indxF) ==0); % live cells that will proliferate or migrate. Engaged ones can't act

newTUcells_fromProl=[];
prolTUcells=[];
n_TUcells_atprol=length(TUcells);

for iloop = 1:numel(act) % only for those that will do anything
    currID = act(iloop); % number within stack of currently acting cell
    ngh = m.S(:,m.indxF(currID)); % cells neighborhood
    ngh2 = ngh(ngh>0); % erasing places that were previously occupied
    indO = find(~L(ngh2),1,'first'); %selecting free spot
    
    if ~isempty(indO) % if there is still a free spot
        L(ngh2(indO)) = true; % add cell to grid
        if P(currID) % proliferation happens
            newCell = uint32(ngh2(indO)); % find place for new cell
            newTUcells_fromProl=[newTUcells_fromProl,newCell(1)];
            prolTUcells=[prolTUcells,TUcells(m.indxF(currID))];
            TUcells = [TUcells, newCell(1)]; % add new cell to stack
            if Ps(currID) % symmetric division
               TUprop.isStem = [TUprop.isStem, true];
               TUprop.Pcap = [TUprop.Pcap, TUprop.Pcap(m.indxF(currID))];  
               TUprop.engaged = [TUprop.engaged, 0]; % LH update property vector for engagement
               TUprop.Engcap = [TUprop.Engcap,TUintmax]; % 
            else % asymmetric division
               TUprop.isStem = [TUprop.isStem, false];
               TUprop.Pcap = [TUprop.Pcap, TUprop.Pcap(m.indxF(currID))-1];
               TUprop.engaged = [TUprop.engaged, 0]; % LH  update property vector for engagement
               TUprop.Engcap = [TUprop.Engcap,TUintmax]; % 
               if ~TUprop.isStem(m.indxF(currID)) % reduce proliferation capacity
                TUprop.Pcap(m.indxF(currID)) = TUprop.Pcap(m.indxF(currID))-1;
               end
            end
        else % migration
            L(TUcells(m.indxF(currID))) = false; % freeing spot
            TUcells(m.indxF(currID)) = uint32(ngh2(indO)); % update cell position
        end
    end
end

if ~isempty(del) % remove dead tumor cells
    L(TUcells(m.indxF(del))) = false;      % first, remove from grid
    [TUcells,TUprop] = removeTU_LH_OG(TUcells,TUprop,m.indxF(del)); % second, remove from stack
end 

TUprop.newTU=newTUcells_fromProl;
TUprop.prolTU=prolTUcells;
TUprop.n_TUcells_atprol=n_TUcells_atprol;
end
