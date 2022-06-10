%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, includes code snippets by Jan Poleszczuk


function [L, IMcells, IMprop] =  IM_go_grow_die_2D_LH_OG(IMcells, IMprop, IMpprol, IMpmig, ...
        IMpdeath, IMrwalk, IMkmax,IMpmax,engagementDuration, ChtaxMap, L, nh)
    
m = getAdjacent_2D_LH_OG(L,IMcells,nh); % create masks for adjacent positions

% P, D and Mi are mutually exclusive; Ps and De are dependent on P
[P,D,Mi] = CellWhichAction(m.randI,IMpprol,IMpdeath,IMpmig);
Pa = P & (IMprop.Kcap(m.indxF)<IMkmax) & (IMprop.engaged(m.indxF)<engagementDuration/4); %  Only the ones that have killed may proliferate: further Ps replaced with Pa. And only the ones not engaged or nearly finished with engagement can proliferate
Mis = Mi & (IMprop.engaged(m.indxF)==0); % not allowed to move if engaged
De = Pa & (IMprop.Pcap(m.indxF) == 0); % proliferation capacity exhaution -> Die. LH: If they try and proliferate beyond their allowance, they will die
del =  D | De;% % cells to delete 
act = find((Pa | Mis) & ~del); % indices to the cells that will perform action

newIMcells_fromProl=[];
prolIMcells=[];
n_IMcells_atprol=length(IMcells);

for iloop = 1:numel(act) % only for those that will do anything
    currID = act(iloop); % number within stack of currently acting cell
    ngh = m.S(:,m.indxF(currID)); % cells neighborhood
    ngh2 = ngh(ngh>0); % erasing places that were previously occupied
    indOL = find(~L(ngh2)); %selecting all free spots  
    chemo = ChtaxMap(ngh2(:)); % extract chemotaxis value at neighbors
    chemo = chemo(~L(ngh2)); % block occupied spots   
    if ~isempty(chemo) % use spot with highest chemo value
        chemo = chemo/max(chemo(:)); % normalize
        chemo = (1-IMrwalk) * chemo + IMrwalk * rand(size(chemo));
        [~,cid] = min(chemo); % lowest distance 
        indO = indOL(cid(1));   
        if ~isempty(indO) %if there is still a free spot
            L(ngh2(indO)) = true; % add new cell to grid
            if Pa(currID) % proliferation
                newIMcells_fromProl=[newIMcells_fromProl,uint32(ngh2(indO))];
                prolIMcells=[prolIMcells,IMcells(m.indxF(currID))];
                IMcells = [IMcells uint32(ngh2(indO))]; % add new cell to stack
                IMprop.Pcap(m.indxF(currID)) = IMprop.Pcap(m.indxF(currID))-1; % decrease remaining prol. cap.
                IMprop.Pcap = [IMprop.Pcap, IMprop.Pcap(m.indxF(currID))]; % update property vector for Pmax 
                IMprop.Kcap = [IMprop.Kcap, IMkmax]; % update property vector for remaining kills
                IMprop.engaged = [IMprop.engaged, 0]; % update property vector for engagement
            else % migration
                L(IMcells(m.indxF(currID))) = false; %freeing spot
                IMcells(m.indxF(currID)) = uint32(ngh2(indO));
            end
        end
    end
end

if ~isempty(del) % updating immune cell death
    L(IMcells(m.indxF(del))) = false;     % remove immune cell from grid
    [IMcells,IMprop] = removeIM_LH_OG(IMcells,IMprop,m.indxF(del)); % second, remove from stack
end
 
IMprop.newIM=newIMcells_fromProl;
IMprop.prolIM=prolIMcells;
IMprop.n_IMcells_atprol=n_IMcells_atprol;
end
