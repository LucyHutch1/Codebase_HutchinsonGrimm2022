%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, includes code snippets by Jan Poleszczuk

function [L, IMcells] =  IM_go_2D_LH_OG(IMcells,IMprop, IMpmig, IMrwalk, ChtaxMap, L, nh)

m = getAdjacent_2D_LH_OG(L,IMcells,nh); % create masks for adjacent positions

Mi = (m.randI <= IMpmig) & IMprop.engaged(m.indxF)==0; % indices to cells that will migrate, Can't migrate if engaged
act = find(Mi); % indices to the cells that will perform action

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
        [~,cid] = min(chemo); % lowest distance  (lowest number)
        indO = indOL(cid(1));   
        if ~isempty(indO) %if there is still a free spot
            L(ngh2(indO)) = true; % add new cell to grid
            L(IMcells(m.indxF(currID))) = false; %freeing spot
            IMcells(m.indxF(currID)) = uint32(ngh2(indO));
        end
    end
end

 
end
