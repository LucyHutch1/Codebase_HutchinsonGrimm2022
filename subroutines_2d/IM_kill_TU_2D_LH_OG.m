%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, includes code snippets by Jan Poleszczuk
function [TUcells, TUprop, IMcells, IMprop, L, Lt] = ...
IM_kill_TU_2D_LH_OG(TUcells, TUprop, IMcells, IMprop, L, Lt,IMpkill,nh,ChtaxMap,engagementDuration)

% Adaptations from original script: 
% - Engagement time and engagement integral are both stored. These
% determine when a tumour cell dies or when a T cell becomes apoptotic

% pre-select immune cells that may be close enough to the tumor
candidates = ChtaxMap(IMcells)<=1;
if sum(candidates(:)) % if there are candidates
    % select cells that are going to kill
    K = candidates & (IMprop.engaged==0) & (IMprop.Kcap>0) & (rand(1,length(IMcells))<IMpkill);
    actK = find(K); % cell indices
    if ~isempty(actK) % if there is a cell that is going to kill
    targetIDs = int32(zeros(1,0)); % preallocate
    killerIDs = int32(zeros(1,0)); % preallocate
    % start tumor cell killing, same random order as before
    St = bsxfun(@plus,IMcells(actK),nh.aux(nh.Pms(:,randi(nh.nP,1,length(actK)))));
    % Quick fix periodic BCs (not perfect) 
    St(St<0) = size(L,1)*size(L,2) + St(St<0);
    St(St>size(L,1)*size(L,2)) = -( size(L,1)*size(L,2)- St(St>size(L,1)*size(L,2)));
    % iterate through all immune cells and look at their neighborhood
    for jj = 1:size(St,2) 
        
        neighboorGridIndices = St(randperm(length(nh.aux)),jj);
        areNeighboorGridIndicesInTumorCells = ismember(neighboorGridIndices(:),TUcells(:));
        % Find the proposed tumour cell targets
        neighboorGridIndicesWithTumor = neighboorGridIndices(areNeighboorGridIndicesInTumorCells);
        areTumorCellsInNeighborhood = ismember(TUcells(:),neighboorGridIndicesWithTumor);
        indicesOfTumorCellsInNeighboorhood = find(areTumorCellsInNeighborhood);
        % Remove from instakill if they are already engaged
        tumorCellsInNeighboorHood = TUcells(indicesOfTumorCellsInNeighboorhood);
        areTumorCellsInNeighboorhoodEngaged = TUprop.engaged(indicesOfTumorCellsInNeighboorhood)>0;
        notEngagedTumorCells = tumorCellsInNeighboorHood(~areTumorCellsInNeighboorhoodEngaged);
        
        % if the cell encounters another cell to kill
        if sum(notEngagedTumorCells>0) > 0
            
%            immuneAttacksCounter = immuneAttacksCounter + 1;
            % if more than 1 possible targets then use the first one
            possibleTargets = notEngagedTumorCells; % possible targets
            killwhat = int32(possibleTargets(1)); % kill only the first candidate
            
            assert(ismember(killwhat, TUcells));
            
            targetIDs = [targetIDs, killwhat]; % add target ID to stack
            killerIDs = [killerIDs, IMcells(actK(jj))]; % add killer ID to stack

        end
    end

    % find indices to killed cell and killer cell. If the unlikely case
    % happens that one tumor cell is simultaneously killed by two immune cells,
    % then both will be exhausted
    auxKillTU = ismember(TUcells,targetIDs); % which tumor cells are killed
    auxKillIM = ismember(IMcells,killerIDs); % which immmune cells do kill

    if sum(auxKillTU)>0                 % if killing happens, then update  

        TUprop.engaged(auxKillTU)= engagementDuration; %  start engagement count for tumour cells
        TUprop.Engcap(auxKillTU)= TUprop.Engcap(auxKillTU)-1; % count number of interactions 
        IMprop.Kcap(auxKillIM) = IMprop.Kcap(auxKillIM)-1; % exhaust killers
        IMprop.engaged(auxKillIM) = engagementDuration; % killers are engaged
    end

    
    end % end actual killing filter
end % end candidate filter
end % end function