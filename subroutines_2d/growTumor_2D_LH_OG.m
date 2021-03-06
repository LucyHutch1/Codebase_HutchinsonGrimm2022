%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, includes code snippets by Jan Poleszczuk

function [mySystemALL, finalImage, finalSummary, fcount] = ...
    growTumor_2D_LH_OG(mySystem,cnst)
% growTumor_2D performs the actual agent-based modeling in 2D
%       input are two structures: mySystem, defining the initial state of
%       the system; and cnst, defining some global constants

% START PREPARATIONS -------------------------------------------
% throw all model parameters to workspace
cellfun(@(f) evalin('caller',[f ' = mySystem.params.' f ';']), fieldnames(mySystem.params));
cellfun(@(f) evalin('caller',[f ' = mySystem.grid.' f ';']), fieldnames(mySystem.grid));
if cnst.createNewSystem % create a new (empty) system
    [L, TUcells, IMcells, TUprop, IMprop] = initializeSystem_2D_LH_OG(N,M,TUpmax,TUintmax);
    Ln=0; % placeholders for Ln and Lf 
    Lf=0;
  %^  Ln = false(size(L));    % initialize necrosis map
  %^  Lf = false(size(L));    % initialize fibrosis map
else % use existing system and grow it further
    cellfun(@(f) evalin('caller',[f ' = mySystem.TU.' f ';']), fieldnames(mySystem.TU));
    cellfun(@(f) evalin('caller',[f ' = mySystem.IM.' f ';']), fieldnames(mySystem.IM));
end
% END PREPARATIONS -------------------------------------------
    
% START INITIALIZE AUX VARIABLES  ----------------------------------------
nh = neighborhood_2D(N);   % get neighborhood indices)
fcount = 0;             % frame counter for video export
finalImage = [];        % set empty resulting image
% END INITIALIZE AUX VARIABLES   -----------------------------------------

% START ITERATION
for i = 1:cnst.nSteps % iterate through time steps
   probSeedNecr = 0; 

% START TUMOR CELL ROUND ------------------------------------------------
%^L(Lf) = rand(sum(Lf(:)),1)>stromaPerm; % permeabilize some stroma-filled grid cells
L([IMcells,TUcells]) = true; % ensure that all cells are present on the grid
[TUcells,TUprop] = shuffleTU(TUcells,TUprop);
[L, TUcells, TUprop] = TU_go_grow_die_2D_LH_OG(L, nh, TUcells, TUprop, TUpprol, TUpmig, TUpdeath, TUps,TUintmax);
Lt = updateTumorGrid(L,TUcells); % update tumor grid
% END TUMOR CELL ROUND ---------------------------------------------------

% START MODIFY PARAMETER MAPS --------------------------------------------
[ChtaxMap, HypoxMap] = updateParameterMaps_LH_OG(Lt,Ln,Lf,fillSE,distMaxNecr);
% END MODIFY PARAMETER MAPS

% START IMMUNE CELL ROUND ------------------------------------------------
L([IMcells,TUcells]) = true; % ensure that all cells are present on the grid
IMinfluxcoords=[];
if rand()<=IMinfluxProb % randomly trigger influx
[L,IMcells,IMprop,IMinfluxcoords] = IMinflux_LH_OG(L,IMcells,IMprop,IMpmax,IMkmax,IMinflRate);
end

[IMcells,IMprop] = shuffleIM_LH_OG(IMcells,IMprop); % randomly shuffle immune cells

n_tum1 = length(TUcells);
if numel(IMcells)>0 % if there are any immune cells 
for j = 1:(IMspeed-1) % allow immune cells to move N times per round
  
    L([IMcells,TUcells]) = true; % ensure that all cells are present on the grid
    [L, IMcells] =  IM_go_2D_LH_OG(IMcells,IMprop, IMpmig, IMrwalk, ChtaxMap, L, nh);
    [TUcells, TUprop, IMcells, IMprop, L, Lt] = ... % tumor cell killing by immune cells 
    IM_kill_TU_2D_LH_OG(TUcells, TUprop, IMcells, IMprop, L, Lt,IMpkill,nh,ChtaxMap,engagementDuration);
    IMprop.engaged(IMprop.engaged>0) = IMprop.engaged(IMprop.engaged>0)-1; % un-engage lymphocytes
    TUprop.engaged(TUprop.engaged>0) = TUprop.engaged(TUprop.engaged>0)-1; % un-engage tumour cells
   
%% Remove dead cells due to engagement
 Di = TUprop.Engcap == 0 & TUprop.engaged == 0; % Get rid of killed ones
 if sum(Di)>0
    L(TUcells(Di))=false;
    Lt(TUcells(Di))=false; %  
    [TUcells,TUprop] = removeTU_LH_OG(TUcells,TUprop,Di);
 end

end

% allow immune cells to move once more or to proliferate or die

[L, IMcells, IMprop] =  IM_go_grow_die_2D_LH_OG(IMcells, IMprop, IMpprol, IMpmig, ...
        IMpdeath, IMrwalk, IMkmax,IMpmax,engagementDuration, ChtaxMap, L, nh); 
end % end (if there are any immune cells)

n_tum2 = length(TUcells);

tucells_killed = n_tum1-n_tum2;

%^L(Ln|Lf) = true; % fully turn on necrosis and fibrosis again

orig = sum(sum(L));
L([IMcells,TUcells]) = true; % ensure that all cells are present on the grid
new = sum(sum(L));
%disp(orig-new);


% END IMMUNE CELL ROUND --------------------------------------------------

% START NECROSIS  --------------------------------------------
%{
necrNum = sum(rand(numel(TUcells),1) <= probSeedNecr);
if numel(TUcells)>1 && necrNum>0 % seed necrosis
    seedCoords = randsample(TUcells,necrNum,true,HypoxMap(TUcells));
    necrosisSeeds = false(N,M);
    necrosisSeeds(seedCoords) = true; 
    %disp([num2str(necrNum), ' cell(s) will trigger necrosis']);
    % smooth and expand necrotic seed map
    necrosisSeeds = expandSeedMap(necrosisSeeds,smoothSE,necrFrac);
    seedCoords = find(necrosisSeeds);
    targetIdx = ismember(TUcells,seedCoords); % find indexes of erased tumor cells
    Lt(TUcells(targetIdx)) = false; % remove cell from grid Lt
    Ln(TUcells(targetIdx)) = true; % add to necrosis grid    
    [TUcells,TUprop] = removeTU_LH(TUcells,TUprop,targetIdx); % second, remove from stack
end
%}
% END NECROSIS  ----------------------------------------------

% START FIBROSIS ------------------------------------------------------
%{
fibrosify = IMprop.Kcap<=95 & (rand(1,numel(IMcells))<probSeedFibr); %
if sum(fibrosify) % exchausted immune cells seed fibrosis
    Lfseed = false(size(L)); % preallocate fibrosis seed map
    Lfseed(IMcells(fibrosify)) = true;
    Lfseed = expandSeedMap(Lfseed,smoothSE,fibrFrac); % smooth and expand fibrotic seed map
    Lfseed(TUcells) = false; 
    Lf(Lfseed & ~Ln) = true; % update fibrosis grid
    [IMcells,IMprop] = removeIM(IMcells,IMprop,fibrosify); % remove fibrosifying immune cells
    L(Lf) = true; % update L grid (master grid)
end
%}
% END FIBROSIS ----------------------------------------------------------

% START DRAWING AND SAVING ---------------------
tumorIsGone = (sum(Lt(:))==0);
if (mod(i-1,cnst.drawWhen)==cnst.drawWhen-1) || tumorIsGone % plot status after N epochs 
    fcount = fcount+1;
    % export current state of the system
    [mySystem,currentSummary] = updateSystem_LH_OG(mySystem,TUcells,TUprop,...
        IMcells,IMprop,ChtaxMap,HypoxMap,L,Lt,Ln,Lf,i,cnst,tucells_killed,IMinfluxcoords);
    
    finalSummary{fcount} = currentSummary; %  Moved out of verbose loop so this is always saved
    mySystemALL{fcount}=mySystem; % 
    
    if cnst.verbose % enforce drawing and create image from plot
        visualize_balls_2D_blank(mySystem);
        drawnow, currFrame = getframe(gca);
        finalImage{fcount} = currFrame.cdata;    
    else finalImage = []; finalSummary = [];
    end

end
% END DRAWING AND SAVING--------------------------

% if there are no tumor cells anymore then the game is over
if tumorIsGone
    disp('Immune cells win');
    imWin = i;
    return
end
if debugmode && findInconsistency_LH_OG(Lt,Lf,Ln,TUcells,IMcells,TUprop,IMprop)  
     error('SEVERE ERROR: INCONSISTENCY FOUND');
end % END DEBUG

end % END ITERATION
end % END FUNCTION
        