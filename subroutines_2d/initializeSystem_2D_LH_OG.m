%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, includes code snippets by Jan Poleszczuk

function [L, TUcells, IMcells, TUprop, IMprop] = initializeSystem_2D_LH_OG(N,M,TUpmax, TUintmax)
    
    % START INITIALIZE GRID  -------------------------------------------
    L = false(N,M);
    % END INITIALIZE GRID  ----------------------------------------------
    
    % START INITIALIZE TUMOR CELLS -------------------------------------------
    TUcells = int32([250]); % % first TU cell is centered
    L(TUcells) = true; 			 % place first tumor cell on the grid
    TUprop.isStem = [true];        % set property of first cell: stemness
    TUprop.Pcap = uint8([TUpmax]); % set property of first cell: proliferation capacity
    TUprop.Engcap = uint8([TUintmax]);          %  set property: tumour cells must die after 2 interactions
    TUprop.engaged = uint8([0]); % add properties: engagement in killing 
    % END INITIALIZE TUMOR CELLS ---------------------------------------------

    % START INITIALIZE IMMUNE CELLS (LYMPHOCYTES) -----------------------------
    IMcells = int32([5000]); 	 % preallocate immune cell position vector
    IMprop.Pcap = uint8([10]); % add properties: max proliferation capacity
    IMprop.Kcap = uint8([10]); % add properties: max killing capacity
    IMprop.engaged = uint8([40]); % add properties: engagement in killing
     % END INITIALIZE IMMUNE CELLS --------------------------------------------

end