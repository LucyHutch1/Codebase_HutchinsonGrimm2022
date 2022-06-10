%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, includes code snippets by Jan Poleszczuk

function [mySystem,currentSummary] = updateSystem_LH_OG(mySystem,TUcells,TUprop,...
    IMcells,IMprop,ChtaxMap,HypoxMap,L,Lt,Ln,Lf,i,cnst,tucells_killed,IMinfluxcoords)

    % copy all variables back to mySystem
    mySystem.TU.TUcells = TUcells;
    mySystem.TU.TUprop.isStem = TUprop.isStem;
    mySystem.TU.TUprop.Pcap = TUprop.Pcap;
    mySystem.TU.TUprop.Engcap = TUprop.Engcap; % 
    mySystem.TU.TUprop.engaged = TUprop.engaged; %

    
    mySystem.IM.IMcells = IMcells;
    mySystem.IM.IMprop.Kcap = IMprop.Kcap;
    mySystem.IM.IMprop.Pcap = IMprop.Pcap;  
    mySystem.IM.IMprop.engaged = IMprop.engaged;  
    mySystem.IM.IMinfcoords = IMinfluxcoords;
    
    mySystem.grid.ChtaxMap = ChtaxMap;
    mySystem.grid.HypoxMap = HypoxMap;
    mySystem.grid.Ln = Ln;
    mySystem.grid.Lf = Lf;
    mySystem.grid.L = L;
    mySystem.grid.Lt = Lt;
    mySystem.grid.StepsDone = i;
    
    % create immune grid
    mySystem.grid.Li = false(size(L));
    mySystem.grid.Li(IMcells) = true;
    
    % summarize system
    if ndims(L) == 2
        currentSummary = summarizeSystem_2D(mySystem,cnst);
    elseif ndims(L) == 3
        error('not yet implemented');
    end
    
end
