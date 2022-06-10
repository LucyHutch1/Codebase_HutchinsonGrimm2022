%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, includes code snippets by Jan Poleszczuk

function [IMcells,IMprop] = removeIM_LH_OG(IMcells,IMprop,idx)

    IMcells(idx) = [];           % remove from stack
    IMprop.Pcap(idx) = [];     % remove Pmax property
    IMprop.Kcap(idx) = [];       % remove Kmax property
    IMprop.engaged(idx) = [];       % remove engagement property
    
end

