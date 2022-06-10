%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, inspired by Jan Poleszczuk

function [TUcells,TUprop] = removeTU_LH_OG(TUcells,TUprop,idx)

    TUcells(idx) = [];           % remove from stack
    TUprop.isStem(idx) = [];     % remove stemness property
    TUprop.Pcap(idx) = [];       % remove Pmax property
    TUprop.engaged(idx) = [];    % remove engaged property    
    TUprop.Engcap(idx) = [];    % remove engaged property
end