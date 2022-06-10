%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, includes code snippets by Jan Poleszczuk

function [IMcells,IMprop] = shuffleIM_LH_OG(IMcells,IMprop)

shf = randperm(length(IMcells)); % prepare random shuffling
IMcells = IMcells(shf); % randomly shuffle cells
IMprop.Pcap = IMprop.Pcap(shf); % shuffle Pmax property accordingly
IMprop.Kcap = IMprop.Kcap(shf); % shuffle Kmax property accordingly
IMprop.engaged = IMprop.engaged(shf); % shuffle engaged property accordingly

end