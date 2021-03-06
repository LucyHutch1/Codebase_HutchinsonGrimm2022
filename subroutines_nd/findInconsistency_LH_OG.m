%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, includes code snippets by Jan Poleszczuk

function out = findInconsistency_LH_OG(Lt,Lf,Ln,TUcells,IMcells,TUprop,IMprop)   
disp('checking for inconsistency...');

	% perform basic consistency checks...
   c(1) = numel(TUcells) == sum(Lt(:));
   c(2) = numel(TUcells) == numel(TUprop.Pcap);
   c(3) = numel(TUcells) == numel(TUprop.isStem); 
   c(4) = numel(IMcells) == numel(IMprop.Pcap); 
   c(5) = numel(IMcells) == numel(IMprop.Kcap); 
   c(6) = ~sum(Ln(TUcells));
   c(7) = ~sum(Lt(IMcells));
   c(8) = ~sum(Lf(TUcells)); 
   c(9) = numel(TUcells) == numel(TUprop.engaged);
   c(10) = numel(TUcells)== numel(TUprop.Engcap);
   
   out = any(~c);
   
   if out % error occured
       disp(c);
   else
	   disp('no problem found');
   end
end