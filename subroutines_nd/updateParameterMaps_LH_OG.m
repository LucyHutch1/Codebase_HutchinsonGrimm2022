%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, includes code snippets by Jan Poleszczuk


function [ChtaxMap, HypoxMap] = updateParameterMaps_LH_OG(Lt,Ln,Lf,fillSE,distMaxNecr)
% update chemotaxis map, has to be double for immune cell migration to 
% work locally

dimen = size(Lt,1);
% generate matrix of Lt and 8 Lts around it
Lt_ext = repmat(Lt,3,3);
%ChtaxMap = double(bwdist(Lt,'euclidean'));
ChtaxMap_ext = double(bwdist(Lt_ext,'chessboard')); %

ChtaxMap = ChtaxMap_ext(dimen+1:2*dimen, dimen+1:2*dimen); %crop for the central copy of biopsy grid

% update hypoxia map (same method for dealing with periodics )
%{ 
% commented out hypoxia 
tem = ~imdilate(Lt|Ln|Lf,fillSE);
tem_ext = repmat(tem,3,3);
HypoxMap_ext = min(double(bwdist(tem_ext,'euclidean')),distMaxNecr) / distMaxNecr;
HypoxMap = HypoxMap_ext(dimen+1:2*dimen, dimen+1:2*dimen);
%}
%HypoxMap = min(double(bwdist(~imdilate(Lt|Ln|Lf,fillSE),'euclidean')),distMaxNecr) / distMaxNecr;

% HypoxMap is essentially a distance map for the distance from the tumor edge
% within the tumor. All values larger than distMaxNecr are cut at this
% threshold. Then, the mask is normalized to 0 ... 1 with 1 being distMaxNecr

HypoxMap=false(dimen);
end