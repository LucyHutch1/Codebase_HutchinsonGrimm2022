%% Edited by Oliver Grimm and Lucy Hutchinson 2021
%
% Original: 2016-2017, created by JN Kather, includes code snippets by Jan Poleszczuk

function sysOut = runSimulation(MySystem,randomSeed,par_all,par_names,saveOutput,saveImage,saveVideo,timepoint,vary_params)    
format compact %, clc - clc clears the command window
addpath('./subroutines_2d/'); % include all functions for the 2D model
addpath('./subroutines_nd/'); % include generic subroutines


[tileDim_x, tileDim_y] = size(MySystem.grid.L);
[sysTempl, cnst] = getSystemParams_LH_OG('2D',[tileDim_x tileDim_y]); 
sysTempl.grid = MySystem.grid; 
sysTempl.TU = MySystem.TU;
sysTempl.IM = MySystem.IM;

% overwrite parameters for optimization
if vary_params == 1
    for i = 1:length(par_all)
        eval(['sysTempl.params.',par_names{i},'=',num2str(par_all(i))]);
    end
    sysTempl.IM.IMprop.Pcap(sysTempl.IM.IMprop.Pcap >= sysTempl.params.IMpmax) = sysTempl.params.IMpmax;
end
    
% overwrite parameters for sensitivity analysis
if vary_params == 2
    eval(['sysTempl.params.',par_names,'=',num2str(par_all)]);
end

cnst.createNewSystem = false; %set to true if not using input tile

cnst.nSteps   = timepoint*2;  % how many iterations in the first place (1 time step is 12 hours) 
cnst.drawWhen = 2;    % draw after ... iterations DEFAULT 2
cnst.verbose = false;
%cnst.drawsteps = true;


% Change random seed for each stochastic run
sysTempl.params.initialSeed = randomSeed;

%% Update output file name

if vary_params
    [sysOut, ~, ~, ~] = growTumor_2D_LH_OG(sysTempl,cnst);
else
    [sysOut, lastFrame, summary, fcount] = growTumor_2D_LH_OG(sysTempl,cnst);
end

format compact
end % Function