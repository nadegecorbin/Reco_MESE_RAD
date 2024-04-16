% config_file 

function write_param(output_folder,step)
%% Where will processed data be written ? 
Param.output_folder=output_folder;
%% Where can we find the toolboxes? 
Param.path_to_spm='/opt/code/docker_save/SPM/SPM_r8091'
Param.path_to_EPGX='/workspace_QMRI/PROJECTS_DATA/2021_RECH_NC_MESE_RADIAL/Reco_MESE_RAD/EPG-X'
Param.path_to_sdc3='/workspace_QMRI/PROJECTS_DATA/2021_RECH_NC_MESE_RADIAL/Reco_MESE_RAD/sdc3_nrz_11aug'

%% Which options of reconstruction? 

switch step
    case 0

Param.DoEchoes=1; 
Param.DoTRImages=0;
Param.Moco=0;

    case 1
% reconstruct echoes ? 
Param.DoEchoes=0;

% reconstruct intermediate images ? 
Param.DoTRImages=1;


% perform motion correction? 
Param.Moco=0;

    case 2

Param.DoEchoes=1; 
Param.DoTRImages=0;
Param.Moco=1;
end


% save in the current folder 
save Param Param
end