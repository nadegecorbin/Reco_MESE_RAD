# Example code to reconstruct and retrospectively correct for motion. 

To be used with data acquired with the 3D radial Multi-Echo Spin echo sequence.

One example dataset with and without deliberate motion is provided in ISMRMRD format and available at: 
10.5281/zenodo.10908742. Reconstructed maps are also provided. 

The example code can be used to reconstruct the motion corrupted dataset with and without motion correction.


The related paper "Whole-brain T2 mapping with radial sampling and retrospective motion correction at 3T"
is currently under submission. 

*Nadege Corbin, 
Centre de Resonance Magnetique des Systemes Biologiques 
UMR 5536 CNRS/Bordeaux University 
2024.04.16*



# Requirements


This reconstruction code requires multiple softwares and toolbox

- MATLAB (tested with R2019b) 
- Gadgetron (https://github.com/gadgetron/gadgetron)
- SPM (https://github.com/spm/)
- EPG-X (https://github.com/mriphysics/EPG-X)
- Anima (https://github.com/Inria-Empenn/Anima-Public)
- sdc3_nrz_11aug (https://github.com/ISMRM/mri_unbound/tree/master/sdc3) 



# Step 1: copy files in the right place


- create the folder `[path to gadgetron]/+gadgetron/+custom/+MESE_RAD`
- copy `handle_connection_MESE_RAD_shared.m` there


- if it does not exist, create the folder ```
[path to gadgetron]/+gadgetron/+custom/+utils ```
- copy `apply_traj_kdata_clean.m` there
- copy `unif_traj.m` there
- copy `write_basis_shared.m` there 


- copy `MESE_RAD_shared_6ch.xml  MESE_RAD_shared_12ch.xml`  into `/usr/local/share/gadgetron/config`



# Step 2: Update path 

update the path to SPM, EPG-X and sdc3_nrz_11aug in the matlab file `write_param.m`


# Step 3: start gadgetron

Open a terminal 
go to the "Reco_MESE_RAD" folder
start gadgetron with the  command: `gadgetron -p 9002`

Open a second terminal:
go to the "Reco_MESE_RAD" folder
run example script with the following command : `./run_example.sh [h5_file_name] [output_folder]`





