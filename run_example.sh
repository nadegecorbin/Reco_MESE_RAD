
echo "run example of retrospective motion correction"  

# Raw data in ISMRMRD format
raw_data=$1

#output_folder
output_folder=$2

matlab -batch "write_param('$output_folder/Nomoco',0)"
gadgetron_ismrmrd_client -p 9002 -f $raw_data -c MESE_RAD_shared_12ch.xml

echo "reconstruct intermediate images"
matlab -batch "write_param('$output_folder/Moco',1)"
gadgetron_ismrmrd_client -p 9002 -f $raw_data -c MESE_RAD_shared_6ch.xml


matlab -batch "write_param('$output_folder/Moco',2)"
gadgetron_ismrmrd_client -p 9002 -f $raw_data -c MESE_RAD_shared_12ch.xml


# For T2 map reconstruction with the Anima toolbox, uncomment those lines
#Echoes=$output_folder/Nomoco/RecoEchoes/Echoes4D.nii      
#python3 animaT2relaxometry.py -i $Echoes --mono-out $output_folder/Nomoco/RecoEchoes/Anima_T2map.nii --tr-value 1000 -e 7 --no-brain-masking
#Echoes=$output_folder/Moco/RecoEchoes/Echoes4D.nii     
#python3 animaT2relaxometry.py -i $Echoes --mono-out $output_folder/Moco/RecoEchoes/Anima_T2map.nii --tr-value 1000 -e 7 --no-brain-masking




    
      
