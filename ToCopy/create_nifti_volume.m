%%*******************************************************%%
%%% Create Nifti file
%%%
%%% Create nifti file with the appropriate header to match Siemens
%%% reconstruction
%%%     - data : 3D volume to save in nifti format
%%%     - file_name: file name of the output nifti file
%%%     - header: mini header
%%%     - res: resolution of the data
%%%     - matrix_size : matrix dimension of the data
%%%
%%% Note : Only available for Sagittal and transversal and 
%%%
%%% NC 08.04.20
%%%********************************************************



function create_nifti_volume(data,file_name,header,res,matrix_size)

addpath('/opt/code/docker_save/SPM/spm12')

% Initialize nifti file
N=nifti;
dat=file_array;
dat.dtype='FLOAT64';
dat.fname='dummy.img';
N.dat=dat;

% update file name
N.dat.fname=file_name;


%% Fill mat field 
% for now, I have no solution to fill the header with a global method that
% would work with any orientation 
% Besides, I don't have a way to select which orientation is actually used 
% The parameter "orientation" is then hard-coded for now but more work is required to
% make the code more general

orientation='Tra';


if strcmp(orientation,'Sag')==1

disp(" The orientation of the FOV is assumed to be sagittal")
%%  Sagittal 
% Fill mat file
MatTransf=[res(3) 0 0;0 0 -res(1);0 res(2) 0];
Mrot=[header.read_dir', header.phase_dir', header.slice_dir'];

Mat33=MatTransf*Mrot;


N.mat(1:3,1:3)=Mat33;
N.mat(1,4)=-((matrix_size(3)/2+0.5)*res(3)+header.position(1));
N.mat(2,4)=(matrix_size(2)/2+1)*res(2)-header.position(2);
N.mat(3,4)=-(matrix_size(1)/2*res(1)-header.position(3));
N.mat(4,4)=1;


% transform data to match Siemens output
% permute
data=permute(data,[2 1 3]);
%flip one dimension
data=data(end:-1:1,end:-1:1,:);
% shift from 1 pixel the second dimension
data=circshift(data,-1,2);

% set dimension of matrix
N.dat.dim=size(data);

% fill dat field
N.dat(:,:,:)=data;

end

if strcmp(orientation,'Tra')==1
%% Transversal
% fill mat file
%% Transversal
% fill mat file
MatTransf=[res(1) 0 0;0 res(2) 0 ;0 0 res(3)];
Mrot=[header.read_dir', header.phase_dir', header.slice_dir'];

Mat33=MatTransf*Mrot;

N.mat(1:3,1:3)=Mat33;

N.mat(1,4)=(matrix_size(2)/2+1)*res(1)-header.position(1);
N.mat(2,4)=-(matrix_size(1)/2*res(1)+header.position(2));
N.mat(3,4)=-((matrix_size(3)/2+0.5)*res(3)-header.position(3));
N.mat(4,4)=1;

% transform data to match SIemens outputedit 
% permute
% flip one dimension
data=permute(data,[1 2 3]);
data=data(:,:,end:-1:1);
% shift from one voxel the second dimension
data=circshift(data,-1,2);



% set dimension of matrix
N.dat.dim=size(data);
% fill dat field
N.dat(:,:,:)=data;
end

% create nifti file
create(N);
end