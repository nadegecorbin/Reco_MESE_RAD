
function handle_connection_MESE_RAD_shared(connection)

load('Param.mat')
output_folder=Param.output_folder;

Moco=Param.Moco;
DoTRImages=Param.DoTRImages;
DoEchoes=Param.DoEchoes;

setenv('OMP_NUM_THREADS', '10')

addpath(genpath(Param.path_to_spm))
addpath(genpath(Param.path_to_sdc3))
addpath(genpath(Param.path_to_EPGX))

%% Get Motion file to apply and regulCmd
disp("handle_connection was called.")

next_acquisition = @connection.next;


sRO = connection.header.encoding.encodedSpace.matrixSize.x*2;
sPE = connection.header.encoding.encodedSpace.matrixSize.y;
sSE = connection.header.encoding.encodedSpace.matrixSize.z;


nContrast = connection.header.encoding.encodingLimits.contrast.maximum+1;
nSet = connection.header.encoding.encodingLimits.set.maximum+1;
nRad = connection.header.encoding.encodingLimits.kspace_encoding_step_1.maximum+1;


matrix_size = [connection.header.encoding.reconSpace.matrixSize.x...
    connection.header.encoding.reconSpace.matrixSize.x...
    connection.header.encoding.reconSpace.matrixSize.x
    ];

res = [connection.header.encoding.encodedSpace.fieldOfView_mm.x/connection.header.encoding.encodedSpace.matrixSize.x ...
    connection.header.encoding.encodedSpace.fieldOfView_mm.x/connection.header.encoding.encodedSpace.matrixSize.x ...
    connection.header.encoding.encodedSpace.fieldOfView_mm.x/connection.header.encoding.encodedSpace.matrixSize.x]; % z is wrong !

acquisition = next_acquisition(); % Call input function to produce the next acquisition.

nCh = size(acquisition.data.data,2);

RefCount=acquisition.ref.count;
Set_idx = acquisition.data.header.set(RefCount+1:end) + 1;

rad_idx=acquisition.data.header.kspace_encode_step_1(RefCount+1:end)+1;


TE_idx = acquisition.data.header.contrast(RefCount+1:end) + 1;

% store data in a list
kdata=reshape(acquisition.data.data(:,:,RefCount+1:end),[sRO,nCh,nContrast,nRad]);


%% prepare parameter for reconstruction
kdata = permute(kdata,[5 1 4 2 3]);
kdata=reshape(kdata,[1 sRO nRad nCh 1 nContrast]);


%% generate traj

traj_bart = gadgetron.custom.utils.unif_traj(sRO,nContrast*nRad*2);

traj_bart=traj_bart(:,:,1:size(traj_bart,3)/2);


%% measure dcf ->


nTR=connection.header.encoding.trajectoryDescription.userParameterLong(3).value;

%% traj_bart per Echo
clear traj_bart_allec
for ec=1:nContrast

    traj_bart_ec=traj_bart(:,:,rad_idx(ec:nContrast:end));

    traj_bart_allec(:,:,:,1,1,ec)=traj_bart_ec;
end


%% --------------------------------------------
% Intermediate image recostruction 
%----------------------------------------------

if DoTRImages

    %% measure dcf ->

    doDCF_tr=1;
    PICS_tr=1;
    nT=nRad/nTR;
    oldsRO=sRO;
    sRO=sRO/2;
    kdata=kdata(:,oldsRO/2+[-(sRO/2-1):sRO/2],:,:,:);



    %% sensitivity coil map with low res

    for ec=1
        traj_bart_ec=traj_bart_allec(:,oldsRO/2+[-(sRO/2-1):sRO/2],:,1,1,ec);
        im_nufft = bart('nufft -i ', traj_bart_ec,kdata(:,:,:,:,1));
        kspace= bart('fft 7',im_nufft);

        % func = sprintf('caldir 24');
        dimref=[24 24 24];
        func = sprintf('ecalib -r %d %d %d -m 1 -d 3 -c0', dimref(1), dimref(2), dimref(3));
        sensitivity_coil_map = bart(func, kspace);
        disp('---- Sensitivity map estimated ----')
    end



    for tr=1:nRad/nTR/nT
        clear DCF_bart
        clear traj_bart_alltr
        clear kdata_alltr
        for t=1:nT
            temp= traj_bart_allec(:,:,(tr-1)*nTR*nT+(t-1)*nTR+[1:nTR],1,1,:);
            temp=permute(temp,[1 2 6 3 4 5 ] );
            temp=reshape(temp,[3,oldsRO,nContrast*nTR 1 1]);
            traj_bart_alltr(:,:,:,1,1,t)=temp(:,oldsRO/2+[-(sRO/2-1):sRO/2],:,1,1);
            kdata_tr=kdata(:,:,(tr-1)*nTR*nT +(t-1)*nTR+[1:nTR],:,:);
            kdata_tr=permute(kdata_tr,[1 2 5 3 4]);
            kdata_alltr(:,:,:,:,1,t)=reshape(kdata_tr,[1 sRO nContrast*nTR nCh]);

        end
        for t=1:nT

            % reshape DCF to 0.5
            traj_bart_tr= traj_bart_alltr(:,:,:,1,1,t);
            traj_dcf = traj_bart_tr/(2*max(abs(traj_bart_tr(:))));

            % reshape to 3,proj*RO
            traj_dcf = reshape(traj_dcf,3,[]);

            disp('   SDC params:');
            numIter = 25;
            effMtx  = size(kdata_tr,2);
            osf     = 2;
            verbose = 1;

            disp('   start SDC calc')
            DCF_tr = sdc3_MAT(traj_dcf,numIter,effMtx,verbose,osf);
            DCF_tr= reshape(DCF_tr,1,[],nContrast*nTR);
            DCF_bart(:,:,:,1,1,t) = DCF_tr;

        end


        DCF_bart = repmat(DCF_bart,[1 1 1 nCh 1 1]);
        writecfl('DCF_bart',sqrt(DCF_bart));
        writecfl('kdata',kdata_alltr);
        writecfl('sens_coil',sensitivity_coil_map);




        %% reconstruct with PICS

        writecfl('traj_bart',traj_bart_alltr);
        regulCmd= 'bart pics --lowmem -i50 -e -d5 -S -R G:32:0:1 -t traj_bart -p DCF_bart kdata sens_coil im_tr';


        system(regulCmd);
        im_tr(:,:,:,(tr-1)*nT+[1:nT])=squeeze(readcfl('im_tr'));

        % save info in struct
        info_tr.nTR=nTR;
        info_tr.reco=regulCmd;
        info_tr.PICS=PICS_tr;
        info_tr.DCF=doDCF_tr;
        info_tr.nT=nT;
        info_tr.Moco=Moco;

        info_tr.comment='';
        if Moco==1
            info_tr.Moco=rp_file;
        end

    end



    % Update Output folder
    full_output_folder= [output_folder '/RecoTRImages'];
    disp(['Nifti files will be created and saved in : ' full_output_folder])
    mkdir(full_output_folder)

    save(fullfile(full_output_folder,'info_tr.mat'),'info_tr');
    im=im_tr;
    Name='TRImages';
    nifticorrect=0;
    for k=1:size(im,4)

        % create nifti file
        nifti_file_name= [full_output_folder '/' Name '_' sprintf('%.4d',k) '.nii'];
        sRO = connection.header.encoding.encodedSpace.matrixSize.x*2;

        createNifti(abs(im(:,:,:,k)),nifti_file_name);

        disp("create_and_save_nifti setup...")
    end


    %% estimate motion


    % Create 4D data of TR images
    ListFiles=cellstr(spm_select('FPList',full_output_folder,'^TR'));
    Output='Orig_TRImages'
    create4Ddata(ListFiles,Output)

    % Realign data
    matlabbatch{1}.spm.spatial.realign.estimate.data = {ListFiles};
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estimate.eoptions.weight = '';
    spm_jobman('run', matlabbatch);
    clear matlabbatch

end

%% --------------------------------------------
% Motion correction
%----------------------------------------------

if Moco==1

    rp_file_folder= [output_folder '/RecoTRImages'];
    rp_file=fullfile(rp_file_folder,'rp_TRImages_0001.txt')
    fromrpfile=load(fullfile(rp_file_folder,'info_tr.mat'))
    nTR=fromrpfile.info_tr.nTR;
    ratio=2 % voxel twice larger in the TRImages
    [kdata traj_bart_allec]=gadgetron.custom.utils.correct_traj_kdata_clean(kdata, traj_bart_allec, rp_file, nRad, nTR, nContrast, sRO,ratio) ;
end




%% --------------------------------------------
% Echo reconstruction 
%----------------------------------------------
if DoEchoes==1

    clear DCFr
    clear DCF
    for ec=1:nContrast

        traj_bart_ec=traj_bart_allec(:,:,:,1,1,ec);

        % reshape DCF to 0.5

        traj_dcf = traj_bart_ec/(2*max(abs(traj_bart_ec(:))));

        % reshape to 3,proj*RO
        traj_dcf = reshape(traj_dcf,3,[]);

        disp('   SDC params:')
        numIter = 25
        effMtx  = size(kdata,2);
        osf     = 2
        verbose = 1

        disp('   start SDC calc')
        DCF(:,ec)= sdc3_MAT(traj_dcf,numIter,effMtx,verbose,osf);
        DCFr(:,:,ec)=repmat(DCF(:,ec),1,nCh);


    end

    % sensitivity coil map

    ec=1
    traj_bart_ec=traj_bart_allec(:,:,:,1,1,ec);
    im_nufft = bart('nufft -a',traj_bart_ec,kdata(:,:,:,:,:,ec).*reshape(DCFr(:,:,ec),size(kdata(:,:,:,:,:,ec))));
    kspace= bart('fft 7',im_nufft);


    dimref=[24 24 24];
    func = sprintf('ecalib -r %d %d %d -m 1 -d 3 -c0.1', dimref(1), dimref(2), dimref(3));
    sensitivity_coil_map = bart(func, kspace);
    disp('---- Sensitivity map estimated ----')

    writecfl('sensitivity_coil_map',sensitivity_coil_map)



    DCF_bart = reshape(DCF,1,[],nRad,1,1,nContrast);

    %basis file
    TE=connection.header.sequenceParameters.TE;
    nT2=200;
    nB1=40;
    filename='basisfile';
    nComp=12;
    compl=1;
    loga=1;
    gadgetron.custom.utils.write_basis_shared(TE,nT2,filename,nComp,nB1,compl,loga);



    writecfl('traj_bart_ec',traj_bart_allec);
    writecfl('DCF_bart',sqrt(DCF_bart))
    writecfl('kdata',kdata)
    writecfl('sensitivity_coil_map',sensitivity_coil_map)


    regulCmd='bart pics --lowmem -i50 -e -d5 -S -R W:7:0:0.0002 -B basisfile -t traj_bart_ec -p DCF_bart kdata sensitivity_coil_map output';


    system(regulCmd);
    im_ec=squeeze((readcfl('output')));


    % back to original space
    basis_final=readcfl('basisfile');
    coeff_maps=readcfl('output');
    im_ec=squeeze(bart('fmac -s 64', basis_final, coeff_maps));


    delete('output.hdr')
    delete('output.cfl')
    delete('kdata.hdr')
    delete('kdata.cfl')
    delete('sensitivity_coil_map.hdr')
    delete('sensitivity_coil_map.cfl')
    delete('traj_bart_ec.hdr')
    delete('traj_bart_ec.cfl')

    % save info in struct
    info_ec.TE=connection.header.sequenceParameters.TE;
    info_ec.TR=connection.header.sequenceParameters.TR;
    info_ec.FA=connection.header.sequenceParameters.flipAngle_deg;
    info_ec.reco=regulCmd;
    info_ec.Moco=Moco;

    info_ec.comment='';

    info_ec.senscoil=func;
    info_ec.nCh=nCh;

    if Moco==1
        info_ec.Moco=rp_file;
    end



    info_ec.basisfile.nComp=nComp;
    info_ec.basisfile.compl=compl;
    info_ec.basisfile.loga=loga;
    info_ec.basisfile.nT2=nT2;
    info_ec.basisfile.nB1=nB1;

end


% Update Output folder
full_output_folder= [output_folder  '/RecoEchoes'];
disp(['Nifti files will be created and saved in : ' full_output_folder])
mkdir(full_output_folder)
disp(['create' full_output_folder])

save(fullfile(full_output_folder,'info_ec.mat'),'info_ec');
im=im_ec;
Name='EchoImage';
for k=1:size(im,4)

    % create nifti file
    nifti_file_name= [full_output_folder '/' Name '_' sprintf('%.4d',k) '.nii'];
    sRO = connection.header.encoding.encodedSpace.matrixSize.x*2;
    gadgetron.custom.utils.create_nifti_volume(abs(im(:,:,:,k)),nifti_file_name,reference_header(acquisition),[res(1) res(2) res(3)],[sRO sRO sRO])
    disp("create_and_save_nifti setup...")
end

ListFiles=cellstr(spm_select('FPList',full_output_folder,'^EchoImage'));
Output='Echoes4D'
create4Ddata(ListFiles,Output)

end

%% suppport functions
function reference = reference_header(acquisition)
% We pick the first header from the header arrays - we need it to initialize the image meta data.
reference = structfun(@(arr) arr(:, 1)', acquisition.data.header, 'UniformOutput', false);
end



function create4Ddata(ListFiles,Output)
matlabbatch{1}.spm.util.cat.vols = ListFiles;
matlabbatch{1}.spm.util.cat.name = Output;
matlabbatch{1}.spm.util.cat.dtype = 0;
matlabbatch{1}.spm.util.cat.RT = NaN;

spm_jobman('run', matlabbatch);
clear matlabbatch
end

