%% COrrect kspace and trajectory 

function [kdata_new traj_bart_ec_new]=correct_traj_kdata_clean(kdata, traj_bart_ec, rp_file, nRad, nTR, nContrast, nRO,ratio) 
%% load motion parameters estimated from Image  TR 
if (nargin < 8)
    ratio=2;
end
M=load(rp_file);


%% multiply by 2 the translations because motion estimated on TR images for which the resolution is twice smaller
M(:,1:3)=ratio*M(:,1:3);

for k=1:nRad/nTR-length(M)
    M=[M;M(end,:)];
end
%%% interpolation for each projection because we have an information only
%%% every one TR
M=interp1(1:nTR:nRad,M,1:nRad); 
[ind]=find(isnan(M),1);

for k=0:size(M,1)-ind
    M(ind+k,:)=M(ind-1,:);
end

%% create transformation matrix 
for k=1:size(M,1)
    Mmat(:,:,k)=spm_matrix(M(k,:));
end

%% transform into the reference of the k-space 
for k=1:size(Mmat,3)
    Pass=[-1 0 0 0 ; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    alignMats(:,:,k)=Pass*Mmat(:,:,k)*inv(Pass);
end


%% kspace trajectory correction
for ec=1:nContrast
    
    for r=1:nRad
        
        %%% rotation
        traj_bart_ec_new(:,:,r,1,1,ec)=(alignMats(1:3,1:3,r))*traj_bart_ec(:,:,r,1,1,ec);
        
        %%% translation
        kdata_new(1,:,r,:,1,ec)=kdata(1,:,r,:,1,ec).*exp((-2*1i*pi/nRO)*(alignMats(1,4,r)*traj_bart_ec_new(1,:,r,1,1,ec) + (alignMats(2,4,r))*traj_bart_ec_new(2,:,r,1,1,ec) + (alignMats(3,4,r))*traj_bart_ec_new(3,:,r,1,1,ec)));
        
    end
end

disp('modif done')

end