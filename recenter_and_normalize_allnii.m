%here I combined the re_center_allnii.m and normalize_allnii_recentered.m together
clear
%direc =  'F:\PET_Huashan\allnii';
%direc =  'C:\MMD68\ncnii';
direc =  'C:\Users\admin\Desktop\PET_process\HUASHAN\';


target_files='*.gz';
cd(direc)
name = sprintf('%s', target_files);
allfilenames = dir(name);
n = length(allfilenames);
for k=1:n
k;
v = spm_vol(allfilenames(k).name);
img = spm_read_vols(v);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   add codes to make sure v.mat is a 
    %   diag matrix
    %      or
    %   v.mat(:,4)=[-v.dim'/2; 1].* [sqrt(sum(v.mat(1:3,1:3).^2)) 1]';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v.mat(:,4)=[-v.dim'/2; 1].* diag(v.mat);
    tmp = v.fname;
%     tmp = strrep(tmp,'NII','Recentered'); %destination folder
    v.fname = strcat(tmp(1:end-4),'_recentered.nii'); 
%cd F:\PET_Huashan\allnii_recentered
    spm_write_vol(v,img);
%cd(direc)
    disp(strcat(allfilenames(k).name,' done.'));
   
end


%clear
%direc =  'F:\PET_Huashan\allnii_recentered';
path = 'C:\Program Files\MATLAB\spm12\spm12\tpm\TPM.nii';
voxel_size = [3 3 3];
target_files='*_recentered.nii';
cd(direc)
name = sprintf('%s', target_files);
allfilenames = dir(name);
n = length(allfilenames);
for k=1:n
k;
% v = spm_vol(allfilenames(k).name);
% img = spm_read_vols(v);
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {strcat(allfilenames(k).name,',1');};
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {strcat(allfilenames(k).name,',1');};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {path};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                             78 76 85];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = voxel_size;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'norm';
spm_jobman('run',matlabbatch);
end


