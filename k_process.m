 clear;
spm_jobman('initcfg')
%#ok<*NBRAK>

%%

% Process pipeline:
% script_process_data！！>PET_non_dynamic_analysis！！>generate_hist


BIS_path = 'C:\yale\bioimagesuite30\';
BET_path = 'C:\yale\bioimagesuite_extra\bin\';
BET_params = [0.5 0.0;0.3 0.2]; % [-f functional first pass, -g functional first pass;-f functional second pass, -g functional second pass];

MNI_file = 'C:\Users\admin\Documents\MATLAB\GARTH\Gemany\Shared_Templates\MNI_3mm_origtemplate\ROI\MNI_T1_3mm_graywhite.nii';
ROI_directory = 'C:\Users\admin\Documents\MATLAB\GARTH\Gemany\Shared_Templates\MNI_3mm_origtemplate\ROI\';
brain_regions_file = ['C:\Users\admin\Documents\MATLAB\GARTH\Gemany\Shared_Templates\brain_regions\BrainRegionsMask_3mm.mat'];
% experiment_folder = 'C:\Users\admin\Desktop\PET_process\adjust_code_Garth\sub001\';
% experiment_folder = 'C:\Users\admin\Desktop\PET_process\HUASHAN\original\';
experiment_folder = 'C:\Users\admin\Desktop\template\';

experiment_folder_fs = flip_slashes(experiment_folder);
MNI_file_fs = flip_slashes(MNI_file);
% Some subjects have upside down functional data, some subjects have errors
% when their registration is done.
% subjects_upside_down = [21 22]; % Removed as has been already fixed, doesn't need to be done again
subjects_upside_down = [];
% subjects_special_transforms = [7 20];
subjects_special_transforms = [];
% special_transforms_folder = 'C:\Gemany\SpecialTransformations\';

% Whether the MNI brain is a GM-WM map (true) or grayscale shades (false)
mni_is_gm_wm = true;

nonlinear_fdg_mni_reg = false;

%% All subs
direc = experiment_folder;
target_files='*.nii';
cd(direc)
name = sprintf('%s', target_files);
allfilenames = dir(name);
n = length(allfilenames);
 for k=1:n

    v = spm_vol(allfilenames(k).name);
%     img = spm_read_vols(v);
    fdg_file=v.fname;
    
    % Set command line parameters
% BET
cmd_bet_set_path = ['path=' BET_path];
cmd_bet_set_output = 'SET FSLOUTPUTTYPE=NIFTI';
% BIS
cmd_bis_set_path = [BIS_path 'setpaths.bat'];

experiment_folder_fs = flip_slashes(experiment_folder);
MNI_file_fs = flip_slashes(MNI_file);
% % nii_folder_fs = flip_slashes(nii_folder);
bis_standard_header = [cmd_bis_set_path '&' BIS_path(1:2) '&cd ' BIS_path 'bis_algorithm&'];
ROI_directory_fs = flip_slashes(ROI_directory);
if exist('fdg_template','var')
    fdg_template_fs = flip_slashes(fdg_template);
    disp('Using custom template fot FDG...');
else
    fdg_template = MNI_file;
    fdg_template_fs = MNI_file_fs;
end

%     fdg_new_size = [3 3 3];
%     cmd_bis_fdg_resample = ['bisexec bis_resampleimage.tcl --blurmode 0 --vox1 ' num2str(fdg_new_size(1)) ' --vox2 ' num2str(fdg_new_size(2)) ' --vox3 ' num2str(fdg_new_size(3)) ' --out ' experiment_folder_fs 'rs_' fdg_file ' --inp ' experiment_folder_fs fdg_file];
%     system([bis_standard_header cmd_bis_fdg_resample]);
%     fdg_file=['rs_' fdg_file];
    
     if nonlinear_fdg_mni_reg
         cmd_bis_compute_pet_mni = ['bisexec bis_nonlinearintensityregister.tcl --metric CC --iterations 50 --out ' experiment_folder_fs num2str(k) 'nonlinear_pet_to_mni.grd --inp ' fdg_template_fs ' --inp2 ' experiment_folder_fs fdg_file];
         cmd_bis_normalize_pet = ['bisexec bis_resliceimage.tcl --out ' experiment_folder_fs 'nn_ar_' fdg_file ' --inp ' MNI_file_fs ' --inp2 ' experiment_folder_fs fdg_file ' --inp3 ' experiment_folder_fs num2str(k) 'nonlinear_pet_to_mni.grd'];
     else
         cmd_bis_compute_pet_mni = ['bisexec bis_linearintensityregister.tcl --mode affine --metric CC --iterations 50 --out ' experiment_folder_fs num2str(k) 'linear_pet_to_mni.matr --inp ' fdg_template_fs ' --inp2 ' experiment_folder_fs fdg_file];
        cmd_bis_normalize_pet = ['bisexec bis_resliceimage.tcl --out ' experiment_folder_fs 'n_ar_' fdg_file ' --inp ' MNI_file_fs ' --inp2 ' experiment_folder_fs fdg_file ' --inp3 ' experiment_folder_fs [num2str(k) 'linear_pet_to_mni.matr']];
     end
     
     system([bis_standard_header cmd_bis_compute_pet_mni]);
     system([bis_standard_header cmd_bis_normalize_pet]);
    
    % Move the SPM transforms
    star_dot_mat_dir = dir([experiment_folder '*.mat']);
    if ~isempty(star_dot_mat_dir)
       movefile([experiment_folder '*.mat'],[experiment_folder 'spm_transforms'],'f');
    end


    all_varnames = {'fdg'};
    native_filenames = {[experiment_folder fdg_file]};


% Save & move native segmented data
 save_str = ['save([experiment_folder ''' fdg_file '_native_data.mat'']'];
for file_saving_index = 1:(length(all_varnames)+1)
    % load variable
    eval([all_varnames{file_saving_index} ' = load_untouch_nii(native_filenames{file_saving_index});']);
    % Move file
    copyfile(native_filenames{file_saving_index},[experiment_folder 'native'],'f');
    % Record to save variable later
    save_str = [save_str ',''' all_varnames{file_saving_index} '''']; %#ok<AGROW>
    % Display
    disp(['Loaded ' all_varnames{file_saving_index} ' in native space']);
end
save_str = [save_str ',''-v7.3'');'];
eval(save_str);
% Clear
for file_clearing_index = 1:length(all_varnames)
    eval(['clear ' all_varnames{file_clearing_index}]);
end


    % ********* Data in MNI space
    mni_filenames = {[experiment_folder 'n_ar_' fdg_file]};


% Save & move data in MNI space
 save_str = ['save([experiment_folder ''' fdg_file '_mni_data.mat'']'];
for file_saving_index = 1:(length(all_varnames)+1)
    % load variable
    eval([all_varnames{file_saving_index} ' = load_untouch_nii(mni_filenames{file_saving_index});']);
    % Move file
    movefile(mni_filenames{file_saving_index},[experiment_folder 'mni'],'f');
    % Record to save variable later
    save_str = [save_str ',''' all_varnames{file_saving_index} '''']; %#ok<AGROW>
    % Display
    disp(['Loaded ' all_varnames{file_saving_index} ' in mni space']);
end
save_str = [save_str ',''-v7.3'');'];
eval(save_str);
% Clear
for file_clearing_index = 1:length(all_varnames)
    eval(['clear ' all_varnames{file_clearing_index}]);
end
 end
