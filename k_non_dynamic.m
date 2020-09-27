% spaces_to_run = {'individual' 'mni'};
clear
spaces_to_run = {'mni'};

%#ok<*NBRAK>

for all_spaces_used = spaces_to_run
    space_used = all_spaces_used{:};
    all_class=['CVD','MMD_hemorrhagic','MMD_ischemic','MMD_pedisch','HUASHAN\original','GERMANY\original','MMD_pedhem'];
   % subject_class ='CVD';
    % subject_class ='MMD_ischemic';
    %subject_class ='MMD_pedisch';
    subject_class ='HUASHAN\original';
    list=load('C:\Users\admin\Desktop\PET_process\HUASHAN\original\mean_maps\allfilenames.mat');
    names=list.allfilenames;
    %% Setup parameters
    
    % space_used = 'individual';
    % space_used = 'mni';
    

%     ec_subjects = [1:4];
%     eo_subjects = [12:13];
    
    TR = 2;
    
%     experiment_parent_folder = '/Volumes/MIA/';
%     output_folder = '/Volumes/MIA/normal_healthyadult_30min/mean_maps/';
%     brain_regions_file = '/Volumes/MIA/GERMANY/Shared_Templates/brain_regions/BrainRegionsMask_2mm.mat';

%     experiment_folder = 'C:\Users\admin\Desktop\PET_process\adjust_code_Garth\sub001\';

experiment_folder = ['C:\Users\admin\Desktop\PET_process\' subject_class '\spm_transforms\'];
output_folder = ['C:\Users\admin\Desktop\PET_process\' subject_class '\mean_maps\'];
    brain_regions_file = ['C:\Users\admin\Desktop\PET_process\mask\BrainRegionsMask_3mm.mat'];

    
    
    if ~exist(output_folder,'dir')
        mkdir(output_folder);
    end
    
    size_mni = [61, 73, 61];
%     size_bold = [64 64 35];

    
    % Note: this study only used isotropic voxels
    switch space_used
        case 'mni'
            voxel_size = 3;
        case 'individual'
            voxel_size = 3;
    end
    

    % Spatial blurring, [fdg;bold] amounts
    fwhm = [12;8];
    blur_sigma = round(fwhm/(2*sqrt(2*log(2))));
    % Spatial blurring / [size sigma]
    spatial_blurring_mm = [blur_sigma*2 blur_sigma];
    
    switch space_used
        case 'mni'
            spatial_blurring = spatial_blurring_mm ./ voxel_size;
        case 'individual'
            spatial_blurring = [];
    end
    
    blur_fdg = false;
    
    % 6 motion parameters + WB, WM and CSF
    num_nuisance_signals = 9;
    
    % Regressions to run.  First one will always be no regression, second
    % one will always be all regressors used.  Further ones can be defined
    % by the variable below if it is nonempty.
    % Order of regressors: global_sig,wm_sig,csf_sig,6xmotion_data
    regressions_to_try = [...
        0 1 1 1 1 1 1 1 1;...
        0 1 0 1 1 1 1 1 1;...
        0 0 1 1 1 1 1 1 1;...
        0 0 0 1 1 1 1 1 1;...
        ];
    num_extra_regression_tests = size(regressions_to_try,1);


 %% Run data
        
    load_preexisting = false;
    
    switch space_used
        case 'individual'
            size_current = size_bold;
        case 'mni'
            size_current = size_mni;
    end
    
    

    
    if strcmpi(space_used,'mni')
        % Brain regions
        load(brain_regions_file);
        gm_regions_map_nocerebellum = ~ismember(mymask,[0,96,97,98]);
    end
    
    % Load in all datas
%     direc = experiment_folder;
% target_files='*gz_mni_data.mat';
% cd(direc)
% name = sprintf('%s', target_files);
% allfilenames = dir(name);
n = length(names);

    if load_preexisting
        load([output_folder 'temp_working_data.mat']); %#ok<UNRCH>
        start_index = subject_index;
        % start_index = 3;
        clear('all_fcd_masked');
    else
        % Allocate
        all_fdg_masked = zeros([size_current,n],'single');
        start_index = 1;
    end
    
    
    for k=1:n
        work_folder = [experiment_folder];
%         target_file=[names(k).name '_mni_data.mat'];
        target_file=['PALS_B12_Brodmann.nii_mni_data.mat'];
        load([work_folder target_file]);
%         load([experiment_folder 'motion_data.mat']);
        
        if strcmpi(space_used,'individual')
            load([experiment_folder 'individual_roi.mat']);
            mymask = double(flip(mymask_mm.img,2)) ./ 211;
            % Keep only regions with full confidence, not border regions
            mymask = (mymask .* (mymask == floor(mymask)));
            gm_regions_map_nocerebellum = ~ismember(mymask,[0,96,97,98]);
        end
        
%         % WB Mask
%         % Mask and save
%         WB_mask = double(bold_wbmask.img);
%         WB_mask(WB_mask ~= 0) = 1;
%         WB_mask(WB_mask == 0) = NaN;
%         all_wb_masks(:,:,:,subject_index) = single(WB_mask);
%         WB_mask_nonans = WB_mask;
%         WB_mask_nonans(isnan(WB_mask_nonans)) = 0;
        
        disp('Getting mean, fdg');
        
        % FDG
        switch space_used
            case 'individual'
                current_fdg = double(flip(fdg.img,2));
            case 'mni'
                current_fdg = double(fdg.img);
        end
        
        % Blur FDG if specified
        if blur_fdg && ~isempty(spatial_blurring)
            current_fdg = smoothn(current_fdg,'gaussian',spatial_blurring(1,1),spatial_blurring(1,2));
        end
        
        xmask=graymask+whitemask;
        newmask=permute(xmask,[2,1,3]);

        current_fdg = current_fdg .* newmask;
        all_fdg_masked(:,:,:,k) = single(current_fdg);

      
%         % Mean
%         all_mean_bold_masked(:,:,:,subject_index) = single(mean(bold_dat,4));
       
    end
        disp(['Finished subject ' num2str(k)]);
        % Save at each loop in case it crashes before finished
        save([output_folder 'fdg_all_' 'HUASHAN'  '.mat'], ...
            'k','all_fdg_masked',...
            '-v7.3')


end
load('C:\Users\admin\Desktop\PET_process\ogi6f.mat')
data = zeros(size(all_fdg_masked(:,:,:,1)),'single');

for i = [1 3:k]
    data = data + all_fdg_masked(:,:,:,i);
end
all_fdg_masked_mean=data./19;
% all_fdg_gwm.Geo=all_fdg_masked_mean;
% p=permute(all_fdg_masked_mean,[2,1,3]);
p_r=p./54400;
% % p_r=p./28800;
glcox=p_r.*ogi6f;
% p1=p_r(:);
% p1(p1<=0)=[];
p2=glcox(:);
p2(p2<=0)=[];
% % figure;hist(p1(:),0:0.001:0.6)
figure;hist(p2(:),0:0.001:0.6)
% % figure;image_series(p_r)
% % colorbar
% % figure;image_series(glcox)
% % colorbar