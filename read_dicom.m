% List of open inputs
% DICOM Import: DICOM files - cfg_files
% DICOM Import: Output directory - cfg_files
nrun = 1; % enter the number of runs here
jobfile = {'C:\Users\admin\Desktop\PET_process\code\read_dicom_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    inputs{1, crun} = cellstr(imgfile); % DICOM Import: DICOM files - cfg_files
    inputs{2, crun} = cellstr('C:\Users\admin\Desktop\PET_process\epil'); % DICOM Import: Output directory - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});
