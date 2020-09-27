clear;
direc = 'G:\2017DICOM\';
target_files='*PET*';
cd(direc)
name = sprintf('%s', target_files);
allfilenames = dir(name);
n = length(allfilenames);
for i=9:n
    sub=allfilenames(i).name;
    cd([direc sub])
    name1 = sprintf('%s', '*2017*');
    allfilenames1 = dir(name1);
    sub1=allfilenames1(1).name;
    cd([direc sub '\' sub1])
    name2 = sprintf('%s', target_files);
    allfilenames2 = dir(name2);
    sub2=allfilenames2(1).name;
    pathname=[direc sub '\' sub1 '\' sub2];
    
    
    sdir=dir([pathname '\' '*.IMA']);%ѡȡIMA
    imgfile='';
    for i=1:length(sdir)
        imgfile=[imgfile;pathname '\' sdir(i).name];
    end

    read_dicom
end