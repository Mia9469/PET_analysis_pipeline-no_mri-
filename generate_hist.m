% load('C:\Users\admin\Desktop\PET_process\HUASHAN\original\mean_maps\fdg_all_HUASHAN.mat')
load('C:\Users\admin\Desktop\PET_process\ogi6f.mat')
data = zeros(size(all_fdg_masked(:,:,:,1)),'single');
for i = 1:k
    data = data + all_fdg_masked(:,:,:,i);
end
all_fdg_masked_mean=data./k;
p=permute(all_fdg_masked_mean,[2,1,3]);
 p_r=p./54400;  %HUASHAN
% p_r=p./28800;   %GERMAN
glcox=p_r.*ogi6f;
% p1=p_r(:);
% p1(p1<=0)=[];
p2=glcox(:);
p2(p2<=0)=[];
% figure;hist(p1(:),0:0.001:0.6)
figure;hist(p2(:),0:0.001:0.6)
% figure;image_series(p_r)
% colorbar
figure;image_series(glcox,0.55)
colorbar