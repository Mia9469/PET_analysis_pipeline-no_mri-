clear;close all
% load('fdg_all_GERMAN.mat')
load('C:\Users\admin\Desktop\PET_process\HUASHAN\original\mean_maps\fdg_all_HUASHAN.mat')
load('C:\Users\admin\Desktop\PET_process\ogi6f.mat')

%% plot all
% y=zeros(271633,1);
% for i = 1:k
%     figure();
%     p=permute(all_fdg_masked(:,:,:,i),[2,1,3]);
%     p_r=p./54400.*ogi6f;
%     p2=p_r(:);
% %     p2(p2<=0)=[];
% %     [nh,xh]=hist(p2(:),0:0.001:0.6);
%     plot(p2);%hold on;
%     ylim([0,1])
%     y=y+p2;
% end
% figure();plot(y./k)
% ylim([0,1])

%% plot hist
data = zeros(size(all_fdg_masked(:,:,:,1)),'single');
for i = [1 3 6:13 15 17:20]%1:k
    data = data + all_fdg_masked(:,:,:,i);
    p=permute(all_fdg_masked(:,:,:,i),[2,1,3]);
    p_r=p./54400;  %HUASHAN
%     p_r=p./28800;   %GERMAN
    glcox=p_r.*ogi6f;
    % p1=p_r(:);
    % p1(p1<=0)=[];
    p2=glcox(:);
    p2(p2<=0)=[];
    % figure;hist(p1(:),0:0.001:0.6)
    figure;hist(p2(:),0:0.001:0.6)
    figure;image_series(glcox,0.55)
    colorbar
%     savefig([num2str(i) '.fig'])
end
all_fdg_masked_mean=data./15;%k;
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
% savefig( 'mean.fig')
%% plot ind br_std/mean
% load('id_82br_stdpetglcx.mat')
% load('id_82br_mpetglcox.mat')
% br_focus=1:82;
% MMD_hemo=MMDhemo_19_std(:,br_focus)./MMD_hemorrhagic_19(:,br_focus);
% MMD_isch=MMDisch_44_std(:,br_focus)./MMD_ischemic_44(:,br_focus);
% CVD=CVD_14_std(:,br_focus)./CVD_14(:,br_focus);
% Hc=Hc_20_std(:,br_focus)./Hc_20(:,br_focus);
% 
% data=MMD_hemo;
% s = size(data);
% for i = 1:s(1)
%     figure();bar(Hc(i,:));
% end