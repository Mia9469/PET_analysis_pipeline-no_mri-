clear;
close all;
load('fdg_all_HUASHAN.mat')
% load('fdg_all_MMD_ischemic.mat')
load('C:\Users\admin\Desktop\PET_process\ogi6f.mat')
l_r_compare=1;

data = zeros(size(all_fdg_masked(:,:,:,1)),'single');
for n = 1:k
    data=all_fdg_masked(:,:,:,n);
    k1=54400;
    bregion=PET_brain_region(82) ;
close all;
% err=zeros(1,length(bregion.numb));
l_r_wb=zeros(73,61,61);
for i=[1:11 42:52]%1:length(bregion.numb)%[12 17 48 49 54 77]
             x=bregion.regions{i};  
%              volum41(i)=sum(s1d_3d(x)).*27e-3;%.*27e-3;% 1mm=1e-3 m; 1mm^3=1e-9 m^3; 1m^3=1e3 L=1e-6 mL; so 1mm^3=1e-3 mL 
%              petglcx41=s1d_3d(permute(data,[2,1,3]).*ogi6f./k1.*bregion.regions{i},0); petglcx41=double(petglcx41);
             ver{i}=double(permute(data,[2,1,3]).*ogi6f./k1.*bregion.regions{i});
            %find those volxels without any neuron inside and remove them to avoide artifact on sparseness
%              petglcx41(find(petglcx41==0))=[]; 
             if l_r_compare 
                 if i>41
%                      l_r_ver=l_r_ver+ver{i}-ver{i-41};
                     l_r_ver1=(ver{i}-ver{i-41});
                     l_r_ver2=(ver{i}+ver{i-41});
                     l_r_ver2(isnan( l_r_ver2)==1)=0;
                     l_r_ver1(isinf( l_r_ver1)==1)=0;
                     maxdiff(i)=max(max(max(l_r_ver1./l_r_ver2)));
                     l_r_wb=l_r_wb+l_r_ver1./l_r_ver2;
                 end
             end             

%                 figure();
%                 hist( petglcx41(:),0:0.001:0.6);
%                 hold on;
%                 [nh,xh]=hist( petglcx41c(:),0:0.001:0.6);
%                 plot(xh,nh,'g') 
%                 title(i)
                %savefig([num2str(i) '.fig'])
            
            %us voxel instead of individual cell
%               mpetglcx41(i)=mean(petglcx41);   stdpetglcx41(i)=std(petglcx41);
            
end

Hc{1,n}=l_r_wb;
Hc{2,n}=maxdiff;
% figure;image_series(CVD{1,1},0.5);colorbar
end