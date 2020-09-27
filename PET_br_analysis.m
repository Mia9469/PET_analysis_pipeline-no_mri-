clear;close all
%% load files
load C:\Users\admin\Desktop\PET_process\mask\BrainRegionsMask_3mm.mat
load('C:\Users\admin\Desktop\PET_process\ogi6f.mat')
load('C:\Users\admin\Desktop\PET_process\mask\cellpopulation1mm3mm.mat', 'ebrain3mm')
%% set test parameters
compare_hist=1;
plot_comparison=0;
corr_analysis=1;
l_r_compare=0;
 compute_sparseness=0;
% n=6;

%% set masks&regions
     gm=graymask;
     wm=whitemask;
     mask_regions=mymask;
    mask_regions=permute(mask_regions,[2,1,3]);
    gm=permute(gm,[2,1,3]);
    wm=permute(wm,[2,1,3]);
BroadmannL=[1;4;5;6;7;8;9;10;11;13;17;18;19;20;21;22;23;24;25;30;31
32;34;36;37;38;39;40;41;44;45;46;47;48;49;50;51;52;53;54;55];
Li=BroadmannL;
BroadmannR=BroadmannL+100;
OtherRegions=[96; 97; 98]; 

%% load & adjust pet mean data
load('all_fdg_gwm.mat')
load('mfdg_epilepsia.mat')
% load('mfdg_epilepsia')
% load('C:\Users\admin\Desktop\PET_process\Denmark_PET_data\denmark.mat')
% load('HUASHAN_MEAN_15.mat')
% data{1}=all_fdg_masked_mean;
% load('GERMANY_MEAN_21.mat')
% data{2}=all_fdg_masked_mean;

% data{1}= all_fdg_gwm.mHcontrol;
% data{1}=glcx;
% data{1}=all_fdg_gwm.mMMD_hemorrhagic;
data{1}=(all_fdg_gwm.mMMD_ischemic+all_fdg_gwm.mMMD_hemorrhagic)/2;
% data{1}=mean_Fdg.epilepsia;
% data{1}=all_fdg_gwm.G22;
data{2}=all_fdg_gwm.mCVD;

% k=54400;
k=ones(1,2);
for i=1:2
%     if data{i}==all_fdg_gwm.G22%Geo
%       k(i)=28800;
    if data{i}==data{2}
        k(i)=28800;
    else
      k(i)=54400;
    end
    data{i}=permute(data{i},[2,1,3]).*ogi6f;
%     data{i}=permute(data{i},[2,1,3]);
end
% data{1}=permute(data{1},[2,1,3]).*ogi6f;
%% calcultation in brain regions
% epil_diff=[5     8     9    46    49];

bregion=PET_brain_region(82) ;
close all;
err=zeros(1,length(bregion.numb));
for i=1:length(bregion.numb)
             x=bregion.regions{i};  
             volum41(i)=sum(s1d_3d(x)).*27e-3;%.*27e-3;% 1mm=1e-3 m; 1mm^3=1e-9 m^3; 1m^3=1e3 L=1e-6 mL; so 1mm^3=1e-3 mL 
%              ver{i}=double(data{1}./k(1).*bregion.regions{i});
             petglcx41=s1d_3d(data{1}./k(1).*bregion.regions{i},0); petglcx41=double(petglcx41);
%              Ns2=s1d_3d(ebrain3mm.gmcellpop.*bregion.regions{i},0); Ns2=double(Ns2);     Nsc2(i)=sum(Ns2);
%              glcx_cell=s1d_3d(data{1}./k(1)./ebrain3mm.gmcellpop.*bregion.regions{i },0);
%             %find those volxels without any neuron inside and remove them to avoide artifact on sparseness
%              n0=(find(Ns2==0)); Ns2(n0)=[]; glcx_cell(n0)=[]; petglcx41(n0)=[]; 
             petglcx41(find(petglcx41==0))=[];  sumpetglcx41(i)=sum(petglcx41);
             mpetglcx41(i)=mean(petglcx41); stdpetglcx41(i)=std(petglcx41); 
             
            if compare_hist
             [nhh,xh]=hist( petglcx41(:),0:0.001:0.6);
             nhh=petglcx41(:); 
%              epil(:,i)=nhh';
                %load control data
                petglcx41c=s1d_3d(data{2}./k(2).*bregion.regions{i},0); petglcx41c=double(petglcx41c);
                verc{i}=double(data{2}./k(2).*bregion.regions{i});
                petglcx41c(find(petglcx41c==0))=[]; 
                                
                [nh,xh]=hist( petglcx41c(:),0:0.001:0.6);
                nh=petglcx41c(:);
%                 Hc(:,i)=nh';
                if  plot_comparison == 1  
                figure();
                plot(nh,'k')
                hold on;
                plot(nhh,'g')
                title(i)
%                 savefig([num2str(i) '.fig'])
               end
                mpetglcx41c(i)=mean(petglcx41c);   stdpetglcx41c(i)=std(petglcx41c);
         
%                 %compute significant difference
%                [err(i),p(i)]=ttest(petglcx41(:),petglcx41c(:));
%                 if find_diff
%                     [nh,xh]=hist( abs(petglcx41(:)-petglcx41c(:)),0:0.001:0.1);
%                     if max(nh>=n)==1
%                        err(i)=max(xh(nh>=n));
%                     end
%                 end
            end  
            
%% sparseness
            if compute_sparseness
                cell=0;
                if cell==1
                mpetglcx41_cell(i)=sumpetglcx41(i)/Nsc2(i);
                sparseness(i)= (1-(sum(petglcx41)/Nsc2(i))^2/sum(glcx_cell.^2.*Ns2/Nsc2(i)))/(1-1/Nsc2(i));
                else
            %us voxel instead of individual cell
                  stdpetglcx41(i)=std(petglcx41); n=size(petglcx41);
                sp(i)= (1-sum(petglcx41)^2/n(1)/sum(petglcx41.^2))/(1-1/n(1));
                end
            end
%               [rho(i),pval(i)]=corr(nhh',nh');

%% left_right half compare
             if l_r_compare 
                 if i>41
                     l_r_ver=l_r_ver+(ver{i}-ver{i-41})./(ver{i}+ver{i-41});
%                      l_r_verc{i-41}=(verc{i}-verc{i-41})./(verc{i}+verc{i-41});
                 end
             end
%% correlation analysis             
             if corr_analysis
                 s1=size(nh);s2=size(nhh);
                 xcf1(i) = max(abs(crosscorr(nh',nhh',min(s1(1),s2(1))-1))); 
%                  n=0;
%                  xcfx=[];
%                  xcfy=[];
%                  xcfz=[];
%                  for j=1:61
%                      if sum(sum(ver(:,:,j)))~=0
%                          n = n+1;
%                          xcfx(n) = corr2(ver(:,:,j),verc(:,:,j));           
%                      end
%                  end 
%                  ver1=permute(ver,[1,3,2]);
%                  ver1c=permute(verc,[1,3,2]);
%                  for j=1:61
%                      if sum(sum(ver1(:,:,j)))~=0
%                          n = n+1;
%                          xcfy(n) = corr2(ver1(:,:,j),ver1c(:,:,j));           
%                      end
%                  end 
%                   ver2=permute(ver,[3,2,1]);
%                   ver2c=permute(verc,[3,2,1]);
%                  for j=1:73
%                      if sum(sum(ver2(:,:,j)))~=0
%                          n = n+1;
%                          xcfz(n) = corr2(ver2(:,:,j),ver2c(:,:,j));           
%                      end
%                  end 
%                  xcf2x(i) = mean(xcfx);
%                  xcf2y(i) = mean(xcfy);
%                  xcf2z(i) = mean(xcfz);
%                  
%                  xcf2(i) = mean([xcf2x(i),xcf2y(i),xcf2z(i)]);
             end
             
%                     mpetglcx41(i)=sum(npetglcx41{i}(:,1).*npetglcx41{i}(:,2))./Nsc2(i); %mean glcx value for each broadmann region
%                     stdpetBAglc41(i)=sqrt(sum(Ns2.*(petglcx41-mpetglcx41(i)).^2)./Nsc2(i));

end
% % glcBA41=[[1:41]',volum41',mpetglcx41',stdpetglcx41']; 
%glcBA82=[[1:82]',volum41',mpetglcx41',stdpetglcx41']; 
% figure();
% bar(mpetglcx41,'c');  
% hold on;  
% errorbar(mpetglcx41,stdpetglcx41,'k','LineStyle','none'); 
% hist_corr={};
% hist_corr.CVD_IMMD=xcf1;
% 
figure();plot(xcf1,'ko');hold on;
plot([0,90],[0.8,0.8],'r-','LineWidth',2);
ylim([0,1])
xlabel('brain region')
ylabel('corr hist')
title('CVD&HMMD')

% 
% figure();plot(xcf2,'ko');
% ylim([0,1])
% xlabel('brain region')
% ylabel('corr image')
% title('Epilepsy&Hc')


%% find brain regions with significant difference
% br=1:82;
% br(err==1)

%% save the values in a cell 

% l_r_verx.CVD=l_r_ver;
% l_r_verx.Hc=l_r_verc;
% mean_std_ox82(1).G22=mpetglcx41c;
% mean_std_ox82(2).G22=stdpetglcx41c;
% mean_std_ox(1).CVD=mpetglcx41;
% mean_std_ox(2).CVD=stdpetglcx41;
% mean_std_ox(1).G22=mpetglcx41c;
% mean_std_ox(2).G22=stdpetglcx41c;
% mean_std_ox(1).Hcontrol=mpetglcx41;
% mean_std_ox(2).Hcontrol=stdpetglcx41;
% mean_std_ox(1).epilesia=mpetglcx41c;
% mean_std_ox(2).epilesia=stdpetglcx41c;
% epil_hcontrol_ox41(1).Hcontrol=mpetglcx41;
% epil_hcontrol_ox41(2).Hcontrol=stdpetglcx41;
% epil_hcontrol_ox41(1).Epil=mpetglcx41c;
% epil_hcontrol_ox41(2).Epil=stdpetglcx41c;

%% plot glucose/volume
% figure;
% w_volum=sum(volum41);
% w_petglcx41=sum(sumpetglcx41);
%  P=polyfit(volum41/w_volum,sumpetglcx41/w_petglcx41,1);
%  x=linspace(0,0.1);
%  plot(x,P(1)*x+P(2),'k-');hold on;
%  plot(volum41/w_volum,sumpetglcx41/w_petglcx41,'ko');
%  title(['Epilepsy ksum=' num2str(P(1))] )
%  xlabel('brain region volume(proportion)')
%  ylabel('metabolism rate of whole region(proportion)')
%  
%  figure;
%  K=polyfit(volum41/w_volum,mpetglcx41,1);
%  x=linspace(0,0.1);
%  plot(x,K(1)*x+K(2),'k-');hold on;
%  plot(volum41/w_volum,mpetglcx41,'ko');
%  ylim([0 0.5])
%  title(['Epilepsy kmean=' num2str(K(1))])
%  xlabel('brain region volume(proportion)')
%  ylabel('mean-value of metabolism rate/pixel')


%% sparseness-cv
% cv=stdpetglcx41./mpetglcx41;%cv_cell=stdpetglcx41_cell./mpetglcx41_cell;
% 
% % for i =1:82;
% %     plot(cv(i),sp(i),'o');hold on;pause(0.1);
% % end
% 
% 
% x=linspace(0.01,1,100);
% 
% figure
% plot(x,x.^2./(1+x.^2),'k-');hold on;
% plot(cv,sp,'o');
% xlim([0 1])
% ylim([0 1])
% 
% xlabel('c.v.')
% ylabel('sparseness')
% title('sparseness-c.v.')
% 
% 
% 
% 
