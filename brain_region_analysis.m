load D:\bigbrain\mask\BrainRegionsMask_3mm.mat
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

bregion41=PET_brain_region(41) ;
for i=1:length(bregion41.numb)
             x=bregion41.regions{i};    
            volum41(i)=sum(s1d_3d(x)).*27e-3;%.*27e-3;% 1mm=1e-3 m; 1mm^3=1e-9 m^3; 1m^3=1e3 L=1e-6 mL; so 1mm^3=1e-3 mL 
            mod2f41=s1d_3d(dk_ebrain.rategmmap.*bregion41.regions{i},0); mod2f41=double(mod2f41);
            mod2glcx41=s1d_3d(dk_ebrain.modbg.*bregion41.regions{i},0); mod2glcx41=double(mod2glcx41);
            petglcx41=s1d_3d(dk_ebrain.petgm.*bregion41.regions{i},0); petglcx41=double(petglcx41);
            pet1glcx41=s1d_3d(ebrainf119.eclose_gm.*bregion41.regions{i},0); pet1glcx41=double(pet1glcx41);
         %petglcx41 is same as pet1glcx41.....double conform them because
         %generated from mod1 and mod2 programs, but using same set of PET data
            mod1glcx41=s1d_3d(ebrainf119.glcg3mm.*bregion41.regions{i},0); mod1glcx41=double(mod1glcx41);
            Ns2=s1d_3d(dk_ebra in.gmcellpop.*bregion41.regions{i},0); Ns2=double(Ns2);     Nsc2(i)=sum(Ns2);
            %find those volxels without any neuron inside and remove them to avoide artifact on sparseness
            n0=(find(Ns2==0)); Ns2(n0)=[];mod2glcx41(n0)=[]; mod1glcx41(n0)=[]; petglcx41(n0)=[]; pet1glcx41(n0)=[]; mod2f41(n0)=[];
                nmod1glcx41{i}=double([mod1glcx41,Ns2]); % model1: how many voxels for each given glcx value for a broadmann region
                nmod2f41{i}=double([mod2f41,Ns2]); %  model2: how many voxels for each given activity rate value for a broadmann region
                nmod2glcx41{i}=double([mod2glcx41,Ns2]); % model2: how many voxels for each given glcx value for a broadmann region
                npet1glcx41{i}=double([pet1glcx41,Ns2]); %  in pet data, how many voxels for each given glcx value for a broadmann region
                npetglcx41{i}=double([petglcx41,Ns2]); %  in pet data, how many voxels for each given glcx value for a broadmann region
                    mmod1glcx41(i)=sum(nmod1glcx41{i}(:,1).*nmod1glcx41{i}(:,2))./Nsc2(i); %mean glcx value for each broadmann region
                    stdmod1BAglc41(i)=sqrt(sum(Ns2.*(mod1glcx41-mmod1glcx41(i)).^2)./Nsc2(i));
          
                 [nfs1,xfs1]=hist(mod2f41,0:0.1:25);xhist1=[xfs1',nfs1'];
                      fhist41{i}=xhist1; %voxel-based glcx distribution
                      sx41(i)=sparseness_hist(xhist1);
                  
            %us voxel instead of individual cell
              mmod1glcx41(i)=mean(mod1glcx41);  stdmod1glcx41(i)=std(mod1glcx41); 
              mmod2glcx41(i)=mean(mod2glcx41);  stdmod2glcx41(i)=std(mod2glcx41); 
              mpetglcx41(i)=mean(pet1glcx41);   stdpetglcx41(i)=std(pet1glcx41); 

              mmod2f41(i)=sum(nmod2f41{i}(:,1).*nmod2f41{i}(:,2))./Nsc2(i); %mean glcx value for each broadmann region
                         stdmod2BAf41(i)=sqrt(sum(Ns2.*(mod2f41-mmod2f41(i)).^2)./Nsc2(i));
      
             [ng1,xg1]=hist(glcxsig,0:0.001:0.5);
             ghist41{i}=[xg1',ng1'];  %voxel-based glcx distribution
       
                         %                     mmod2glcx41(i)=sum(nmod2glcx41{i}(:,1).*nmod2glcx41{i}(:,2))./Nsc2(i); %mean glcx value for each broadmann region
%                     stdmod2BAglc41(i)=sqrt(sum(Ns2.*(mod2glcx41-mmod2glcx41(i)).^2)./Nsc2(i));
%                     mpetglcx41(i)=sum(npetglcx41{i}(:,1).*npetglcx41{i}(:,2))./Nsc2(i); %mean glcx value for each broadmann region
%                     stdpetBAglc41(i)=sqrt(sum(Ns2.*(petglcx41-mpetglcx41(i)).^2)./Nsc2(i));
%                     mpet1glcx41(i)=sum(npet1glcx41{i}(:,1).*npet1glcx41{i}(:,2))./Nsc2(i); %mean glcx value for each broadmann region
%                     stdpet1BAglc41(i)=sqrt(sum(Ns2.*(pet1glcx41-mpet1glcx41(i)).^2)./Nsc2(i));
end
glcBA41=[[1:41]',volum41',mpetglcx41',stdpetglcx41',mmod1glcx41',stdmod1glcx41', mmod2glcx41',stdmod2glcx41',mmod2f41',stdmod2BAf41',sx41'];

