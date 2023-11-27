clear all;clc;
%generate statistical correlation pictures figure
globalsignal='nGSR';
thres_type='zPos_Bin_';
Distance_Range='000-Inf';
group_fix={''};
LoadRoot=strcat('H:\all_subjects_res\Results\');
SaveRoot=strcat('H:\all_subjects_res\statistical_res\',globalsignal,'\');
excelroot='H:\all_subjects_res\all_subj_info.xlsx';
[behavioraldata,subjlist]=xlsread(excelroot);
age=behavioraldata(:,1);
group=behavioraldata(:,2);%do not consider firstep/medicate effect
FCRI=behavioraldata(:,3);
medicated=behavioraldata(:,5);
subjlist=subjlist(2:length(behavioraldata)+1);%remove variable name
maskroot=strcat('E:\Matlab toolboxes\SeeCAT\templates\GMwithoutCER_61x73x61.nii');%only calculate FCS on GM mask
[maskdata,~,~,Header]=rest_to4d(maskroot);
Header.dt=[16,0];
[X,Y,Z]=size(maskdata);
maskdata=reshape(maskdata,1,[]);
maskdata=logical(maskdata);

%hyperparameter
predefine_rect= 1.0e+03 *[0.7 1.3575 0.5680 0.1360];%store xmin ymin width height of color bar location
white_box=255.*ones(predefine_rect(4),predefine_rect(3),3);
img_height=1500;
img_width=2000;
white_img=255.*ones(img_height,img_width,3);
display_px=2;
display_pn=0;
display_nx=-2;
display_nn=0;
plot_to_surface=0;

N = 1000;

color_num = 64;
alpha = 1;
show_colorbar = true;
ref_txt_path = 'E:\Matlab toolboxes\plotBrainSurfaceRes\GMwithoutCER_61x73x61_UniqueVoxelValue.txt';
ref_nii_path = 'E:\Matlab toolboxes\plotBrainSurfaceRes\GMwithoutCER_61x73x61_UniqueVoxelValue.nii';
display_range = [display_nx,display_nn,display_pn,display_px];

[~,~,~,Header]=rest_to4d('E:\MATLAB toolboxes\SeeCAT\templates\GMwithoutCER_61x73x61.nii');
Header.dt=[16,0];
T_thres=0.001;%no need to set T treshold allready thresholded in generation


%load data and provide center level statistical result
hc_subjlist=subjlist(find(group==1));
bc_subjlist=subjlist(find(group==2));
hc_all_voxelDeg=[];
bc_all_voxelDeg=[];
for temp_subj=1:length(subjlist)
    %disp(['Current processing subject: ',subjlist{temp_subj} ' Range: ' Distance_Range{tempRange} ' Frequency: ' all_freq{temp_freq}])
    if strcmp(globalsignal,'GSR')
        smoothfix='global';
    else
        smoothfix='';
    end
    if group(temp_subj)==1
        subjPath=strcat(LoadRoot,'VoxelDeg_FunImgARWSD',smoothfix,'CFB_SnPM/',thres_type,Distance_Range,'/Group1/',thres_type,Distance_Range,'_',subjlist{temp_subj},'.nii');
        [subjAllVolume,~,~,Header]=rest_to4d(subjPath);
        subjAllVolume=reshape(subjAllVolume,1,[]);
        subjAllVolume(find(isnan(subjAllVolume)))=0;%correct NaN values after seeCAT
        hc_all_voxelDeg=[hc_all_voxelDeg;subjAllVolume(find(maskdata==1))];
    elseif group(temp_subj)==2
        subjPath=strcat(LoadRoot,'VoxelDeg_FunImgARWSD',smoothfix,'CFB_SnPM/',thres_type,Distance_Range,'/Group2/',thres_type,Distance_Range,'_',subjlist{temp_subj},'.nii');
        [subjAllVolume,~,~,Header]=rest_to4d(subjPath);
        subjAllVolume=reshape(subjAllVolume,1,[]);
        subjAllVolume(find(isnan(subjAllVolume)))=0;%correct NaN values after seeCAT
        bc_all_voxelDeg=[bc_all_voxelDeg;subjAllVolume(find(maskdata==1))];  
    else           
    end
end

%% save MDD and NC group level results
%% mean map
BC_meanmap=zeros(1,length(maskdata));
BC_meanmap(find(maskdata==1))=mean(hc_all_voxelDeg);
BC_meanmap(find(isnan(BC_meanmap)))=0;
BC_map_scale=BC_meanmap;%use original mean of all subject zscore value
BC_map_scale(find(isnan(BC_map_scale)))=0;
BC_map_scale=reshape(BC_map_scale,[X,Y,Z]);
tempsavename=strcat(SaveRoot,globalsignal,'_',thres_type,'VoxelDeg_NC_meanmap.nii');
Header.fname=tempsavename;%keep Header info same
rest_Write4DNIfTI(BC_map_scale,Header,tempsavename);
if exist(tempsavename,'file')
    temp_data=rest_to4d(tempsavename);
    save_path=strcat(tempsavename(1:end-4),'_pic.nii');
    rest_Write4DNIfTI(temp_data,Header,save_path);
    load('Cfg.mat');
    %set display range
    EC.lot.view=3;
    EC.msh.alpha=alpha;
    EC.vol.display=1;
    EC.vol.color_map=100;
    EC.vol.nn=display_nn;
    EC.vol.nx=display_nx;
    EC.vol.pn=display_pn;
    EC.vol.px=display_px;
    save('Cfg.mat','EC');
    savename = strcat(save_path(1:end-4),'.jpg');
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',save_path,'Cfg.mat',savename);
    close(gcf);
end
        
NC_meanmap=zeros(1,length(maskdata));
NC_meanmap(find(maskdata==1))=mean(bc_all_voxelDeg);
NC_meanmap(find(isnan(NC_meanmap)))=0;
NC_map_scale=NC_meanmap;%use original mean of all subject zscore value
NC_map_scale(find(isnan(NC_map_scale)))=0;
NC_map_scale=reshape(NC_map_scale,[X,Y,Z]);
tempsavename=strcat(SaveRoot,globalsignal,'_',thres_type,'VoxelDeg_BC_meanmap.nii');
Header.fname=tempsavename;%keep Header info same
rest_Write4DNIfTI(NC_map_scale,Header,tempsavename);
if exist(tempsavename,'file')
    temp_data=rest_to4d(tempsavename);
    save_path=strcat(tempsavename(1:end-4),'_pic.nii');
    rest_Write4DNIfTI(temp_data,Header,save_path);
    load('Cfg.mat');
    %set display range
    EC.lot.view=3;
    EC.msh.alpha=alpha;
    EC.vol.display=1;
    EC.vol.color_map=100;
    EC.vol.nn=display_nn;
    EC.vol.nx=display_nx;
    EC.vol.pn=display_pn;
    EC.vol.px=display_px;
    save('Cfg.mat','EC');
    savename = strcat(save_path(1:end-4),'.jpg');
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',save_path,'Cfg.mat',savename);
    close(gcf);
end


