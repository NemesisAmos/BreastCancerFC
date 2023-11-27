clear all;clc;
%generate statistical correlation pictures figure
globalsignal='nGSR';
thres_type={'zPos_Bin'};
group_fix={''};
data_root=strcat('H:\all_subjects_res\statistical_res\',globalsignal,'\');
vol_size=[61,73,61];
%hyperparameter
predefine_rect= 1.0e+03 *[0.7 1.3575 0.5680 0.1360];%store xmin ymin width height of color bar location
white_box=255.*ones(predefine_rect(4),predefine_rect(3),3);
img_height=1500;
img_width=2000;
white_img=255.*ones(img_height,img_width,3);
display_px=5;
display_pn=2;
display_nx=-5;
display_nn=-2;
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

for temp_group=1:length(group_fix)
    jagsaw_col_pic_num=1;
    jawsaw_row_pic_num=ceil(length(thres_type)/jagsaw_col_pic_num);
    for temp_row=1:jawsaw_row_pic_num
        eval(['Functional_voxel_corrmap_row',num2str(temp_row),'=[];']);
    end
    Functional_voxel_corrmap=[];
    % Functional_voxel_corrmap=255*ones(img_height*jawsaw_row_pic_num,img_width*jagsaw_col_pic_num,3);
    for temp_thres=1:length(thres_type)
        temp_data=zeros(vol_size);
        temp_pos_path=strcat(data_root,'FCS_',thres_type{temp_thres},'_group_diff',group_fix{temp_group},'\',globalsignal,'_FCS_significant_positive_group_diff_SnPM.nii');
        if exist(temp_pos_path,'file')
            temp_pos_data=rest_to4d(temp_pos_path);
            temp_data=temp_data-temp_pos_data;%fits Breast cancer-NC
             save_path=strcat(temp_pos_path(1:end-4),'_pic.nii');
             rest_Write4DNIfTI(-temp_pos_data,Header,save_path);
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
        temp_neg_path=strcat(data_root,'FCS_',thres_type{temp_thres},'_group_diff',group_fix{temp_group},'\',globalsignal,'_FCS_significant_negative_group_diff_SnPM.nii');
        if exist(temp_neg_path,'file')
            temp_neg_data=rest_to4d(temp_neg_path);
            temp_data=temp_data+temp_neg_data;%fits Breast cancer-NC
            save_path=strcat(temp_neg_path(1:end-4),'_pic.nii');
            rest_Write4DNIfTI(temp_neg_data,Header,save_path);
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
        temp_row=ceil(temp_thres/jagsaw_col_pic_num);
        temp_col=mod(temp_thres,jagsaw_col_pic_num);
        if temp_col==0
            temp_col=jagsaw_col_pic_num;
        end
        if length(unique(temp_data))>1 && length(find(abs(unique(temp_data))>=T_thres))>0
            temp_save_path=strcat(data_root,'FCS_',thres_type{temp_thres},'_group_diff',group_fix{temp_group},'\',globalsignal,'_FCS_significant_group_diff_SnPM.nii');
            Header.fname=temp_save_path;
            rest_Write4DNIfTI(temp_data,Header,temp_save_path);
            if plot_to_surface
                output_txt_path = amos_vol2surf_txt(temp_save_path, ref_txt_path, ref_nii_path);
                savename = strcat(temp_save_path(1:end-4),'_flatsurf.jpg');
                amos_fullview_brainsurface(output_txt_path, savename, color_num, color_map, alpha, show_colorbar, display_range);
            else
                load('Cfg.mat');
                %set display range
                EC.lot.view=3;
                EC.msh.alpha=alpha;
                EC.vol.display=1;
                EC.vol.color_map=100;
%                 EC.vol.CMt=100;
                EC.vol.nn=display_nn;
                EC.vol.nx=display_nx;
                EC.vol.pn=display_pn;
                EC.vol.px=display_px;
                save('Cfg.mat','EC');
                savename = strcat(temp_save_path(1:end-4),'.jpg');
                BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',temp_save_path,'Cfg.mat',savename);
                close(gcf);
            end
            temp_img=imread(savename);
            temp_img(predefine_rect(2)+1:predefine_rect(2)+predefine_rect(4),predefine_rect(1)+1:predefine_rect(1)+predefine_rect(3),:)=white_box;
            eval(['Functional_voxel_corrmap_row',num2str(temp_row),'=[Functional_voxel_corrmap_row',num2str(temp_row),32,'temp_img];']);
        else
            eval(['Functional_voxel_corrmap_row',num2str(temp_row),'=[Functional_voxel_corrmap_row',num2str(temp_row),32,'white_img];']);
        end
        
        if temp_col<jagsaw_col_pic_num
            add_col_num=jagsaw_col_pic_num-temp_col;
            for temp_add_col=1:add_col_num
                eval(['Functional_voxel_corrmap_row',num2str(temp_row),'=[Functional_voxel_corrmap_row',num2str(temp_row),32,'white_img];']);
            end
        end
    end
    
    for temp_row=1:jawsaw_row_pic_num
        eval(['Functional_voxel_corrmap=[Functional_voxel_corrmap; Functional_voxel_corrmap_row',num2str(temp_row),'];']);
    end
    %save result figure
    if plot_to_surface
        Figure_savename=strcat(data_root,'Figure',group_fix{temp_group},'_all_voxel_FCS_diffmap_flatsurf.jpg');
    else
        Figure_savename=strcat(data_root,'Figure',group_fix{temp_group},'_all_voxel_FCS_diffmap.jpg');
    end
    imwrite(Functional_voxel_corrmap,Figure_savename);
end