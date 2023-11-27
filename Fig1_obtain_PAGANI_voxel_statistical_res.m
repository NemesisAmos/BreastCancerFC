clear all;clc;
%% hyperparameters
thres_type={'zPos_Bin'};
sparsity_thres='1.000%';
globalsignal='nGSR';
global_fix='';
thres_p=0.05;
if strcmp(globalsignal,'GSR')
    global_fix='global';
end
DataRoot='H:\all_subjects_res\';
fMRI_brainsize=[61,73,61];
all_nobal_parameters={'bc','cp','deg','ec','eff','modu','pc'};
all_global_parameters={'smallWorld','gamma','lambda','Lp','cp'};
[GMmaskdata,~,~,Header]=rest_to4d('E:\MATLAB toolboxes\SeeCAT\templates\GMwithoutCER_61x73x61.nii');
Header.dt=[16,0];
GMmaskdata=reshape(GMmaskdata,1,[]);
GMmaskdata=logical(GMmaskdata);
GMmaskindexs=find(GMmaskdata==1);
%% load demogrqhic information from EXCEL data, EXCEL sequence needs to be the same in GRETNA runnning folder
[data,txt]=xlsread(strcat(DataRoot,'all_subj_info.xlsx'));
Subject_IDs=txt(2:end,1);
age=data(:,1);
group=data(:,2);
kick_out_IDs={''};
kick_out_indexs=[];
for temp_ID=1:length(kick_out_IDs)
    temp_drop_ind=find(strcmp(Subject_IDs,kick_out_IDs{temp_ID}));
    kick_out_indexs=[kick_out_indexs;temp_drop_ind];
end
group(kick_out_indexs,:)=[];
age(kick_out_indexs,:)=[];
subj_Num=length(Subject_IDs)-length(kick_out_indexs);
%% calculate and save nodal results
% load FCS data
% for temp_thres=1:length(thres_type)
%     for temp_nodal_feature=1:length(all_nobal_parameters)
%         eval(['all_subject_',all_nobal_parameters{temp_nodal_feature},'=zeros(subj_Num,fMRI_brainsize(1)*fMRI_brainsize(2)*fMRI_brainsize(3));'])
%     end
%     
%     for temp_subject=1:subj_Num
%         temp_subject_root=strcat(DataRoot,'ForPAGANI\FunImgARWSD',global_fix,'CFB\unweighted\');
%         %load nodal results
%         temp_data_allpath=spm_select('FPlist',temp_subject_root,strcat(Subject_IDs{temp_subject},'_4DVolume_45381_spa',sparsity_thres,'_cor.*\.nii'));
%         for temp_path_ind=1:size(temp_data_allpath,1)
%             temp_data_path=temp_data_allpath(temp_path_ind,:);
%             for temp_nodal_feature=1:length(all_nobal_parameters)
%                 if ~isempty(findstr(temp_data_path,strcat(all_nobal_parameters{temp_nodal_feature},'.nii')))
%                     temp_subject_data=rest_to4d(temp_data_path);
%                     disp(length(find(temp_subject_data~=0)))
%                     eval(['all_subject_',all_nobal_parameters{temp_nodal_feature},'(temp_subject,:)=reshape(temp_subject_data,1,[]);'])
%                 end
%             end
%         end
%     end
%     
%     %add GM mask, only calculate voxels within GM mask
%     for temp_voxel=1:length(GMmaskindexs)
%         for temp_nodal_feature=1:length(all_nobal_parameters)
%             eval(['temp_diff_stats=gretna_glm(all_subject_',all_nobal_parameters{temp_nodal_feature},'(:,GMmaskindexs(temp_voxel)),[group age],','''t''',',1);']);
%             eval(['all_subject_',all_nobal_parameters{temp_nodal_feature},'_diff_Tmap(GMmaskindexs(temp_voxel))=temp_diff_stats.t;'])
%         end
%     end
%     
%     % reshape and save voxel wise difference result for SeeCAT multiple comparsion correction
%     for temp_index=1:length(all_nobal_parameters)
%         temp_feature=all_nobal_parameters{temp_index};
%         eval(['all_subject_',all_nobal_parameters{temp_nodal_feature},'_diff_Tmap=reshape(all_subject_',all_nobal_parameters{temp_nodal_feature},'_diff_Tmap,fMRI_brainsize);']);
%         temp_diff_Tmap_savename=strcat(DataRoot,'statistical_res\',globalsignal,'\',globalsignal,'_',thres_type,'_Wei_000-Inf_',temp_feature,'_voxel_diff_Tmap_spar',sparsity_thres,'.nii');
%         Header.fname=temp_diff_Tmap_savename;
%         eval(['rest_Write4DNIfTI(all_subject_',all_nobal_parameters{temp_nodal_feature},'_diff_Tmap,Header,temp_diff_Tmap_savename);']);
%     end
% end
%% calculate and save global results
% load FCS data
for temp_thres=1:length(thres_type)
    for temp_global_feature=1:length(all_global_parameters)
        eval(['all_subject_global',all_global_parameters{temp_global_feature},'=zeros(subj_Num,1);'])
    end
    
    for temp_subject=1:subj_Num
        temp_subject_root=strcat(DataRoot,'ForPAGANI\FunImgARWSD',global_fix,'CFB\unweighted\');
        %load global results
        temp_data_allpath=spm_select('FPlist',temp_subject_root,strcat(Subject_IDs{temp_subject},'_4DVolume_45381_spa',sparsity_thres,'_cor.*\.txt'));
        for temp_path_ind=1:size(temp_data_allpath,1)
            temp_data_path=temp_data_allpath(temp_path_ind,:);
            for temp_global_feature=1:length(all_global_parameters)
                if ~isempty(findstr(temp_data_path,strcat(all_global_parameters{temp_global_feature},'.txt')))
                    temp_txt_data=dlmread(temp_data_path);
                    disp(Subject_IDs{temp_subject})
                    disp(temp_data_allpath)
                    eval(['all_subject_global',all_global_parameters{temp_global_feature},'(temp_subject)=temp_txt_data(1);'])
                end
            end
        end
    end
    
     %generate bar plot for all global features 
    bp = BarPlot('ylabel', 'Value');
    cmap = [240,79,0;1,111,176]./255;
    
    
    all_subject_global_diff_T=zeros(1,length(all_global_parameters));
    all_subject_global_diff_P=zeros(1,length(all_global_parameters));
     
    for temp_global_feature=1:length(all_global_parameters)
        temp_feature=all_global_parameters{temp_global_feature};
        eval(['temp_diff_stats=gretna_glm(all_subject_global',temp_feature,',[group age],','''t''',',1);']);
        all_subject_global_diff_T(temp_global_feature)=temp_diff_stats.t;
        all_subject_global_diff_P(temp_global_feature)=temp_diff_stats.p;
        % create the group and the 2 bars within group
        g = bp.addGroup(temp_feature);
        eval(['temp_feature_values=all_subject_global',temp_feature,';']);
        temp_feature_values=temp_feature_values./max(temp_feature_values);
        temp_hc_mean=mean(temp_feature_values(find(group==1)));
        temp_hc_std=std(temp_feature_values(find(group==1)));
        temp_bc_mean=mean(temp_feature_values(find(group==2)));
        temp_bc_std=std(temp_feature_values(find(group==2)));      
        b1 = g.addBar('HC', temp_hc_mean, 'error', temp_hc_std, 'FaceColor', cmap(1, :));
        b2 = g.addBar('BC', temp_bc_mean, 'error', temp_bc_std, 'FaceColor', cmap(2, :));
        % add the bridge connecting bar 1 and bar 2
        if all_subject_global_diff_P(temp_global_feature)<thres_p
            g.addBridge('*', b1, b2, 'FontSize', 12);
        end
        
        
    end
    bp.render();
    temp_savename=strcat(DataRoot,'statistical_res\',globalsignal,'\',globalsignal,'_',thres_type{temp_thres},'_spar',sparsity_thres,'_global_diff_barplot.jpg');
    saveas(gcf,temp_savename);
    close(gcf);
    
    temp_savename=strcat(DataRoot,'statistical_res\',globalsignal,'\',globalsignal,'_',thres_type{temp_thres},'_spar',sparsity_thres,'_global_diff_pvalues.mat');
    save(temp_savename,'all_subject_global_diff_P','all_subject_global_diff_T');
    
  
    
end

