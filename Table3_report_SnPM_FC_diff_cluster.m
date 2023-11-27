clear all;clc;
%generate statistical correlation pictures figure 5
%hyperparameter
globalsignal='nGSR';
zscore_fix={'z'};
thres_type={'Pos_Wei_'};
group_type={''};
freq='001-010';
Distance_Range='20-180';
data_root=strcat('H:/all_subjects_res/');
ClusterConnectivityCriterion=26;
maskroot=strcat('E:/MATLAB toolboxes/SeeCAT/templates/GMwithoutCER_61x73x61.nii');%only calculate FCS on GM mask
[maskdata,~,~,Header]=rest_to4d(maskroot);
Yeomaskroot=strcat('E:/Matlab Projects/BCT/2016_01_16_BCT/yeo/Yeo_JNeurophysiol11_MNI152/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');%only calculate FCS on GM mask
Yeomaskdata=rest_to4d(Yeomaskroot);
Yeolabelnames={'Visual','Smatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontal Parietal','Default Mode'};
excelpath=strcat(data_root,'all_subj_info.xlsx');
[behavioraldata,subjlist]=xlsread(excelpath);
group=behavioraldata(:,2);%do not consider firstep/medicate effect
age=behavioraldata(:,1);%do not consider firstep/medicate effect
subjlist=subjlist(2:end,1);%remove variable name

VoxelPThreshold=0.001;
for temp_group=1:length(group_type)
    for temp_thres=1:length(thres_type)
        for temp_fix=1:length(zscore_fix)
            all_FCScluster_belongings={};
            all_FCcluster_belongings={};
            all_cluster_peakcoord=[];
            all_cluster_size=[];
            all_cluster_PeakT=[];
            all_cluster_metrics_map_to_Yeo={};
            all_cluster_cohensd={};
            %kick out subject with missing behavioral data
            FCS_temp_path=strcat(data_root,'/statistical_res/',globalsignal,'/FCS_',zscore_fix{temp_fix},thres_type{temp_thres},'group_diff',group_type{temp_group},'/',globalsignal,'_FCS_significant_group_diff_SnPM.nii');
            [FCS_data,tempVoxelSize,~,Header]=rest_to4d(FCS_temp_path);
            VoxelSize=tempVoxelSize;
            FCS_cluster_belongings=amos_ClusterReport(FCS_data,Header,ClusterConnectivityCriterion,6,1);
            for temp_FCScluster_ind=1:length(FCS_cluster_belongings)
                temp_FCScluster=FCS_cluster_belongings{temp_FCScluster_ind};
                %prepare data for cohens d calculation
                all_subj_maskedvoxelFC=zeros(length(subjlist),length(find(maskdata==1)));
                for temp_subj=1:length(subjlist)
                    if strcmp(globalsignal,'nGSR')
                        globalind='';
                    else
                        globalind='global';
                    end
                    subjPath=strcat(data_root,'/Results/Freq',freq,'_FC_FunImgARWSD',globalind,'CFB_',Distance_Range,'/', zscore_fix{temp_fix},thres_type{temp_thres},'_FC_Seed_',temp_FCScluster,'/FC_Seed_',temp_FCScluster,'_',subjlist{temp_subj},'.nii');
                    subjAllVolume=rest_to4d(subjPath);
                    subjAllVolume=reshape(subjAllVolume,1,[]);
                    all_subj_maskedvoxelFC(temp_subj,:)=subjAllVolume(find(maskdata==1));
                end
                all_subj_voxelFC=zeros(length(subjlist),length(maskdata));
                all_subj_voxelFC(:,find(maskdata==1))=all_subj_maskedvoxelFC;
                FC_temp_path=strcat(data_root,'/statistical_res/',globalsignal,'/',zscore_fix{temp_fix},thres_type{temp_thres},'_',zscore_fix{temp_fix},'FC_Seed_',temp_FCScluster,group_type{temp_group},'/FC_significant_group_diff_SnPM.nii');
                FC_data=rest_to4d(FC_temp_path);
                [FC_cluster_belongings,~,FC_cluster_sub,FC_cluster_peakpos,FC_cluster_peakcoords,FC_cluster_sizes]=amos_ClusterReport(FC_data,Header,ClusterConnectivityCriterion,6,1);
                for temp_FCcluster_ind=1:length(FC_cluster_belongings)
                    max_overlap_netowrk='None';
                    max_overlap_count=0;
                    for tempnetwork=1:max(unique(Yeomaskdata))
                        tempnetworkind=find(Yeomaskdata==tempnetwork);
                        temp_overlap_count=length(intersect(FC_cluster_sub{temp_FCcluster_ind},tempnetworkind));
                        if temp_overlap_count>max_overlap_count
                            max_overlap_count=temp_overlap_count;
                            max_overlap_netowrk=Yeolabelnames{tempnetwork};
                        end
                    end
                    all_FCScluster_belongings=[all_FCScluster_belongings;temp_FCScluster];
                    all_FCcluster_belongings=[all_FCcluster_belongings;FC_cluster_belongings{temp_FCcluster_ind}];
                    all_cluster_peakcoord=[all_cluster_peakcoord;cell2mat(FC_cluster_peakcoords(temp_FCcluster_ind)')];
                    all_cluster_size=[all_cluster_size;FC_cluster_sizes(temp_FCcluster_ind)'];
                    if mean(FC_data(FC_cluster_sub{temp_FCcluster_ind}))>0
                        all_cluster_PeakT=[all_cluster_PeakT;max(FC_data(FC_cluster_sub{temp_FCcluster_ind}))];
                    else
                        all_cluster_PeakT=[all_cluster_PeakT;min(FC_data(FC_cluster_sub{temp_FCcluster_ind}))];
                    end
                    all_cluster_metrics_map_to_Yeo=[all_cluster_metrics_map_to_Yeo;max_overlap_netowrk];
                    % add cohen's d computation
                    temp_cluster_allsubj_voxelFC=all_subj_voxelFC(:,FC_cluster_sub{temp_FCcluster_ind});
                    all_cluster_cohensd=[all_cluster_cohensd;mean(computeCohen_d(temp_cluster_allsubj_voxelFC(find(group==1),:),temp_cluster_allsubj_voxelFC(find(group==2),:)),2)];
                end
                test_break=1;
            end
            %generate report table at the same time
            all_cluster_size=cell2mat(all_cluster_size)*prod(VoxelSize);
            variablenames={'Number','PeakT','x', 'y', 'z', 'Region', 'Size_mm', 'NetworkBelongings','FCScluster','Cohensd'};
            indexs=1:length(all_FCcluster_belongings);
            indexs=indexs';
            finaltab=table(indexs,all_cluster_PeakT,all_cluster_peakcoord(:,1),all_cluster_peakcoord(:,2),all_cluster_peakcoord(:,3),all_FCcluster_belongings,all_cluster_size,all_cluster_metrics_map_to_Yeo,all_FCScluster_belongings,all_cluster_cohensd,'VariableNames',variablenames);
            savename=strcat(data_root,'statistical_res/',globalsignal,'/Tabel3_Freq',freq,'_',zscore_fix{temp_fix},thres_type{temp_thres},Distance_Range,group_type{temp_group},'_All_FC_significant_group_diff_cluster_table.xlsx');
            writetable(finaltab,savename);
        end
    end
end


