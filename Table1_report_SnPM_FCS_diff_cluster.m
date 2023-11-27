clear all;clc;
%generate statistical correlation pictures figure 5
%hyperparameter
globalsignal='nGSR';
zscore_fix={'z'};
thres_type={'Pos_Bin'};
group_fix={''};
data_root='H:/all_subjects_res/';
ClusterConnectivityCriterion=26;
maskroot=strcat('E:/MATLAB toolboxes/SeeCAT/templates/GMwithoutCER_61x73x61.nii');%only calculate FCS on GM mask
[maskdata,~,~,Header]=rest_to4d(maskroot);
Yeomaskroot=strcat('E:/Matlab Projects/BCT/2016_01_16_BCT/yeo/Yeo_JNeurophysiol11_MNI152/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');%only calculate FCS on GM mask
Yeomaskdata=rest_to4d(Yeomaskroot);
Yeolabelnames={'Visual','Smatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontal Parietal','Default Mode'};
excelpath=strcat(data_root,'all_subj_info.xlsx');
[behavioraldata,subjlist]=xlsread(excelpath);
subjlist= subjlist(2:end,:);
group=behavioraldata(:,2);%do not consider firstep/medicate effect
medicated=behavioraldata(:,5);
VoxelPThreshold=0.001;
IsTwoTailed=true;
dLh=0;
for temp_group=1:length(group_fix)
    for temp_thres=1:length(thres_type)
        for temp_fix=1:length(zscore_fix)
            all_FCScluster_belongings={};
            all_cluster_metrics={};
            all_cluster_peakcoord=[];
            all_cluster_size=[];
            all_cluster_PeakR=[];
            all_cluster_p=[];
            all_cluster_metrics_map_to_Yeo={};
            all_cluster_cohensd={};
            temp_path=strcat(data_root,'statistical_res/',globalsignal,'/FCS_',zscore_fix{temp_fix},thres_type{temp_thres},'_group_diff',group_fix{temp_group},'/',globalsignal,'_FCS_significant_group_diff_SnPM.nii');
            VoxelSize=[];
            if logical(exist(temp_path,'file'))
                [data,tempVoxelSize,~,Header]=rest_to4d(temp_path);
                VoxelSize=tempVoxelSize;
                [cluster_belongings,~,cluster_sub,~,cluster_peakcoords,cluster_sizes]=amos_ClusterReport(data,Header,ClusterConnectivityCriterion,6,1);
                cluster_sizes=cell2mat(cluster_sizes);
                for temp_cluster=1:length(cluster_belongings)
                    %needs to replace with snpm p values, check snpm !!!
                    ClusterData = data(cluster_sub{temp_cluster});
                    [ClusterP,~] = GRF_get_cluster_P_and_Mean(maskdata,ClusterData,VoxelPThreshold,dLh,IsTwoTailed);
                    
                    all_subj_maskedvoxelFCS=zeros(length(subjlist),length(find(maskdata==1)));
                    for temp_subj=1:length(subjlist)
                        if strcmp(globalsignal,'nGSR')
                            globalind='';
                        else
                            globalind='global';
                        end
                        subjPath=strcat(data_root,'/Results/VoxelDeg_FunImgARWSD',globalind,'CFB/', zscore_fix{temp_fix},thres_type{temp_thres},'_000-Inf/',zscore_fix{temp_fix},thres_type{temp_thres},'_000-Inf_',subjlist{temp_subj},'.nii');
                        subjAllVolume=rest_to4d(subjPath);
                        subjAllVolume=reshape(subjAllVolume,1,[]);
                        all_subj_maskedvoxelFCS(temp_subj,:)=subjAllVolume(find(maskdata==1));
                    end
                    all_subj_voxelFCS=zeros(length(subjlist),length(maskdata));
                    all_subj_voxelFCS(:,find(maskdata==1))=all_subj_maskedvoxelFCS;               
                    %add Yeo atlas parcellation on report
                    max_overlap_netowrk='None';
                    max_overlap_count=0;
                    for tempnetwork=1:max(unique(Yeomaskdata))
                        tempnetworkind=find(Yeomaskdata==tempnetwork);
                        temp_overlap_count=length(intersect(cluster_sub{temp_cluster},tempnetworkind));
                        if temp_overlap_count>max_overlap_count
                            max_overlap_count=temp_overlap_count;
                            max_overlap_netowrk=Yeolabelnames{tempnetwork};
                        end
                    end
                    all_FCScluster_belongings=[all_FCScluster_belongings;cluster_belongings{temp_cluster}];
                    all_cluster_metrics=[all_cluster_metrics;['FCS_groupdiff_',zscore_fix{temp_fix},thres_type{temp_thres}]];
                    all_cluster_peakcoord=[all_cluster_peakcoord;cell2mat(cluster_peakcoords(temp_cluster)')];
                    if mean(data(cluster_sub{temp_cluster}))>0
                        all_cluster_PeakR=[all_cluster_PeakR;max(data(cluster_sub{temp_cluster}))];
                    else
                        all_cluster_PeakR=[all_cluster_PeakR;min(data(cluster_sub{temp_cluster}))];
                    end
                    all_cluster_size=[all_cluster_size;cluster_sizes(temp_cluster)'];
                    all_cluster_metrics_map_to_Yeo=[all_cluster_metrics_map_to_Yeo;max_overlap_netowrk];
                    all_cluster_p=[all_cluster_p;ClusterP];
                    % add cohen's d computation
                    temp_cluster_allsubj_voxelFCS=all_subj_voxelFCS(:,cluster_sub{temp_cluster});
                    all_cluster_cohensd=[all_cluster_cohensd;mean(computeCohen_d(temp_cluster_allsubj_voxelFCS(find(group==1),:),temp_cluster_allsubj_voxelFCS(find(group==2 & medicated==0),:)),2)];
                end
            end
            % arrange data to table format
            all_cluster_size=all_cluster_size*prod(VoxelSize);
            variablenames={'Number', 'Region', 'Metric', 'x', 'y', 'z', 'PeakR', 'Size_mm', 'Pvalue','NetworkBelongings','Cohensd'};
            indexs=1:length(all_FCScluster_belongings);
            indexs=indexs';
            finaltab=table(indexs,all_FCScluster_belongings,all_cluster_metrics,all_cluster_peakcoord(:,1),...
                all_cluster_peakcoord(:,2),all_cluster_peakcoord(:,3),all_cluster_PeakR,all_cluster_size,all_cluster_p,all_cluster_metrics_map_to_Yeo,all_cluster_cohensd,'VariableNames',variablenames);
            savename=strcat(data_root,'statistical_res/',globalsignal,'/Table1_FCS_',zscore_fix{temp_fix},thres_type{temp_thres},group_fix{temp_group},'_significant_groupdiff_behavioral_cluster_info.xlsx');
            writetable(finaltab,savename);
        end
    end
end
