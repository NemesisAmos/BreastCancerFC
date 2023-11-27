clear all;clc;
%generate statistical correlation pictures figure 5
%hyperparameter
globalsignal='GSR';
thres_type={'zPos_Bin_'};
behavioral_features={'FCRI','FCRtotal'};
data_root=strcat('H:/breast_cancer/statistical_res/',globalsignal,'/');
ClusterConnectivityCriterion=26;
maskroot=strcat('E:/MATLAB toolboxes/SeeCAT/templates/GMwithoutCER_61x73x61.nii');%only calculate FCS on GM mask
[maskdata,~,~,Header]=rest_to4d(maskroot);
Yeomaskroot=strcat('E:/Matlab Projects/BCT/2016_01_16_BCT/yeo/Yeo_JNeurophysiol11_MNI152/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');%only calculate FCS on GM mask
Yeomaskdata=rest_to4d(Yeomaskroot);
Yeolabelnames={'Visual','Smatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontal Parietal','Default Mode'};
all_FCScluster_belongings={};
all_cluster_peakintensity=[];
all_cluster_metrics={};
all_cluster_peakcoord=[];
all_cluster_size=[];
all_cluster_meanR=[];
all_cluster_p=[];
all_cluster_metrics_map_to_Yeo={};
VoxelPThreshold=0.001;
IsTwoTailed=true;
dLh=0;
for temp_thres=1:length(thres_type)
    for temp_behavioral=1:length(behavioral_features)
        temp_path=strcat(data_root,'FCS_',thres_type{temp_thres},'corr_with_',behavioral_features{temp_behavioral},'/FCS_significant_corr_with_',behavioral_features{temp_behavioral},'.nii');
        if exist(temp_path,'file')
            [data,VoxelSize,~,Header]=rest_to4d(temp_path);
            [cluster_belongings,cluster_peakintensity,cluster_sub,~,cluster_peakcoords,cluster_sizes]=amos_ClusterReport(data,Header,ClusterConnectivityCriterion,6,1);
            cluster_sizes=cell2mat(cluster_sizes);
            for temp_cluster=1:length(cluster_belongings)
                %needs to replace with snpm p values, check snpm !!!
                ClusterData = data(cluster_sub{temp_cluster});
                [ClusterP,~] = GRF_get_cluster_P_and_Mean(maskdata,ClusterData,VoxelPThreshold,dLh,IsTwoTailed);
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
                all_cluster_peakintensity=[all_cluster_peakintensity;cluster_peakintensity(temp_cluster)];
                all_cluster_metrics=[all_cluster_metrics;['FCS_corr_with_',behavioral_features{temp_behavioral}]];
                all_cluster_peakcoord=[all_cluster_peakcoord;cell2mat(cluster_peakcoords(temp_cluster)')];
                all_cluster_meanR=[all_cluster_meanR;mean(data(cluster_sub{temp_cluster}))];
                all_cluster_size=[all_cluster_size;cluster_sizes(temp_cluster)'];
                all_cluster_metrics_map_to_Yeo=[all_cluster_metrics_map_to_Yeo;max_overlap_netowrk];
                all_cluster_p=[all_cluster_p;ClusterP];
            end
        end
    end
end
% arrange data to table format
all_cluster_size=all_cluster_size*prod(VoxelSize);
variablenames={'Number', 'Region', 'Metric', 'x', 'y', 'z', 'Peak_R', 'Size_mm', 'Pvalue','NetworkBelongings'};
indexs=1:length(all_FCScluster_belongings);
indexs=indexs';
finaltab=table(indexs,all_FCScluster_belongings,all_cluster_metrics,all_cluster_peakcoord(:,1),...
    all_cluster_peakcoord(:,2),all_cluster_peakcoord(:,3),all_cluster_peakintensity,all_cluster_size,all_cluster_p,all_cluster_metrics_map_to_Yeo,'VariableNames',variablenames);
savename=strcat(data_root,'Table2_FCS_',thres_type{temp_thres},'significant_corr_with_behavioral_cluster_info.xlsx');
writetable(finaltab,savename);