clear all;clc;
%generate statistical correlation pictures figure 5
%hyperparameter
globalsignal='nGSR';
thres_type={'zPos_Bin_'};
freq='001-010';
Distance_Range='20-180';
zscore_fix='z';
behavioral_features={'FCRI','FCRtotal'};
data_root=strcat('H:/breast_cancer/');
ClusterConnectivityCriterion=26;
maskroot=strcat('E:/MATLAB toolboxes/SeeCAT/templates/GMwithoutCER_61x73x61.nii');%only calculate FCS on GM mask
[maskdata,~,~,Header]=rest_to4d(maskroot);
Yeomaskroot=strcat('E:/Matlab Projects/BCT/2016_01_16_BCT/yeo/Yeo_JNeurophysiol11_MNI152/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');%only calculate FCS on GM mask
Yeomaskdata=rest_to4d(Yeomaskroot);
Yeolabelnames={'Visual','Smatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontal Parietal','Default Mode'};
excelpath=strcat(data_root,'breast_cancer_subj_info.xlsx');
[behavioraldata,subjlist]=xlsread(excelpath);
behavioraldata=behavioraldata(:,[3,4,6,7]);%match with behavioral
subjlist=subjlist(2:end,1);%remove variable name
VoxelPThreshold=0.001;
for temp_thres=1:length(thres_type)
    for temp_behavioral_ind=1:length(behavioral_features)
        all_FCScluster_belongings={};
        all_FCcluster_belongings={};
        all_cluster_metrics={};
        all_cluster_PeakR=[];
        all_cluster_metrics_map_to_Yeo={};
        all_cluster_cohensd={};
        %kick out subject with missing behavioral data
        temp_behavioral_data=behavioraldata(:,temp_behavioral_ind);
        all_group=zeros(1,length(temp_behavioral_data));
        switch temp_behavioral_ind
            case 1
                all_group(find(temp_behavioral_data>=13))=1;
            case 2
                all_group(find(behavioraldata(:,1)>=13))=1;
            case 3
                all_group(find(temp_behavioral_data>=10))=1;
            case 4
                all_group(find(temp_behavioral_data>=5))=1;
        end
        temp_subjlist=subjlist(find(temp_behavioral_data~=-1));
        all_group=all_group(find(temp_behavioral_data~=-1));
        FCS_temp_path=strcat(data_root,'/statistical_res/',globalsignal,'/FCS_',thres_type{temp_thres},'corr_with_',behavioral_features{temp_behavioral_ind},'/FCS_significant_corr_with_',behavioral_features{temp_behavioral_ind},'.nii');
        if exist(FCS_temp_path,'file')
            [FCS_data,~,~,Header]=rest_to4d(FCS_temp_path);
            FCS_cluster_belongings=amos_ClusterReport(FCS_data,Header,ClusterConnectivityCriterion,6,1);
            for temp_FCScluster_ind=1:length(FCS_cluster_belongings)
                temp_FCScluster=FCS_cluster_belongings{temp_FCScluster_ind};
                %prepare data for cohens d calculation
                all_subj_maskedvoxelFC=zeros(length(temp_subjlist),length(find(maskdata==1)));
                for temp_subj=1:length(temp_subjlist)
                    if strcmp(globalsignal,'nGSR')
                        globalind='';
                    else
                        globalind='global';
                    end
                    subjPath=strcat(data_root,'/Results/Freq',freq,'_FC_FunImgARWSD',globalind,'CFB_',Distance_Range,'/', behavioral_features{temp_behavioral_ind},'_',zscore_fix,'FC_Seed_',temp_FCScluster,'/',...
                        zscore_fix,'FC_Seed_',temp_FCScluster,'_',temp_subjlist{temp_subj},'.nii');
                    subjAllVolume=rest_to4d(subjPath);
                    subjAllVolume=reshape(subjAllVolume,1,[]);
                    all_subj_maskedvoxelFC(temp_subj,:)=subjAllVolume(find(maskdata==1));
                end
                all_subj_voxelFC=zeros(length(temp_subjlist),length(maskdata));
                all_subj_voxelFC(:,find(maskdata==1))=all_subj_maskedvoxelFC;
                FC_temp_path=strcat(data_root,'/statistical_res/',globalsignal,'/FC_',thres_type{temp_thres},temp_FCScluster,'_corr_with_',behavioral_features{temp_behavioral_ind},'/FC_significant_corr_with_',behavioral_features{temp_behavioral_ind},'.nii');
                FC_data=rest_to4d(FC_temp_path);
                [FC_cluster_belongings,~,FC_cluster_sub,~,~,~]=amos_ClusterReport(FC_data,Header,ClusterConnectivityCriterion,6,1);
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
                    all_cluster_metrics=[all_cluster_metrics;['FC_corr_with_',behavioral_features{temp_behavioral_ind}]];
                    all_cluster_PeakR=[all_cluster_PeakR;max(FC_data(FC_cluster_sub{temp_FCcluster_ind}))];
                    all_cluster_metrics_map_to_Yeo=[all_cluster_metrics_map_to_Yeo;max_overlap_netowrk];
                    % add cohen's d computation
                    temp_cluster_allsubj_voxelFC=all_subj_voxelFC(:,FC_cluster_sub{temp_FCcluster_ind});
                    all_cluster_cohensd=[all_cluster_cohensd;mean(computeCohen_d(temp_cluster_allsubj_voxelFC(find(all_group==1),:),temp_cluster_allsubj_voxelFC(find(all_group==2),:)),2)];
                end
            end
            
            %seperate report by positive and negative FCS clusters, radar map
            %also needs to be grouped by pos and negative clusters, can check
            %network belonging by circos plot
            pos_corr_indexs=find(all_cluster_PeakR>0);
            pos_corr_cluster_networks=all_cluster_metrics_map_to_Yeo(pos_corr_indexs);
            pos_corr_network_map_counts=zeros(length(Yeolabelnames),1);
            
            neg_corr_indexs=find(all_cluster_PeakR<0);
            neg_corr_cluster_networks=all_cluster_metrics_map_to_Yeo(neg_corr_indexs);
            neg_corr_network_map_counts=zeros(length(Yeolabelnames),1);
            for temp_network=1:length(Yeolabelnames)
                isExist = strcmp(pos_corr_cluster_networks, Yeolabelnames{temp_network});
                indexs = find(isExist);
                pos_corr_network_map_counts(temp_network)=length(indexs);
                
                isExist = strcmp(neg_corr_cluster_networks, Yeolabelnames{temp_network});
                indexs = find(isExist);
                neg_corr_network_map_counts(temp_network)=length(indexs);
            end
            %generate report table at the same time
            variablenames={'Number','Region','PeakR','Metric','NetworkBelongings','FCScluster','Cohensd'};
            indexs=1:length(all_FCcluster_belongings);
            indexs=indexs';
            finaltab=table(indexs,all_FCcluster_belongings,all_cluster_PeakR,all_cluster_metrics,all_cluster_metrics_map_to_Yeo,all_FCScluster_belongings,all_cluster_cohensd,'VariableNames',variablenames);
            savename=strcat(data_root,'statistical_res/',globalsignal,'/Table4_','Freq',freq,zscore_fix,Distance_Range,'_',behavioral_features{temp_behavioral_ind},'_All_FC_significant_correlation_cluster_table.xlsx');
            writetable(finaltab,savename);
            
            test_break=1;
        end
    end
end





