clear all;clc;
%generate statistical correlation pictures figure 5
%hyperparameter
globalsignal='nGSR';
thres_type={'zPos_Bin_'};
group_type={''};
freq='001-010';
Distance_Range='20-180';
behavioral_features={'FCRI','FCRtotal'};
data_root=strcat('H:/breast_cancer/');
ClusterConnectivityCriterion=26;
%load excel data
excelpath=strcat(data_root,'breast_cancer_subj_info.xlsx');
[behavioraldata,subjlist]=xlsread(excelpath);
age=behavioraldata(:,1);
medicated=behavioraldata(:,5);
behavioraldata=behavioraldata(:,[3,4,6,7]);%match with behavioral
subjlist=subjlist(2:end,1);%remove variable name
%load subject signficant cluster peak FCS data
for temp_group=1:length(group_type)
    for temp_thres=1:length(thres_type)
        for temp_behavioral_ind=1:length(behavioral_features)
            temp_behavioral_data=behavioraldata(:,temp_behavioral_ind);
            temp_subjlist=subjlist(find(temp_behavioral_data~=-1 & medicated~=1));
            temp_behavioral_data=temp_behavioral_data(find(temp_behavioral_data~=-1 & medicated~=1));
            all_cluster_all_subject_cluster_avg_FCS_corr=[];
            all_cluster_all_subject_cluster_peak_FCS=[];
            FCS_temp_path=strcat(data_root,'/statistical_res/',globalsignal,'/FCS_',thres_type{temp_thres},'corr_with_',behavioral_features{temp_behavioral_ind},group_type{temp_group},'/FCS_significant_corr_with_',behavioral_features{temp_behavioral_ind},'.nii');
            if exist(FCS_temp_path,'file')
                [FCS_data,~,~,Header]=rest_to4d(FCS_temp_path);
                [FCS_cluster_belongings,FCS_cluster_peakintensity,FCS_cluster_sub,FCS_cluster_peakpos]=amos_ClusterReport(FCS_data,Header,ClusterConnectivityCriterion,6,1);
                for temp_FCScluster_ind=1:length(FCS_cluster_belongings)
                    temp_FCScluster=FCS_cluster_belongings{temp_FCScluster_ind};
                    all_subj_cluster_peak_FCS=zeros(length(temp_subjlist),1);
                    all_subj_cluster_avg_FCS=zeros(length(temp_subjlist),1);
                    for temp_subject=1:length(temp_subjlist)
                        if strcmp(globalsignal,'nGSR')
                            globalind='';
                        else
                            globalind='global';
                        end
                        subjPath=strcat(data_root,'/Results/VoxelDeg_FunImgARWSD',globalind,'CFB/',thres_type{temp_thres}, '000-Inf/',thres_type{temp_thres},'000-Inf_',temp_subjlist{temp_subject},'.nii');
                        subjAllVolume=rest_to4d(subjPath);
                        subjAllVolume=reshape(subjAllVolume,1,[]);
                        temp_cluster_inds=FCS_cluster_sub{temp_FCScluster_ind};
                        all_subj_cluster_peak_FCS(temp_subject)=subjAllVolume(temp_cluster_inds(FCS_cluster_peakpos{temp_FCScluster_ind}));
                    end
                    all_cluster_all_subject_cluster_peak_FCS=[all_cluster_all_subject_cluster_peak_FCS all_subj_cluster_peak_FCS];
                    temp_cluster_inds=FCS_cluster_sub{temp_FCScluster_ind};
                    all_subj_cluster_avg_FCS(temp_subject)=mean(FCS_data(temp_cluster_inds));
                    all_cluster_all_subject_cluster_avg_FCS_corr=[all_cluster_all_subject_cluster_avg_FCS_corr all_subj_cluster_avg_FCS];
                end
                save(strcat(data_root,'/statistical_res/',globalsignal,'/FCS_correlate_with',behavioral_features{temp_behavioral_ind},group_type{temp_group},'_cluster_avg_FCS.mat'),'FCS_cluster_belongings','all_cluster_all_subject_cluster_peak_FCS',...
                    'all_cluster_all_subject_cluster_avg_FCS_corr');
                %scatter plot FCS with behavioral, needs to plot multiple cluster data into same scatter plot
                tmp_X=temp_behavioral_data;
                rtmp_X=amos_regress(tmp_X,age(find(medicated~=1)));
                tmp_Y=all_cluster_all_subject_cluster_peak_FCS;
                rtmp_Y=amos_regress(tmp_Y,age(find(medicated~=1)));
                [r_tmp,p_tmp]=corr(rtmp_X,rtmp_Y,'type','pearson');
                all_cluster_corr_R = r_tmp;
                sorted_p_tmp=sort(p_tmp);
                [~,rank]=ismember(p_tmp,sorted_p_tmp);
                p_tmp = p_tmp*size(all_cluster_all_subject_cluster_peak_FCS,2)./rank;%FDR correction
                all_cluster_corr_P = p_tmp;    
                for temp_cluster=1:length(FCS_cluster_belongings)
                    tmp_belongingnames = FCS_cluster_belongings{temp_cluster};
                    tmp_belongingnames = regexprep(tmp_belongingnames, '_', ' ');%replace _ with space
                    FCS_cluster_belongings{temp_cluster} = tmp_belongingnames;
                    tmp_X=temp_behavioral_data;
                    rtmp_X=amos_regress(tmp_X,age(find(medicated~=1)));
                    tmp_Y=all_cluster_all_subject_cluster_peak_FCS(:,temp_cluster);
                    rtmp_Y=amos_regress(tmp_Y,age(find(medicated~=1)));
                    figure;
                    gretna_plot_regression(rtmp_X,rtmp_Y);
                    xlabel(behavioral_features{temp_behavioral_ind},'FontWeight','Bold');
                    ylabel(tmp_belongingnames,'FontWeight','Bold');
                    saveas(gcf,strcat(data_root,'/statistical_res/',globalsignal,'/FCS_',thres_type{temp_thres},tmp_belongingnames,'_correlate_with',behavioral_features{temp_behavioral_ind},group_type{temp_group},...
                        'R_',num2str(all_cluster_corr_R(temp_cluster)),'_P_',num2str(all_cluster_corr_P(temp_cluster)),'.jpg'));
                    close(gcf);
                end               
            end
        end
    end
end
