clear all;clc;
%generate FC difference map for SeeCAT thresholding
%hyperparameter
zscore_fix={'z'};
Distance_Range={'20-180'};
Network_fix='Pos_Wei_';
all_freq={'001-010'};
globalsig='GSR\';
LoadRoot=strcat('D:\BaiduNetdiskDownload\breast_cancer_T1_fMRI_40\');
excelroot='D:\BaiduNetdiskDownload\breast_cancer_T1_fMRI_40\breast_cancer_subj_info.xlsx';
[behavioraldata,subjlist]=xlsread(excelroot);
behavioraldata=behavioraldata(:,[3,4,6,7]);%do not consider firstep/medicate effect
behavioral={'FCRI','FCRtotal','GAD','PHQ'};
subjlist=subjlist(2:length(behavioraldata)+1);%remove variable name
maskroot=strcat('E:/Matlab toolboxes/SeeCAT/templates/GMwithoutCER_61x73x61.nii');%only calculate FCS on GM mask
maskdata=rest_to4d(maskroot);
[X,Y,Z]=size(maskdata);
maskdata=reshape(maskdata,1,[]);
maskdata=logical(maskdata);
ClusterConnectivityCriterion=26;
cluster_size_thres=0;%do not use additional threshold
Header = spm_vol(maskroot);
Header.dt(1) = 16;
for temp_behavioral=1:length(behavioral)
    group=zeros(size(behavioraldata,1),1);
    switch temp_behavioral
        case 1
            group(find(behavioraldata(:,1)<=13))=1;
            group(find(behavioraldata(:,1)>13))=2;
        case 2
            group(find(behavioraldata(:,1)<=13))=1;
            group(find(behavioraldata(:,1)>13))=2;
        case 3
            group(find(behavioraldata(:,3)<=4))=1;
            group(find(behavioraldata(:,3)>4))=2;
        case 4
            group(find(behavioraldata(:,4)<=10))=1;
            group(find(behavioraldata(:,4)>10))=2;            
    end;
    age=behavioraldata(:,1);
    mod=[group age];
    for tempfix=1:length(zscore_fix)
        for tempRange=1:length(Distance_Range)
            tempPath=strcat(LoadRoot,'statistical_res\',globalsig,'\FCS_',Network_fix(1:3),'_corr_with_',behavioral{temp_behavioral},'\FCS_significant_corr_with_',behavioral{temp_behavioral},'.nii');%only calculate FCS on ANCOVA mask
            data=rest_to4d(tempPath);
            [cluster_belongings,cluster_peakintensity,cluster_sub,cluster_peakpos,cluster_peakcoords,cluster_sizes]=amos_ClusterReport(data,Header,ClusterConnectivityCriterion);
            cluster_sizes=cell2mat(cluster_sizes);
            large_cluster=find(cluster_sizes>cluster_size_thres);%only report large clusters
            for temp_freq=1:length(all_freq)
                %change to seecat processed data
                all_voxelFC=zeros(length(subjlist),length(find(maskdata==1)),length(large_cluster));
                for temp_subj=1:length(subjlist)
                    for tempcluster=1:length(large_cluster)
                        tempclustername=cluster_belongings{tempcluster};
                        if strcmp(globalsig,'nGSR\')
                            globalind='';
                        else
                            globalind='global';
                        end
                        subjPath=strcat(LoadRoot,'\Results\Freq',all_freq{temp_freq},'_FC_FunImgARWSD',globalind,'CFB_',Distance_Range{tempRange},'\', behavioral{temp_behavioral},'_',zscore_fix{tempfix},'FC_Seed_',tempclustername,'\',...
                            zscore_fix{tempfix},'FC_Seed_',tempclustername,'_',subjlist{temp_subj},'.nii');
                        subjAllVolume=rest_to4d(subjPath);
                        subjAllVolume=reshape(subjAllVolume,1,[]);
                        all_voxelFC(temp_subj,:,tempcluster)=subjAllVolume(find(maskdata==1));
                    end
                end
                
                for tempcluster=1:length(large_cluster)                  
                    %save cluster level results
                    temp_clustername=cluster_belongings{large_cluster(tempcluster)};
                    %diff map
                    %obtain mdd vs nc FC difference map for current cluster,regress age and sex
                    stats = gretna_glm(all_voxelFC(:,:,tempcluster),mod,'t',1);%match prof Xia's definition 
                    
                    Group_diff_tmap=zeros(1,length(maskdata));
                    Group_diff_tmap(find(maskdata==1))=stats.t;
                    Group_diff_tmap(find(isnan(Group_diff_tmap)))=0;
                    Group_diff_tmap=reshape(Group_diff_tmap,[X,Y,Z]);
                    Header.fname=strcat(LoadRoot,'statistical_res\',globalsig,'Freq',all_freq{temp_freq},'_',zscore_fix{tempfix},Network_fix,Distance_Range{tempRange},'_',behavioral{temp_behavioral},'_Cluster_',temp_clustername,'_Group_FCdiff_Tmap.nii');
                    spm_write_vol(Header,Group_diff_tmap);
                                      
                end
            end
        end
    end
end