clear all;clc;
%generate FC difference map for SeeCAT thresholding
%hyperparameter
zscore_fix={'z'};
Distance_Range={'20-180'};
thres_type={'Pos_Bin_'};
group_type={''};
all_freq={'001-010'};
globalsig='nGSR';
LoadRoot=strcat('H:\all_subjects_res\');
excelroot='H:\all_subjects_res\all_subj_info_new.xlsx';
[behavioraldata,subjlist]=xlsread(excelroot);
group=behavioraldata(:,4);%do not consider firstep/medicate effect
subjlist=subjlist(2:length(behavioraldata)+1);%remove variable name
maskroot=strcat('E:\Matlab toolboxes\SeeCAT\templates\GMwithoutCER_61x73x61.nii');%only calculate FCS on GM mask
maskdata=rest_to4d(maskroot);
[X,Y,Z]=size(maskdata);
maskdata=reshape(maskdata,1,[]);
maskdata=logical(maskdata);
ClusterConnectivityCriterion=26;
cluster_size_thres=0;%do not use additional threshold
Header = spm_vol(maskroot);
Header.dt(1) = 16;
for temp_group=1:length(group_type)
    for temp_thres=1:length(thres_type)
        for tempfix=1:length(zscore_fix)
            for tempRange=1:length(Distance_Range)
                tempPath=strcat(LoadRoot,'statistical_res\',globalsig,'\FCS_',zscore_fix{tempfix},thres_type{temp_thres},'group_diff',group_type{temp_group},'\',globalsig,'_FCS_significant_group_diff_SnPM.nii');%only calculate FCS on ANCOVA mask
                data=rest_to4d(tempPath);
                [cluster_belongings,cluster_peakintensity,cluster_sub,cluster_peakpos,cluster_peakcoords,cluster_sizes]=amos_ClusterReport(data,Header,ClusterConnectivityCriterion,6,1);
                cluster_sizes=cell2mat(cluster_sizes);
                large_cluster=find(cluster_sizes>cluster_size_thres);%only report large clusters
                for temp_freq=1:length(all_freq)
                    %change to seecat processed data
                    all_voxelFC=zeros(length(subjlist),length(find(maskdata==1)),length(large_cluster));
                    for temp_subj=1:length(subjlist)
                        for tempcluster=1:length(large_cluster)
                            tempclustername=cluster_belongings{tempcluster};
                            if strcmp(globalsig,'nGSR')
                                globalind='';
                            else
                                globalind='global';
                            end
                            subjPath=strcat(LoadRoot,'\Results\Freq',all_freq{temp_freq},'_FC_FunImgARWSD',globalind,'CFB_',Distance_Range{tempRange},'\', zscore_fix{tempfix},thres_type{temp_thres},'_',zscore_fix{tempfix},'FC_Seed_',tempclustername,'\');
                            subjFileName=strcat(zscore_fix{tempfix},'FC_Seed_',tempclustername,'_',subjlist{temp_subj},'.nii');
                            %move file to group folder based on group index
                            if group(temp_subj)==1
                                destDir=strcat(subjPath,'\Group1');
                            else
                                destDir=strcat(subjPath,'\Group2');
                            end
                            
                            if ~exist(destDir,'dir')
                                mkdir(destDir);
                            end
                            tempFileName=strcat(subjPath,subjFileName);
                            if exist(tempFileName,'file')
                                movefile(tempFileName,strcat(destDir));
                            else
                                disp(['File ',tempFileName,' not found....'])
                            end
                        end
                    end
                end
            end
        end
    end
end