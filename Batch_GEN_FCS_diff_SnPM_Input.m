clear all;clc;
%generate FC difference map for SeeCAT thresholding
%hyperparameter
zscore_fix={'','z'};
Distance_Range={'000-Inf'};
thres_type={'Pos_Wei_','Pos_Bin_'};
globalsig='nGSR';
LoadRoot=strcat('H:\all_subjects_res\');
excelroot='H:\all_subjects_res\all_subj_info_new.xlsx';
[behavioraldata,subjlist]=xlsread(excelroot);
group=behavioraldata(:,2);%do not consider firstep/medicate effect
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

for temp_thres=1:length(thres_type)
    for tempfix=1:length(zscore_fix)
        for tempRange=1:length(Distance_Range)
            for temp_subj=1:length(subjlist)
                if strcmp(globalsig,'nGSR')
                    globalind='';
                else
                    globalind='global';
                end
                subjPath=strcat(LoadRoot,'\Results\','VoxelDeg_FunImgARWSD',globalind,'CFB\', zscore_fix{tempfix},thres_type{temp_thres},Distance_Range{tempRange},'\');
                subjFileName=strcat(zscore_fix{tempfix},thres_type{temp_thres},Distance_Range{tempRange},'_',subjlist{temp_subj},'.nii');
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
                end
            end
        end
    end
end
