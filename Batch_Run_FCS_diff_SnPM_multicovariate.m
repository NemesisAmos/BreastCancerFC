clear all;clc;
%generate FC difference map for SeeCAT thresholding
%hyperparameter
zscore_fix={'z'};
Distance_Range={'000-Inf'};
thres_type={'Pos_Bin_'};
group_type={''};
globalsig='GSR';
LoadRoot=strcat('H:\all_subjects_res\');
excelroot='H:\all_subjects_res\all_subj_info_new.xlsx';
[behavioraldata,subjlist]=xlsread(excelroot);
age=behavioraldata(:,1);
edu=behavioraldata(:,2);
group=behavioraldata(:,4);%do not consider firstep/medicate effect
FCRI=behavioraldata(:,5);
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
                if strcmp(globalsig,'nGSR')
                    globalind='';
                else
                    globalind='global';
                end
                FCSDataPath=strcat(LoadRoot,'\Results\','VoxelDeg_FunImgARWSD',globalind,'CFB\', zscore_fix{tempfix},thres_type{temp_thres},Distance_Range{tempRange},'\');
                %construct data
                %% SnPM design matrix construction
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.DesignName = '2 Groups: Two Sample T test; 1 scan per subject';
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.DesignFile = 'snpm_bch_ui_TwoSampT';
                output_Dir=strcat(LoadRoot,'statistical_res\',globalsig,'\FCS_',zscore_fix{tempfix},thres_type{temp_thres},'group_diff',group_type{temp_group});
                if ~exist(output_Dir,'dir')
                    mkdir(output_Dir);
                end
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.dir = {output_Dir};
                NC_subjects=subjlist(find(group==1));
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans1=cell(length(NC_subjects),1);
                for temp_subject=1:length(NC_subjects)
                    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans1{temp_subject}=strcat(FCSDataPath,'Group1\',zscore_fix{tempfix},thres_type{temp_thres},Distance_Range{tempRange},'_',NC_subjects{temp_subject},'.nii,1');
                end
                if strcmp(group_type{temp_group},'_FCRI_low')
                    BrestCancer_subjects=subjlist(find(group==2 & FCRI<13));
                elseif strcmp(group_type{temp_group},'_FCRI_high')
                    BrestCancer_subjects=subjlist(find(group==2 & FCRI>=13));
                else
                    BrestCancer_subjects=subjlist(find(group==2));
                end
                
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans2=cell(length(BrestCancer_subjects),1);
                for temp_subject=1:length(BrestCancer_subjects)
                    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans2{temp_subject}=strcat(FCSDataPath,'Group2\',zscore_fix{tempfix},thres_type{temp_thres},Distance_Range{tempRange},'_',BrestCancer_subjects{temp_subject},'.nii,1');
                end

                if strcmp(group_type{temp_group},'_FCRI_low')
                    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(1).c = age(find(FCRI==-1 | FCRI<13));
                    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(2).c = edu(find(FCRI==-1 | FCRI<13));
                elseif strcmp(group_type{temp_group},'_FCRI_high')
                    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(1).c = age(find(FCRI==-1 | FCRI>=13));
                    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(2).c = edu(find(FCRI==-1 | FCRI>=13));
                else
                    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(1).c = age;
                    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(2).c = edu;
                end
                        
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(1).cname = 'Age';
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(2).cname = 'Education';
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.nPerm = 10000;
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.vFWHM = [0 0 0];
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.bVolm = 1;
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.ST.ST_later = -1;
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.tm.tm_none = 1;
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.im = 1;
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.em = {'H:\GMwithoutCER_61x73x61.nii,1'};
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalc.g_omit = 1;
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.gmsca.gmsca_no = 1;
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.glonorm = 1;
                %% SnPM Compute
                matlabbatch{2}.spm.tools.snpm.cp.snpmcfg = {strcat(output_Dir,'\SnPMcfg.mat')};
                %% SnPM Inference
                matlabbatch{3}.spm.tools.snpm.inference.SnPMmat = {strcat(output_Dir,'\SnPM.mat')};
                matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = 0.01;
                matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.PthC = 0.05;
                matlabbatch{3}.spm.tools.snpm.inference.Tsign = 1;
                matlabbatch{3}.spm.tools.snpm.inference.WriteFiltImg.name = strcat(globalsig,'_FCS_significant_positive_group_diff_SnPM');
                matlabbatch{3}.spm.tools.snpm.inference.Report = 'MIPtable';
                
                matlabbatch{4}.spm.tools.snpm.inference.SnPMmat = {strcat(output_Dir,'\SnPM.mat')};
                matlabbatch{4}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = 0.01;
                matlabbatch{4}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.PthC = 0.05;
                matlabbatch{4}.spm.tools.snpm.inference.Tsign = -1;
                matlabbatch{4}.spm.tools.snpm.inference.WriteFiltImg.name = strcat(globalsig,'_FCS_significant_negative_group_diff_SnPM');
                matlabbatch{4}.spm.tools.snpm.inference.Report = 'MIPtable';
                
                spm('defaults', 'FMRI'); % configure spm
                spm_jobman('run', matlabbatch);
            end
        end
    end
end

