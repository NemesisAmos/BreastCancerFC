clear all;clc;
%generate FC difference map for SeeCAT thresholding
%hyperparameter
zscore_fix={''};
Distance_Range={'000-Inf'};
thres_type={'zPos_Bin_'};
group_type={''};
behavioral={'FCRI','FCRtotal'};
globalsig='nGSR';
LoadRoot=strcat('H:\breast_cancer\');
excelroot='H:\breast_cancer\breast_cancer_subj_info.xlsx';
[behavioraldata,subjlist]=xlsread(excelroot);
age=behavioraldata(:,1);
group=behavioraldata(:,2);%do not consider firstep/medicate effect
FCRI=behavioraldata(:,3);
FCRtotal=behavioraldata(:,4);
medicated=behavioraldata(:,5);
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
    for temp_behavioral=1:length(behavioral)
        for temp_thres=1:length(thres_type)
            for tempfix=1:length(zscore_fix)
                for tempRange=1:length(Distance_Range)
                    if strcmp(globalsig,'nGSR')
                        globalind='';
                    else
                        globalind='global';
                    end
                    FCSDataPath=strcat(LoadRoot,'\Results\VoxelDeg_FunImgARWSD',globalind,'CFB\', zscore_fix{tempfix},thres_type{temp_thres},Distance_Range{tempRange},'\');
                    %construct data
                    %% SnPM design matrix construction
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.DesignName = 'SingleSub: Simple Regression (correlation); single covariate of interest';
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.DesignFile = 'snpm_bch_ui_Corr1S';
                    output_Dir=strcat(LoadRoot,'statistical_res\',globalsig,'\FCS_',zscore_fix{tempfix},thres_type{temp_thres},'corr_with_',behavioral{temp_behavioral},group_type{temp_group});
                    if ~exist(output_Dir,'dir')
                        mkdir(output_Dir);
                    end
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.dir = {output_Dir};
%                     breast_cancer_subjlist=subjlist(find(medicated~=1));
                    breast_cancer_subjlist=subjlist;
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.P=cell(length(breast_cancer_subjlist),1);
                    for temp_subject=1:length(breast_cancer_subjlist)
                        matlabbatch{1}.spm.tools.snpm.des.Corr1S.P{temp_subject}=strcat(FCSDataPath,zscore_fix{tempfix},thres_type{temp_thres},Distance_Range{tempRange},'_',breast_cancer_subjlist{temp_subject},'.nii,1');
                    end    
                    if strcmp(behavioral{temp_behavioral},'FCRI')
                        temp_covint=FCRI;
                    else
                        temp_covint=FCRtotal;
                    end
                                      
%                     if strcmp(group_type{temp_group},'_FCRI_low')
%                         matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov.c = age(find(FCRI<13 & medicated~=1));
%                         matlabbatch{1}.spm.tools.snpm.des.Corr1S.CovInt = temp_covint(find(FCRI<13 & medicated~=1));
%                         matlabbatch{1}.spm.tools.snpm.des.Corr1S.P = matlabbatch{1}.spm.tools.snpm.des.Corr1S.P(find(FCRI<13 & medicated~=1));
%                         matlabbatch{1}.spm.tools.snpm.des.Corr1S.xblock = length(find(FCRI<13 & medicated~=1));
%                     elseif strcmp(group_type{temp_group},'_FCRI_high')
%                         matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov.c = age(find(FCRI>=13 & medicated~=1));
%                         matlabbatch{1}.spm.tools.snpm.des.Corr1S.CovInt = temp_covint(find(FCRI>=13 & medicated~=1));
%                         matlabbatch{1}.spm.tools.snpm.des.Corr1S.P = matlabbatch{1}.spm.tools.snpm.des.Corr1S.P(find(FCRI>=13 & medicated~=1));
%                         matlabbatch{1}.spm.tools.snpm.des.Corr1S.xblock = length(find(FCRI>=13 & medicated~=1));
%                     else
%                         matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov.c = age(find(medicated~=1));
%                         matlabbatch{1}.spm.tools.snpm.des.Corr1S.CovInt = temp_covint(find(medicated~=1));
%                         matlabbatch{1}.spm.tools.snpm.des.Corr1S.xblock = length(FCRI(find(medicated~=1)));
%                     end
                    if strcmp(group_type{temp_group},'_FCRI_low')
                        matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov.c = age(find(FCRI<13));
                        matlabbatch{1}.spm.tools.snpm.des.Corr1S.CovInt = temp_covint(find(FCRI<13));
                        matlabbatch{1}.spm.tools.snpm.des.Corr1S.P = matlabbatch{1}.spm.tools.snpm.des.Corr1S.P(find(FCRI<13));
                        matlabbatch{1}.spm.tools.snpm.des.Corr1S.xblock = length(find(FCRI<13));
                    elseif strcmp(group_type{temp_group},'_FCRI_high')
                        matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov.c = age(find(FCRI>=13));
                        matlabbatch{1}.spm.tools.snpm.des.Corr1S.CovInt = temp_covint(find(FCRI>=13));
                        matlabbatch{1}.spm.tools.snpm.des.Corr1S.P = matlabbatch{1}.spm.tools.snpm.des.Corr1S.P(find(FCRI>=13));
                        matlabbatch{1}.spm.tools.snpm.des.Corr1S.xblock = length(find(FCRI>=13));
                    else
                        matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov.c = age;
                        matlabbatch{1}.spm.tools.snpm.des.Corr1S.CovInt = temp_covint;
                        matlabbatch{1}.spm.tools.snpm.des.Corr1S.xblock = length(FCRI);
                    end
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.cov.cname = 'Age';                 
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.nPerm = 10000;
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.vFWHM = [0 0 0];
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.bVolm = 1;
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.ST.ST_later = -1;
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.masking.tm.tm_none = 1;
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.masking.im = 1;
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.masking.em = {'H:\GMwithoutCER_61x73x61.nii,1'};
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.globalc.g_omit = 1;
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.globalm.gmsca.gmsca_no = 1;
                    matlabbatch{1}.spm.tools.snpm.des.Corr1S.globalm.glonorm = 1;
                    %% SnPM Compute
                    matlabbatch{2}.spm.tools.snpm.cp.snpmcfg = {strcat(output_Dir,'\SnPMcfg.mat')};
                    %% SnPM Inference
                    matlabbatch{3}.spm.tools.snpm.inference.SnPMmat = {strcat(output_Dir,'\SnPM.mat')};
                    matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = 0.01;
                    matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.PthC = 0.05;
                    matlabbatch{3}.spm.tools.snpm.inference.Tsign = 1;
                    matlabbatch{3}.spm.tools.snpm.inference.WriteFiltImg.name = strcat('FCS_significant_pos_corr_with_',behavioral{temp_behavioral});
                    matlabbatch{3}.spm.tools.snpm.inference.Report = 'MIPtable';
                    
                    matlabbatch{4}.spm.tools.snpm.inference.SnPMmat = {strcat(output_Dir,'\SnPM.mat')};
                    matlabbatch{4}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = 0.01;
                    matlabbatch{4}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.PthC = 0.05;
                    matlabbatch{4}.spm.tools.snpm.inference.Tsign = -1;
                    matlabbatch{4}.spm.tools.snpm.inference.WriteFiltImg.name = strcat('FCS_significant_neg_corr_with_',behavioral{temp_behavioral});
                    matlabbatch{4}.spm.tools.snpm.inference.Report = 'MIPtable';
                    
                    spm('defaults', 'FMRI'); % configure spm
                    spm_jobman('run', matlabbatch);
                end
            end
        end
    end
end