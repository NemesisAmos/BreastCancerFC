tic
clear all;clc;
%% hyperparameters
SamplePeriod=0.7;
FilterBands=[0.01,0.1];
run_Parallel= true;
% run_Parallel= false;
FC_thres=0.2;
Lower_dist_thres=20;
Upper_dist_thres=180;
globalsig='nGSR';
WorkDir = 'H:\all_subjects_res\';

MaskName = strcat('E:\Matlab toolboxes\SeeCAT\templates\GMwithoutCER_61x73x61.nii');
thres_type={'zPos_Bin_'};
group_type={''};
if strcmp(globalsig,'nGSR')
    globalind='';
else
    globalind='global';
end
StartFolder = strcat('FunImgARWSD',globalind,'CFB\');
for temp_group=1:length(group_type)
    for temp_thres=1:length(thres_type)
        SeedsPath=strcat('H:\all_subjects_res\statistical_res\',globalsig,'\FCS_',thres_type{temp_thres},'group_diff',group_type{temp_group},'\',globalsig,'_FCS_significant_group_diff_SnPM.nii');
        ClusterConnectivityCriterion=26;
        [data,~,~,Header]=rest_to4d(SeedsPath);
        [cluster_belongings,~,cluster_sub,~,~,~]=amos_ClusterReport(data,Header,ClusterConnectivityCriterion,6,1);
        nSeed = length(cluster_belongings);
        %% start
        if run_Parallel
            %     CoreNum = 8; % 设置CPU核心数量
            %     parpool('local', CoreNum);
            if exist('gcp.m','file')
                try
                    gcp;
                end
            elseif  parpool('size') == 0
                try
                    parpool;
                end
            end
        end
        
        if ispc
            filesep = '\';
        else
            filesep = '/';
        end
        
        for temp_lower_band=1:length(FilterBands)-1                     
            fprintf('Calculating Functional Connectivity Started.\n');
            
            hdr_mask = spm_vol(MaskName);
            [vol_mask,XYZ] = spm_read_vols(hdr_mask);
            
            vol_seed = cell(nSeed,1);
            for i = 1:nSeed
                vol_seed{i} = zeros(hdr_mask.dim);
                vol_seed{i}(cluster_sub{i}) = 1;
            end
            
            Sublist = dir([WorkDir,filesep,StartFolder]);
            if strcmpi(Sublist(3).name,'.DS_Store')
                Sublist(1:3) = [];
            else
                Sublist(1:2) = [];
            end
            mask_ind = reshape(vol_mask>0,1,[]);
            XYZ_mask = XYZ(:,mask_ind);
            
            startFreq =['00',num2str(100*FilterBands(temp_lower_band))];
            startFreq = startFreq(end-2:end);
            endFreq =['00',num2str(100*(FilterBands(temp_lower_band+1)))];
            endFreq = endFreq(end-2:end);
            for temp_dist=1:length(Lower_dist_thres)
                mkdir([WorkDir,filesep,'Results',filesep,'Freq',startFreq,'-',endFreq,'_FC_',StartFolder(1:end-1),'_',num2str(Lower_dist_thres(temp_dist)),'-',num2str(Upper_dist_thres(temp_dist))]);
            end
            
            if run_Parallel
                parfor i = 1:length(Sublist)
                    fprintf(['Calculating ',Sublist(i).name,' Functional Connectivity.\n']);
                    cd([WorkDir,filesep,StartFolder,filesep,Sublist(i).name]);
                    Filename = dir('*.nii');
                    Nii = nifti(Filename.name);
                    volCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
                    maskCourse = volCourse(:,mask_ind);
                    seedCourse = zeros(Nii.dat.dim(1,4),nSeed);
                    for j = 1:nSeed
                        seedCourse(:,j) = mean(volCourse(:,reshape(vol_seed{j}>0,1,[])),2);
                    end
                    
                    [maskCourse_filtered] = y_IdealFilter(maskCourse, SamplePeriod, [FilterBands(temp_lower_band),FilterBands(temp_lower_band+1)]);
                    maskCourse_filtered(find(isnan(maskCourse_filtered)))=0;
                    maskCourse_filtered = detrend(maskCourse_filtered);%add detrend on masked data
                    maskCourse_filtered = maskCourse_filtered - repmat(mean(maskCourse_filtered),[Nii.dat.dim(1,4),1]);
                    maskCourse_filtered = maskCourse_filtered./repmat(std(maskCourse_filtered,0,1),[Nii.dat.dim(1,4),1]);
                    
                    [seedCourse_filtered] = y_IdealFilter(seedCourse, SamplePeriod, [FilterBands(temp_lower_band),FilterBands(temp_lower_band+1)]);
                    seedCourse_filtered(find(isnan(seedCourse_filtered)))=0;
                    seedCourse_filtered = detrend(seedCourse_filtered);%add detrend on masked data
                    seedCourse_filtered = seedCourse_filtered - repmat(mean(seedCourse_filtered),[Nii.dat.dim(1,4),1]);
                    seedCourse_filtered = seedCourse_filtered./repmat(std(seedCourse_filtered,0,1),[Nii.dat.dim(1,4),1]);
                    
                    r = seedCourse_filtered' * maskCourse_filtered ./(Nii.dat.dim(1,4)-1);
                    for temp_dist=1:length(Lower_dist_thres)
                        for j = 1:nSeed
                            XYZ_seed = XYZ(:,reshape(vol_seed{j}>0,1,[]));
                            D = pdist2(XYZ_seed',XYZ_mask');
                            D = mean(D,1);%obtain mean distance for seed
                            r_tmp=r(j,:);
                            r_tmp(find(r_tmp<FC_thres))=0;
                            r_tmp(D<=Lower_dist_thres(temp_dist)|D>Upper_dist_thres(temp_dist)) = 0;
                            r(j,:)=r_tmp;
                        end
                        z = FisherTrans(r);
                        cd([WorkDir,filesep,'Results',filesep,'Freq',startFreq,'-',endFreq,'_FC_',StartFolder(1:end-1),'_',num2str(Lower_dist_thres(temp_dist)),'-',num2str(Upper_dist_thres(temp_dist))]);
                        for j = 1:nSeed
                            hdr_fc = hdr_mask;
                            hdr_fc.fname = ['FC_Seed','_',cluster_belongings{j},'_',Sublist(i).name,'.nii'];
                            hdr_fc.dt(1) = 16;
                            vol_fc = zeros(hdr_fc.dim);
                            vol_fc(mask_ind) = r(j,:);
                            spm_write_vol(hdr_fc,vol_fc);
                            hdr_zfc = hdr_fc;
                            hdr_zfc.fname = ['z',hdr_zfc.fname];
                            vol_zfc = zeros(hdr_zfc.dim);
                            vol_zfc(mask_ind) = z(j,:);
                            spm_write_vol(hdr_zfc,vol_zfc);
                        end
                    end
                end
            else
                for i = 1:length(Sublist)
                    fprintf(['Calculating ',Sublist(i).name,' Functional Connectivity.\n']);
                    cd([WorkDir,filesep,StartFolder,filesep,Sublist(i).name]);
                    Filename = dir('*.nii');
                    Nii = nifti(Filename.name);
                    volCourse = reshape(double(Nii.dat),[Nii.dat.dim(1,1) * Nii.dat.dim(1,2) * Nii.dat.dim(1,3),Nii.dat.dim(1,4)])';
                    maskCourse = volCourse(:,mask_ind);
                    seedCourse = zeros(Nii.dat.dim(1,4),nSeed);
                    for j = 1:nSeed
                        seedCourse(:,j) = mean(volCourse(:,reshape(vol_seed{j}>0,1,[])),2);
                    end
                    
                    [maskCourse_filtered] = y_IdealFilter(maskCourse, SamplePeriod, [FilterBands(temp_lower_band),FilterBands(temp_lower_band+1)]);
                    maskCourse_filtered(find(isnan(maskCourse_filtered)))=0;
                    maskCourse_filtered = detrend(maskCourse_filtered);%add detrend on masked data
                    maskCourse_filtered = maskCourse_filtered - repmat(mean(maskCourse_filtered),[Nii.dat.dim(1,4),1]);
                    maskCourse_filtered = maskCourse_filtered./repmat(std(maskCourse_filtered,0,1),[Nii.dat.dim(1,4),1]);
                    
                    [seedCourse_filtered] = y_IdealFilter(seedCourse, SamplePeriod, [FilterBands(temp_lower_band),FilterBands(temp_lower_band+1)]);
                    seedCourse_filtered(find(isnan(seedCourse_filtered)))=0;
                    seedCourse_filtered = detrend(seedCourse_filtered);%add detrend on masked data
                    seedCourse_filtered = seedCourse_filtered - repmat(mean(seedCourse_filtered),[Nii.dat.dim(1,4),1]);
                    seedCourse_filtered = seedCourse_filtered./repmat(std(seedCourse_filtered,0,1),[Nii.dat.dim(1,4),1]);
                    
                    r = seedCourse_filtered' * maskCourse_filtered ./(Nii.dat.dim(1,4)-1);
                    for temp_dist=1:length(Lower_dist_thres)
                        for j = 1:nSeed
                            XYZ_seed = XYZ(:,reshape(vol_seed{j}>0,1,[]));
                            D = pdist2(XYZ_seed',XYZ_mask');
                            D = mean(D,1);%obtain mean distance for seed
                            r_tmp=r(j,:);
                            r_tmp(find(r_tmp<FC_thres))=0;
                            r_tmp(D<=Lower_dist_thres(temp_dist)|D>Upper_dist_thres(temp_dist)) = 0;
                            r(j,:)=r_tmp;
                        end
                        z = FisherTrans(r);
                        cd([WorkDir,filesep,'Results',filesep,'Freq',startFreq,'-',endFreq,'_FC_',StartFolder(1:end-1),'_',num2str(Lower_dist_thres(temp_dist)),'-',num2str(Upper_dist_thres(temp_dist))]);
                        for j = 1:nSeed
                            hdr_fc = hdr_mask;
                            hdr_fc.fname = ['FC_Seed','_',cluster_belongings{j},'_',Sublist(i).name,'.nii'];
                            hdr_fc.dt(1) = 16;
                            vol_fc = zeros(hdr_fc.dim);
                            vol_fc(mask_ind) = r(j,:);
                            spm_write_vol(hdr_fc,vol_fc);
                            hdr_zfc = hdr_fc;
                            hdr_zfc.fname = ['z',hdr_zfc.fname];
                            vol_zfc = zeros(hdr_zfc.dim);
                            vol_zfc(mask_ind) = z(j,:);
                            spm_write_vol(hdr_zfc,vol_zfc);
                        end
                    end
                end
            end
            fprintf('Calculating Functional Connectivity Finished.\n');
            toc
        end
        
        for temp_lower_band=1:length(FilterBands)-1
            startFreq =['00',num2str(100*FilterBands(temp_lower_band))];
            startFreq = startFreq(end-2:end);
            endFreq =['00',num2str(100*(FilterBands(temp_lower_band+1)))];
            endFreq = endFreq(end-2:end);
            for temp_dist=1:length(Lower_dist_thres)
                temp_path=[WorkDir,filesep,'Results',filesep,'Freq',startFreq,'-',endFreq,'_FC_',StartFolder(1:end-1),'_',num2str(Lower_dist_thres(temp_dist)),'-',num2str(Upper_dist_thres(temp_dist))];
                cd(temp_path);
                for j = 1:nSeed
                    if size(spm_select('fplist',temp_path,['zFC_Seed_',cluster_belongings{j},'.*nii']),1)>0
                        movefile(['zFC_Seed_',cluster_belongings{j},'*.nii'],[WorkDir,filesep,'Results',filesep,'Freq',startFreq,'-',endFreq,'_FC_',StartFolder(1:end-1),'_',num2str(Lower_dist_thres(temp_dist)),'-',num2str(Upper_dist_thres(temp_dist)),filesep,thres_type{temp_thres},'_zFC_Seed_',cluster_belongings{j}]);
                    end
                    if size(spm_select('fplist',temp_path,['FC_Seed_',cluster_belongings{j},'.*nii']),1)>0
                        movefile(['FC_Seed_',cluster_belongings{j},'*.nii'],[WorkDir,filesep,'Results',filesep,'Freq',startFreq,'-',endFreq,'_FC_',StartFolder(1:end-1),'_',num2str(Lower_dist_thres(temp_dist)),'-',num2str(Upper_dist_thres(temp_dist)),filesep,thres_type{temp_thres},'_FC_Seed_',cluster_belongings{j}]);
                    end
                end
                cd(WorkDir);
            end
        end
    end
end