%generate subject level pair-wise FC between all signficant FCS and FC correlate with behavioral clusters
clear all;clc;
%generate statistical correlation pictures figure
data_root=strcat('D:\BaiduNetdiskDownload\breast_cancer_T1_fMRI_40\');
globalsignal='GSR';
thres_type={'Pos'};
behavioral_features={'FCRI','FCRtotal'};
[behavioral_data,all_txt]=xlsread(strcat(data_root,'breast_cancer_subj_info.xlsx'));
subject_list=all_txt(2:end,1);
subject_num=length(subject_list);
Yeomaskroot=strcat('E:\Matlab Projects\BCT\2016_01_16_BCT\yeo\Yeo_JNeurophysiol11_MNI152\Yeo2011_7Networks_MNI152_FreeSurferConformed3mm_LiberalMask.nii');%only calculate FCS on GM mask
Yeomaskdata=rest_to4d(Yeomaskroot);
Yeolabelnames={'Visual','Smatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontal Parietal','Default Mode'};
for temp_thres=1:length(thres_type)
    for temp_behavioral=1:length(behavioral_features)
        for temp_subject=1:subject_num
            %load orig subject data for further processing
            if strcmp(globalsignal,'GSR')
                globalind='global';
            else
                globalind='';
            end
            temp_subject_path=strcat(data_root,'FunImgARWSD',globalind,'CFB\',subject_list{temp_subject},'\',subject_list{temp_subject},'_4DVolume.nii');
            temp_subject_data=rest_to4d(temp_subject_path);
            temp_subject_data=reshape(temp_subject_data,size(temp_subject_data,1)*size(temp_subject_data,2)*size(temp_subject_data,3),[]);
            temp_subject_allclusters=[];
            temp_subject_allcluster_belongings={};
            temp_subject_allcluster_voxel_num=[];
            temp_subject_allcluster_map_to_Yeo={};
            temp_subject_allcluster_map_to_Yeo_index=[];
            all_x=[];
            all_y=[];
            all_z=[];
            all_cluster_num=0;
            
            temp_FCS_path=strcat(data_root,'statistical_res\',globalsignal,'\FCS_Pos_corr_with_',behavioral_features{temp_behavioral},'\FCS_significant_corr_with_',behavioral_features{temp_behavioral},'.nii');
            if logical(exist(temp_FCS_path,'file'))
                ClusterConnectivityCriterion=26;
                [temp_FCS_data,~,~,Header]=rest_to4d(temp_FCS_path);
                [FCS_cluster_belongings,~,FCS_cluster_sub,~,FCS_cluster_peakcoords,~]=amos_ClusterReport(temp_FCS_data,Header,ClusterConnectivityCriterion,6,1);
                for temp_FCS_cluster_ind=1:length(FCS_cluster_belongings)
                    all_cluster_num=all_cluster_num+1;
                    max_overlap_netowrk='None';
                    max_overlap_count=0;
                    max_overlap_index=0;
                    for tempnetwork=1:max(unique(Yeomaskdata))
                        tempnetworkind=find(Yeomaskdata==tempnetwork);
                        temp_overlap_count=length(intersect(FCS_cluster_sub{temp_FCS_cluster_ind},tempnetworkind));
                        if temp_overlap_count>max_overlap_count
                            max_overlap_count=temp_overlap_count;
                            max_overlap_netowrk=Yeolabelnames{tempnetwork};
                            max_overlap_index=tempnetwork;
                        end
                    end
                    temp_subject_allcluster_map_to_Yeo=[temp_subject_allcluster_map_to_Yeo;max_overlap_netowrk];
                    temp_subject_allcluster_map_to_Yeo_index=[temp_subject_allcluster_map_to_Yeo_index;max_overlap_index];
                    temp_subject_allcluster_belongings=[temp_subject_allcluster_belongings;FCS_cluster_belongings{temp_FCS_cluster_ind}];
                    temp_subject_allcluster_voxel_num=[temp_subject_allcluster_voxel_num;length(FCS_cluster_sub{temp_FCS_cluster_ind})];
                    temp_cluster_peakxyz=FCS_cluster_peakcoords{temp_FCS_cluster_ind};
                    all_x=[all_x;temp_cluster_peakxyz(1)];
                    all_y=[all_y;temp_cluster_peakxyz(2)];
                    all_z=[all_z;temp_cluster_peakxyz(3)];
                    %extract FCS cluster mean time series for all subject data
                    temp_subject_temp_cluster_ts=mean(temp_subject_data(FCS_cluster_sub{temp_FCS_cluster_ind},:));
                    temp_subject_allclusters=[temp_subject_allclusters;temp_subject_temp_cluster_ts];
                    temp_FC_path=strcat(data_root,'statistical_res\',globalsignal,'\',thres_type{temp_thres},'_',behavioral_features{temp_behavioral},'_FC_corr_with_',...
                        behavioral_features{temp_behavioral},'\',FCS_cluster_belongings{temp_FCS_cluster_ind},'\FC_significant_corr_with_',behavioral_features{temp_behavioral},'.nii');
                    if logical(exist(temp_FC_path,'file'))
                        temp_FC_data=rest_to4d(temp_FC_path);
                        [FC_cluster_belongings,~,FC_cluster_sub,~,FC_cluster_peakcoords,~]=amos_ClusterReport(temp_FC_data,Header,ClusterConnectivityCriterion,6,1);
                        for temp_FC_cluster_ind=1:length(FC_cluster_belongings)
                            all_cluster_num=all_cluster_num+1;
                            max_overlap_netowrk='None';
                            max_overlap_count=0;
                            max_overlap_index=0;
                            for tempnetwork=1:max(unique(Yeomaskdata))
                                tempnetworkind=find(Yeomaskdata==tempnetwork);
                                temp_overlap_count=length(intersect(FC_cluster_sub{temp_FC_cluster_ind},tempnetworkind));
                                if temp_overlap_count>max_overlap_count
                                    max_overlap_count=temp_overlap_count;
                                    max_overlap_netowrk=Yeolabelnames{tempnetwork};
                                    max_overlap_index=tempnetwork;
                                end
                            end
                            temp_subject_allcluster_map_to_Yeo=[temp_subject_allcluster_map_to_Yeo;max_overlap_netowrk];
                            temp_subject_allcluster_map_to_Yeo_index=[temp_subject_allcluster_map_to_Yeo_index;max_overlap_index];
                            temp_subject_allcluster_belongings=[temp_subject_allcluster_belongings;FC_cluster_belongings{temp_FC_cluster_ind}];
                            temp_subject_allcluster_voxel_num=[temp_subject_allcluster_voxel_num;length(FC_cluster_sub{temp_FC_cluster_ind})];
                            temp_cluster_peakxyz=FC_cluster_peakcoords{temp_FC_cluster_ind};
                            all_x=[all_x;temp_cluster_peakxyz(1)];
                            all_y=[all_y;temp_cluster_peakxyz(2)];
                            all_z=[all_z;temp_cluster_peakxyz(3)];
                            %extract FC cluster mean time series for all subject data
                            temp_subject_temp_cluster_ts=mean(temp_subject_data(FC_cluster_sub{temp_FC_cluster_ind},:));
                            temp_subject_allclusters=[temp_subject_allclusters;temp_subject_temp_cluster_ts];
                        end
                    end
                end
            end
            %remove space, undefined and BA xx in belonging names, and obtain first AAL region name
            for temp_belonging_index=1:all_cluster_num
                pattern = '(\s|\undefined_|\BA.*?_)';
                temp_name = temp_subject_allcluster_belongings{temp_belonging_index};               
                temp_name = regexprep(temp_name, pattern, '');%remove space
                temp_subject_allcluster_belongings{temp_belonging_index}=temp_name;
            end
            %remove small clusters with same names
%             remove_indexs=[];
%             for temp_cluster_belonging_x=1:all_cluster_num
%                 for temp_cluster_belonging_y=1:all_cluster_num
%                     if strcmp(temp_subject_allcluster_belongings{temp_cluster_belonging_x},temp_subject_allcluster_belongings{temp_cluster_belonging_y})
%                         if temp_subject_allcluster_voxel_num(temp_cluster_belonging_x)<temp_subject_allcluster_voxel_num(temp_cluster_belonging_y)
%                             remove_indexs=[remove_indexs;temp_cluster_belonging_x];
%                         end
%                     end
%                 end
%             end
%             remove_indexs=unique(remove_indexs);
%             temp_subject_allcluster_map_to_Yeo(remove_indexs)=[];
%             temp_subject_allcluster_map_to_Yeo_index(remove_indexs)=[];
%             temp_subject_allcluster_belongings(remove_indexs)=[];
%             temp_subject_allcluster_voxel_num(remove_indexs)=[];
%             all_x(remove_indexs)=[];
%             all_y(remove_indexs)=[];
%             all_z(remove_indexs)=[];
%             temp_subject_allclusters(remove_indexs,:)=[];
%             all_cluster_num=all_cluster_num-length(remove_indexs);

            temp_subject_FCmatrix=corr(temp_subject_allclusters');
            if temp_subject==1
                all_subject_FCmatrix=zeros(subject_num,all_cluster_num,all_cluster_num);
                z_all_subject_FCmatrix=zeros(subject_num,all_cluster_num,all_cluster_num);
            end
            all_subject_FCmatrix(temp_subject,:,:)=temp_subject_FCmatrix;
            z_all_subject_FCmatrix(temp_subject,:,:)=FisherTrans(temp_subject_FCmatrix);
        end
        %aveage all subject pair-wise FC results
        mean_all_subject_FCmatrix=mean(all_subject_FCmatrix);
        mean_all_subject_FCmatrix(find(isinf(mean_all_subject_FCmatrix)))=1;
        mean_all_subject_FCmatrix=squeeze(mean_all_subject_FCmatrix);
        
        mean_z_all_subject_FCmatrix=mean(z_all_subject_FCmatrix);
        mean_z_all_subject_FCmatrix(find(isinf(mean_z_all_subject_FCmatrix)))=1;
        mean_z_all_subject_FCmatrix=squeeze(mean_z_all_subject_FCmatrix);
        %use louvain or k-means cluster to identify network level modularities
        allmodules=zeros(100,size(mean_all_subject_FCmatrix,2));
        for temp_node=1:100
            [allmodules(temp_node,:),Q(temp_node,:)]=modularity_louvain_und_sign(mean_all_subject_FCmatrix,'sta');
        end
        [maxvalue,maxindex]=max(Q);
        connection_module=allmodules(maxindex,:); 
        %re-organize ROI based on community index
        [sorted_connection_module,sorted_connection_indexs]=sort(connection_module);
        temp_subject_allcluster_map_to_Yeo=temp_subject_allcluster_map_to_Yeo(sorted_connection_indexs);
        temp_subject_allcluster_map_to_Yeo_index=temp_subject_allcluster_map_to_Yeo_index(sorted_connection_indexs);
        module_mapped_belonging_names=temp_subject_allcluster_belongings(sorted_connection_indexs);
        mean_all_subject_FCmatrix=mean_all_subject_FCmatrix(sorted_connection_indexs,:);
        mean_all_subject_FCmatrix=mean_all_subject_FCmatrix(:,sorted_connection_indexs);
        imagesc(mean_all_subject_FCmatrix);
        saveas(gcf,strcat(data_root,'statistical_res\',globalsignal,'\',behavioral_features{temp_behavioral},'_all_subject_mean_cluster_pair_wise_FC_matrix.jpg'));
        close(gcf);   
        save(strcat(data_root,'statistical_res\',globalsignal,'\',behavioral_features{temp_behavioral},'_all_subject_mean_cluster_pair_wise_FC_results.mat'),'connection_module',...
            'mean_z_all_subject_FCmatrix','z_all_subject_FCmatrix','all_subject_FCmatrix','mean_all_subject_FCmatrix','module_mapped_belonging_names','temp_subject_allcluster_map_to_Yeo');       
        %连接矩阵生成txt或者mat，名字community等信息生成csv 用于 brant circos plot,
        mean_all_subject_FCmatrix(find(mean_all_subject_FCmatrix<0.2))=0;
        dlmwrite(strcat(data_root,'statistical_res\',globalsignal,'\',behavioral_features{temp_behavioral},'_all_subject_mean_cluster_pair_wise_FC.txt'),mean_all_subject_FCmatrix);
        node_num = length(sorted_connection_module);
        node_module_name = cell(node_num,1);
        node_in_module_index = zeros(node_num,1);
        node_in_Yeo_module_index = zeros(node_num,1);  
        node_in_Yeo_module_stack_count = zeros(length(Yeolabelnames),1);
        for temp_node = 1:node_num
            if temp_node ==1
                node_in_module_index(temp_node)=1;
            else
                if (sorted_connection_module(temp_node)-sorted_connection_module(temp_node-1))==0
                    node_in_module_index(temp_node)=node_in_module_index(temp_node-1)+1;
                else
                    node_in_module_index(temp_node)=1;
                end
            end
            node_module_name{temp_node}=['community',num2str(sorted_connection_module(temp_node))];
            if temp_subject_allcluster_map_to_Yeo_index(temp_node)>0
                node_in_Yeo_module_stack_count(temp_subject_allcluster_map_to_Yeo_index(temp_node))=node_in_Yeo_module_stack_count(temp_subject_allcluster_map_to_Yeo_index(temp_node))+1;
                node_in_Yeo_module_index(temp_node)=node_in_Yeo_module_stack_count(temp_subject_allcluster_map_to_Yeo_index(temp_node));
            else
                node_in_Yeo_module_index(temp_node)=0;% not index for None map with Yeo atlas regions
            end
        end
        indexes=1:node_num; indexes=indexes';
        variablenames={'x','y','z','label_old','label','vox_num','index','module','index_module','index_node','Yeo_module','index_Yeo_module','index_node_Yeo_module'};  
        savetable=table(all_x,all_y,all_z,module_mapped_belonging_names,module_mapped_belonging_names,temp_subject_allcluster_voxel_num,indexes,node_module_name,sorted_connection_module',node_in_module_index,...
            temp_subject_allcluster_map_to_Yeo,temp_subject_allcluster_map_to_Yeo_index,node_in_Yeo_module_index,'VariableNames',variablenames);
        writetable(savetable,strcat(data_root,'statistical_res\',globalsignal,'\',behavioral_features{temp_behavioral},'_all_subject_mean_cluster_pair_wise_FC_table.csv'));
    end
end


