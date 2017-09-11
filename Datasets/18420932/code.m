%% Hillenmeyer~Giaever, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hillenmeyer_giaever_2008.pmid = 18420932;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hillenmeyer_giaever_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hom data

% Load the data
[FILENAMES{end+1}, data] = read_data('read_matrix_file','./raw_data/hom.ratio_result_nm.pub');

% Get the list of ORFs
hit_strains = strtok(data.labels_row(:,1),':');

% Get the data itself
hit_data = -data.data;  % taking the opposite because the original values are log2(control/treatment)
    
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% Load the column-to-dataset mapping and average the columns that are
% replicates and map to the same dataset
[FILENAMES{end+1}, datasets_hom] = read_data('textread','./extras/datasets.txt','%d %s','delimiter','\t');
dataset_id = datasets_hom{1};
column_name = datasets_hom{2};
[~,ind1,ind2] = intersect(column_name, data.labels_col);
data.dataset_id = nan(size(data.data,2),1);
data.dataset_id(ind2) = dataset_id(ind1);

[hit_datasetids, hit_data] = grpstats(hit_data', data.dataset_id, {'gname','mean'});
hit_data = hit_data';
hit_datasetids = cellfun(@str2num, hit_datasetids);

inds = find(hit_datasetids==0);
hit_datasetids(inds) = [];
hit_data(:,inds) = [];

%% Het data

% Load the data
[FILENAMES{end+1}, data_het] = read_data('read_matrix_file','./raw_data/het.ratio_result_nm.pub');

% Get the list of ORFs
hit_strains_het = strtok(data_het.labels_row(:,1),':');

% Get the data itself
hit_data_het = -data_het.data;  % taking the opposite because the original values are log2(control/treatment)
    
% Eliminate all white spaces & capitalize
hit_strains_het = clean_orf(hit_strains_het);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_het));
disp(hit_strains_het(inds));  

% If the same strain is present more than once, average its values
[hit_strains_het, hit_data_het] = grpstats(hit_data_het, hit_strains_het, {'gname','mean'});

% Load the column-to-dataset mapping and average the columns that are
% replicates and map to the same dataset
[FILENAMES{end+1}, het_datasets] = read_data('readtable','./extras/datasets_het.txt','delimiter','\t');
dataset_id = het_datasets.Dataset;
column_name = het_datasets.Name;
conditionset_id = het_datasets.Conditionset;

% Delete some of the datasets that couldn't be parsed
inds = find(dataset_id == 0);
dataset_id(inds) = [];
column_name(inds) = [];
conditionset_id(inds) = [];

conditionset_id = cellfun(@str2num, conditionset_id);

[~,ind1,ind2] = intersect(column_name, data_het.labels_col);
data_het.dataset_id = nan(size(data_het.data,2),1);
data_het.dataset_id(ind2) = dataset_id(ind1);

[hit_datasetids_het, hit_data_het] = grpstats(hit_data_het', data_het.dataset_id, {'gname','mean'});
hit_data_het = hit_data_het';
hit_datasetids_het = cellfun(@str2num, hit_datasetids_het);

inds = find(hit_datasetids_het==0);
hit_datasetids_het(inds) = [];
hit_data_het(:,inds) = [];

%% Join the 2 datasets

hit_strains_all = unique([hit_strains; hit_strains_het]);
hit_data_all = nan(length(hit_strains_all), length(datasets.id));
hit_dataset_ids_all = [hit_datasetids; hit_datasetids_het];

[~,ind1,ind2] = intersect(hit_strains, hit_strains_all);
hit_data_all(ind2,1:length(hit_datasetids)) = hit_data(ind1,:);

[~,ind1,ind2] = intersect(hit_strains_het, hit_strains_all);
hit_data_all(ind2,length(hit_datasetids)+1:end) = hit_data_het(ind1,:);

%% Prepare final dataset

[~,ind1,ind2] = intersect(hit_dataset_ids_all, datasets.id);
hit_datasetnames(ind1) = datasets.standard_name(ind2);
hit_datasetnames = hit_datasetnames';

hillenmeyer_giaever_2008.orfs = hit_strains_all;
hillenmeyer_giaever_2008.ph = hit_datasetnames;
hillenmeyer_giaever_2008.data = hit_data_all;
hillenmeyer_giaever_2008.dataset_ids = hit_dataset_ids_all;

%% Save

save('./hillenmeyer_giaever_2008.mat','hillenmeyer_giaever_2008');

%% Print out

fid = fopen('./hillenmeyer_giaever_2008.txt','w');
write_matrix_file(fid, hillenmeyer_giaever_2008.orfs, hillenmeyer_giaever_2008.ph, hillenmeyer_giaever_2008.data);
fclose(fid);

end

