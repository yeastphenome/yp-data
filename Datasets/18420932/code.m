%% Hillenmeyer~Giaever, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hillenmeyer_giaever_2008.pmid = 18420932;

%% Data

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
[FILENAMES{end+1}, datasets] = read_data('textread','./extras/datasets.txt','%d %s','delimiter','\t');
dataset_id = datasets{1};
column_name = datasets{2};
[~,ind1,ind2] = intersect(column_name, data.labels_col);
data.dataset_id = nan(size(data.data,2),1);
data.dataset_id(ind2) = dataset_id(ind1);

[hit_datasetids, hit_data] = grpstats(hit_data', data.dataset_id, {'gname','mean'});
hit_data = hit_data';
hit_datasetids = cellfun(@str2num, hit_datasetids);

inds = find(hit_datasetids==0);
hit_datasetids(inds) = [];
hit_data(:,inds) = [];

% Load the datasetID-to-datasetName mapping
[FILENAMES{end+1}, datasets2] = read_data('textread', './extras/YeastPhenome_18420932_datasets_list.txt','%d %s','delimiter','\t');
dataset_id = datasets2{1};
dataset_name = datasets2{2};
[~,ind1,ind2] = intersect(dataset_id, hit_datasetids);
hit_datasetnames = cell(size(hit_datasetids));
hit_datasetnames(ind2) = dataset_name(ind1);

%% Prepare final dataset

hillenmeyer_giaever_2008.orfs = hit_strains;
hillenmeyer_giaever_2008.ph = hit_datasetnames;
hillenmeyer_giaever_2008.data = hit_data;
hillenmeyer_giaever_2008.dataset_ids = hit_datasetids;

%% Save

save('./hillenmeyer_giaever_2008.mat','hillenmeyer_giaever_2008');

%% Print out

fid = fopen('./hillenmeyer_giaever_2008.txt','w');
write_matrix_file(fid, hillenmeyer_giaever_2008.orfs, hillenmeyer_giaever_2008.ph, hillenmeyer_giaever_2008.data);
fclose(fid);

end

