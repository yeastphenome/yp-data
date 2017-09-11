%% Lis~Bobek, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
lis_bobek_2013.pmid = 23208710;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(lis_bobek_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread', './raw_data/AAC.01439-12_zac999101582so1.xlsx', 'all data');

% Get the list of ORFs
hit_strains = data(2:end, 1);

% Clean up ORFs
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds)); 

% Get data from hits
hit_data = -cell2mat(data(2:end, 2:5)); 

% Average any repeated value
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% Split essentials and non-essentials
hit_strains_all = hit_strains;
hit_data_all = nan(length(hit_strains), size(hit_data,2)*2);

load essential_genes_151215.mat
inds_essential = find(ismember(hit_strains_all, essential_genes));
hit_data_all(inds_essential,5:8) = hit_data(inds_essential,:);

inds_nonessential = find(~ismember(hit_strains_all, essential_genes));
hit_data_all(inds_nonessential,1:4) = hit_data(inds_nonessential,:);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [516 4925 4926 4927 11870 11871 11872 11873]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
lis_bobek_2013.orfs = hit_strains_all;
lis_bobek_2013.ph = hit_data_names;
lis_bobek_2013.data = hit_data_all;
lis_bobek_2013.dataset_ids = hit_data_ids;

%% Save

save('./lis_bobek_2013.mat','lis_bobek_2013');

%% Print out

fid = fopen('./lis_bobek_2013.txt','w');
write_matrix_file(fid, lis_bobek_2013.orfs, lis_bobek_2013.ph, lis_bobek_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(lis_bobek_2013)
end


end
