%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
herrero_moreno_2008.pmid = 19047157;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(herrero_moreno_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data_s] = read_data('xlsread','./raw_data/s_table_1.xlsx', 'Table 1');
[FILENAMES{end+1}, data_r] = read_data('xlsread','./raw_data/stab_3.xlsx', 'Table 1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains_s = data_s(:,[2 5 8 11]);
hit_strains_s = reshape(hit_strains_s, [], 1);
hit_strains_s(cellfun(@isnumeric, hit_strains_s)) = [];

hit_strains_r = data_r(:,2);
hit_strains_r(cellfun(@isnumeric, hit_strains_r)) = [];

% Eliminate all white spaces & capitalize
hit_strains_s = clean_orf(hit_strains_s);
hit_strains_r = clean_orf(hit_strains_r);

% If in gene name form, transform into ORF name
hit_strains_s = translate(hit_strains_s);
hit_strains_r = translate(hit_strains_r);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_s));
disp(hit_strains_s(inds));  

inds = find(~is_orf(hit_strains_r));
disp(hit_strains_r(inds));  
hit_strains_r(inds) = [];

% Merge sensitive and resistant lists
hit_strains = [hit_strains_s; hit_strains_r];
hit_data = [-ones(size(hit_strains_s)); ones(size(hit_strains_r))];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [106];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
herrero_moreno_2008.orfs = hit_strains;
herrero_moreno_2008.ph = hit_data_names;
herrero_moreno_2008.data = hit_data;
herrero_moreno_2008.dataset_ids = hit_data_ids;

%% Save

save('./herrero_moreno_2008.mat','herrero_moreno_2008');

%% Print out

fid = fopen('./herrero_moreno_2008.txt','w');
write_matrix_file(fid, herrero_moreno_2008.orfs, herrero_moreno_2008.ph, herrero_moreno_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(herrero_moreno_2008)
end

end

