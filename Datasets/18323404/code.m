%% Loukin~Saimi, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
loukin_saimi_2008.pmid = 18323404;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(loukin_saimi_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/07-101410SuppTableS1.xlsx');

% Get the list of ORFs and the correponding data 
hit_strains = data(2:end,1);

% Get the data itself
hit_data = data(2:end,3:4);
hit_data = cell2mat(hit_data);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5234; 5233];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
loukin_saimi_2008.orfs = hit_strains;
loukin_saimi_2008.ph = hit_data_names;
loukin_saimi_2008.data = hit_data;
loukin_saimi_2008.dataset_ids = hit_data_ids;

%% Save

save('./loukin_saimi_2008.mat','loukin_saimi_2008');

%% Print out

fid = fopen('./loukin_saimi_2008.txt','w');
write_matrix_file(fid, loukin_saimi_2008.orfs, loukin_saimi_2008.ph, loukin_saimi_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(loukin_saimi_2008)
end

end

