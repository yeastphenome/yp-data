%% Dilda~Hogg, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
dilda_hogg_2008.pmid = 18160327;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(dilda_hogg_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, sensitive] = read_data('xlsread','./raw_data/Tables12.xlsx', 'Sensitivity');

% Get the list of ORFs and translate
sensitive_strains = sensitive(:,1);
sensitive_strains = translate(sensitive_strains);

% If possible, fix the problem (typos, omissions etc.)
sensitive_strains(ismember(sensitive_strains, {'CRS5'})) = {'YOR031W'};

% Get the data
sensitive_data = sensitive(:,2);
sensitive_data = cell2mat(sensitive_data);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(sensitive_strains));
sensitive_strains(inds) = [];
sensitive_data(inds) = [];

%% Now retrieve resistant tables
[FILENAMES{end+1}, resistant] = read_data('xlsread','./raw_data/Tables12.xlsx', 'Resistance');

% Get the list of ORFs and translate
resistant_strains = resistant(:,1);
resistant_strains = translate(resistant_strains);

% Get the data
resistant_data = resistant(:,2);
ind = find(~cellfun(@isnumeric, resistant_data));
resistant_data{ind} = [4000];
resistant_data = cell2mat(resistant_data);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(resistant_strains));
resistant_strains(inds) = [];
resistant_data(inds) = [];

%% Combine data

hit_strains = [resistant_strains; sensitive_strains];
hit_data = [resistant_data; sensitive_data];

% Transform data
hit_data = log(hit_data / 1100);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [1323];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
dilda_hogg_2008.orfs = hit_strains;
dilda_hogg_2008.ph = hit_data_names;
dilda_hogg_2008.data = hit_data;
dilda_hogg_2008.dataset_ids = hit_data_ids;

%% Save

save('./dilda_hogg_2008.mat','dilda_hogg_2008');

%% Print out

fid = fopen('./dilda_hogg_2008.txt','w');
write_matrix_file(fid, dilda_hogg_2008.orfs, dilda_hogg_2008.ph, dilda_hogg_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(dilda_hogg_2008)
end

end

