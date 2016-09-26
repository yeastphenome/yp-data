%% Pan~Boeke, 2010

function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
pan_boeke_2010.pmid = 20660648;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(pan_boeke_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the Hit Data

% Hypersensitive strains
[FILENAMES{end+1}, data] = read_data('textread','./raw_data/TableS1.txt', '%s','delimiter','\n');

% Get the list of ORFs
hit_strains = {};

for i = 1:length(data)
    C = regexp(data{i},' ','split');
    hit_strains = [hit_strains; C{1}];
end

% Get the data itself
hit_data = ones(size(hit_strains));
hit_data = -2 * hit_data;

% Sensitive strains
[FILENAMES{end+1}, data] = read_data('textread','./raw_data/TableS2.txt', '%s','delimiter','\n');

% Get the list of ORFs
for i = 1:length(data)
    C = regexp(data{i},' ','split');
    hit_strains = [hit_strains; C{1}];
end

sens_data = ones(size(data));
sens_data = -1 * sens_data;

% Combine hypersensitive and sensistive data
hit_data = vertcat(hit_data, sens_data);

% Resistant strains
[FILENAMES{end+1}, data] = read_data('textread','./raw_data/TableS3.txt', '%s','delimiter','\n');

% Get the list of ORFs
for i = 1:length(data)
    C = regexp(data{i},' ','split');
    hit_strains = [hit_strains; C{1}];
end

res_data = ones(size(data));

% Combine data
hit_data = vertcat(hit_data, res_data);

% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [109];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
pan_boeke_2010.orfs = hit_strains;
pan_boeke_2010.ph = hit_data_names;
pan_boeke_2010.data = hit_data;
pan_boeke_2010.dataset_ids = hit_data_ids;

%% Save

save('./pan_boeke_2010.mat','pan_boeke_2010');

%% Print out

fid = fopen('./pan_boeke_2010.txt','w');
write_matrix_file(fid, pan_boeke_2010.orfs, pan_boeke_2010.ph, pan_boeke_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(pan_boeke_2010)
end

end

