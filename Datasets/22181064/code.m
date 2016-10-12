%% Rodriguez~Cordero, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
rodriguez_cordero_2012.pmid = 22181064;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(rodriguez_cordero_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('textread','./raw_data/Rodr-guez-Porrata_et_al-2012-Journal_of_Applied_Microbiology.txt', '%s', 'delimiter', '\n');

hit_strains = {};

% Collect the strains
for i = 1:length(data)
    C = strsplit(data{i});
    hit_strains = [hit_strains; C{1}];
end

% Clean hit_strains
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];

% Create data matrix
hit_data = zeros(size(hit_strains)) - 1;


% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [130];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
rodriguez_cordero_2012.orfs = hit_strains;
rodriguez_cordero_2012.ph = hit_data_names;
rodriguez_cordero_2012.data = hit_data;
rodriguez_cordero_2012.dataset_ids = hit_data_ids;

%% Save

save('./rodriguez_cordero_2012.mat','rodriguez_cordero_2012');

%% Print out

fid = fopen('./rodriguez_cordero_2012.txt','w');
write_matrix_file(fid, rodriguez_cordero_2012.orfs, rodriguez_cordero_2012.ph, rodriguez_cordero_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(rodriguez_cordero_2012)
end

end

