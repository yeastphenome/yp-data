%% Schlecht~Stonge, 2014
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
schlecht_stonge_2014.pmid = 24708151;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(schlecht_stonge_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load Hit Strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread', './raw_data/Table.S1.xlsx', 'Table S1');

% Get the list of ORFs
hit_strains = data(2:end, 3);

% Clean up ORFs
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds)); 

% Get data from hits and standardize it
hit_data = data(2:end, 8:2:14);
hit_data(strcmp('NA', hit_data)) = {NaN};
hit_data = cell2mat(hit_data);

% Average any repeated value
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [692; 453; 690; 691];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
schlecht_stonge_2014.orfs = hit_strains;
schlecht_stonge_2014.ph = hit_data_names;
schlecht_stonge_2014.data = hit_data;
schlecht_stonge_2014.dataset_ids = hit_data_ids;

%% Save

save('./schlecht_stonge_2014.mat','schlecht_stonge_2014');

%% Print out

fid = fopen('./schlecht_stonge_2014.txt','w');
write_matrix_file(fid, schlecht_stonge_2014.orfs, schlecht_stonge_2014.ph, schlecht_stonge_2014.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(schlecht_stonge_2014)
end

end