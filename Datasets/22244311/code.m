%% Pir~Oliver, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
pir_oliver_2012.pmid = 22244311;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(pir_oliver_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/12918_2012_852_MOESM3_ESM.xlsx', 'data');

% Get the list of ORFs and the correponding data 
hit_strains = data(2:end,1);

% Find the columns of interest
columnData = data(1,:);
ind = find(~cellfun(@isempty, strfind(columnData, '.FCC')));

% Get the data itself
hit_data = data(2:end,ind);

% Get rid of non-numeric indices
hit_data(~cellfun(@isnumeric, hit_data)) = {NaN};
hit_data = cell2mat(hit_data);

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11859; 11860; 11861; 15991; 15992; 15995; 15996];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
pir_oliver_2012.orfs = hit_strains;
pir_oliver_2012.ph = hit_data_names;
pir_oliver_2012.data = hit_data;
pir_oliver_2012.dataset_ids = hit_data_ids;

%% Save

save('./pir_oliver_2012.mat','pir_oliver_2012');

%% Print out

fid = fopen('./pir_oliver_2012.txt','w');
write_matrix_file(fid, pir_oliver_2012.orfs, pir_oliver_2012.ph, pir_oliver_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(pir_oliver_2012)
end

end