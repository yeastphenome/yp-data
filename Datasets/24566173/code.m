%% Bleackley~Anderson, 2014
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
bleackley_anderson_2014.pmid = 24566173;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(bleackley_anderson_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/full list NaD1 resistance.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,7);

% Get the data itself
hit_data = data(2:end,6);
hit_data = cell2mat(hit_data);

% Transform the data to fit
hit_data = hit_data * -1;
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
hit_strains(ismember(hit_strains, {'YLR287-A'})) = {'YLR287C-A'};
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [757];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
bleackley_anderson_2014.orfs = hit_strains;
bleackley_anderson_2014.ph = hit_data_names;
bleackley_anderson_2014.data = hit_data;
bleackley_anderson_2014.dataset_ids = hit_data_ids;

%% Save

save('./bleackley_anderson_2014.mat','bleackley_anderson_2014');

%% Print out

fid = fopen('./bleackley_anderson_2014.txt','w');
write_matrix_file(fid, bleackley_anderson_2014.orfs, bleackley_anderson_2014.ph, bleackley_anderson_2014.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(bleackley_anderson_2014)
end

end

