%% Botet~Santos, 2007
function FILENAMES = code()
% NOTES = data unnormalized; potentially batch, row/col normalization needed?

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
botet_santos_2007.pmid = 17873082;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(botet_santos_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

treatments = {'SMM, 77 h'; 'SMM, 120 h'; 'Sulfanilamide, 0.1 mg/ml, 77 h';'Sulfanilamide, 0.1 mg/ml, 120 h'};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/1ScreenSULFA&MS&MS+PABA.xlsx', 'DATA');

% Get indices of the data columns
ind_data = 15:18;

hit_strains = data.raw(:,10);
hit_data = data.raw(:,ind_data);

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YER050'})) = {'YER050C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

% Make sure all the data are numbers
inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {NaN};
hit_data = cell2mat(hit_data);

% Normalize by UNTREATED
hit_data(:,3) = hit_data(:,3)./hit_data(:,1);
hit_data(:,4) = hit_data(:,4)./hit_data(:,2);
hit_data(:,1:2) = [];
treatments(1:2) = [];

% Average data for identical ORFs that appear multiple times
[hit_strains,hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [137; 244];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
botet_santos_2007.orfs = hit_strains;
botet_santos_2007.ph = hit_data_names;
botet_santos_2007.data = hit_data;
botet_santos_2007.dataset_ids = hit_data_ids;

%% Save

save('./botet_santos_2007.mat','botet_santos_2007');

%% Print out

fid = fopen('./botet_santos_2007.txt','w');
write_matrix_file(fid, botet_santos_2007.orfs, botet_santos_2007.ph, botet_santos_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(botet_santos_2007)
end

end
