%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
bonangelino_bonifacino_2002.pmid = 12134085;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(bonangelino_bonifacino_2002.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/hit_list.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,4);

% Get the data itself
hit_data = data(:,3); % if the dataset is discrete or binary
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(strcmp(hit_strains, 'YKL054')) = {'YKL054C'};
hit_strains(strcmp(hit_strains, 'YDR525W--')) = {'YDR525W-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

hit_data = strtrim(hit_data);

hit_data(strcmp('W', hit_data)) = {-1};
hit_data(strcmp('M', hit_data)) = {-2};
hit_data(strcmp('S', hit_data)) = {-3};
hit_data(strcmp('W-M', hit_data)) = {-1.5};

hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = 469;

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
bonangelino_bonifacino_2002.orfs = hit_strains;
bonangelino_bonifacino_2002.ph = hit_data_names;
bonangelino_bonifacino_2002.data = hit_data;
bonangelino_bonifacino_2002.dataset_ids = hit_data_ids;

%% Save

save('./bonangelino_bonifacino_2002.mat','bonangelino_bonifacino_2002');

%% Print out

fid = fopen('./bonangelino_bonifacino_2002.txt','w');
write_matrix_file(fid, bonangelino_bonifacino_2002.orfs, bonangelino_bonifacino_2002.ph, bonangelino_bonifacino_2002.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(bonangelino_bonifacino_2002)
end

end

