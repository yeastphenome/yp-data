%% Hoose~Polymenis, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
hoose_polymenis_2012.pmid = 22438835;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hoose_polymenis_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Dataset_S1.xlsx', 'Data and Correlation');

% Get the orfs and clean them up
strains = data(3:end,2);
strains = clean_orf(strains);

hit_data = data(3:end, 4);

% Find anything that doesn't look like an ORF
strains(ismember(strains, {'YELOO1C'})) = {'YEL001C'};
inds = find(~is_orf(strains));
strains(inds) = [];
hit_data(inds) = [];

%% Get data from hits
% Make a data matrix
hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[strains, hit_data] = grpstats(hit_data, strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [515];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

hoose_polymenis_2012.orfs = strains;
hoose_polymenis_2012.ph = hit_data_names;
hoose_polymenis_2012.data = hit_data;
hoose_polymenis_2012.dataset_ids = hit_data_ids;

%% Save

save('./hoose_polymenis_2012.mat','hoose_polymenis_2012');

%% Print out

fid = fopen('./hoose_polymenis_2012.txt','w');
write_matrix_file(fid, hoose_polymenis_2012.orfs, hoose_polymenis_2012.ph, hoose_polymenis_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(hoose_polymenis_2012)
end

end
