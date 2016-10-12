%% Batova~Schuller, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
batova_schuller_2010.pmid = 20202201;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(batova_schuller_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/1471-2164-11-153-S1.xlsx','Sheet1');

% Isolate the hit strains
hit_strains = data(1:3:end);
hit_data = [data(2:3:end) data(3:3:end)];

% Clean hit_strains
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% Make a data matrix
hit_data(strcmp(hit_data,'+')) = {0};
hit_data(strcmp(hit_data,'sl')) = {-1};
hit_data(strcmp(hit_data,'?')) = {-2};

hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [153; 437];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
batova_schuller_2010.orfs = hit_strains;
batova_schuller_2010.ph = hit_data_names;
batova_schuller_2010.data = hit_data;
batova_schuller_2010.dataset_ids = hit_data_ids;

%% Save

save('./batova_schuller_2010.mat','batova_schuller_2010');

%% Print out

fid = fopen('./batova_schuller_2010.txt','w');
write_matrix_file(fid, batova_schuller_2010.orfs, batova_schuller_2010.ph, batova_schuller_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(batova_schuller_2010)
end

end

