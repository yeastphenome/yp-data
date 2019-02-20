%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
holland_avery_2007.pmid = 18088421;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(holland_avery_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/13059_2007_1748_MOESM1_ESM.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(14:end,1);

% Get the data itself
hit_data = data(14:end,4);
inds = find(~cellfun(@isnumeric, hit_data));
hit_strains(inds) = [];
hit_data(inds,:) = [];
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

inds = find(~cellfun(@isstr, hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16258];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
holland_avery_2007.orfs = hit_strains;
holland_avery_2007.ph = hit_data_names;
holland_avery_2007.data = hit_data;
holland_avery_2007.dataset_ids = hit_data_ids;

%% Save

save('./holland_avery_2007.mat','holland_avery_2007');

%% Print out

fid = fopen('./holland_avery_2007.txt','w');
write_matrix_file(fid, holland_avery_2007.orfs, holland_avery_2007.ph, holland_avery_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(holland_avery_2007)
end

end

