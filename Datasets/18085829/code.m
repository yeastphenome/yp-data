%% Alvaro~Rothstein, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
alvaro_rothstein_2007.pmid = 18085829;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(alvaro_rothstein_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/pgen.0030228.st001.xlsx', 'SIA-Rad52 focus screen');

% Get the list of ORFs and the correponding data 
hit_strains = data(9:end,4);

% Get the data itself
hit_data = data(9:end,3);
hit_data(find(~cellfun(@isnumeric, hit_data))) = {NaN};
hit_data = cell2mat(hit_data);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [591];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
alvaro_rothstein_2007.orfs = hit_strains;
alvaro_rothstein_2007.ph = hit_data_names;
alvaro_rothstein_2007.data = hit_data;
alvaro_rothstein_2007.dataset_ids = hit_data_ids;

%% Save

save('./alvaro_rothstein_2007.mat','alvaro_rothstein_2007');

%% Print out

fid = fopen('./alvaro_rothstein_2007.txt','w');
write_matrix_file(fid, alvaro_rothstein_2007.orfs, alvaro_rothstein_2007.ph, alvaro_rothstein_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(alvaro_rothstein_2007)
end

end

