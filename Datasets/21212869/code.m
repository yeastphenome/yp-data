%% Bleackley~MacGillivray, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

bleackley_macgillivray_2011.pmid = 21212869;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(bleackley_macgillivray_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data 

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/metallomicsbleackley raw data.xls', 'rawdata');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data.raw(:,1);

% Get the data itself
hit_data = data.raw(:, 3:2:13); % if the dataset is discrete or binary
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Eliminate anything that is numeric
inds = find(cellfun(@isnumeric, hit_strains) | cellfun(@isempty, hit_strains));
hit_strains(inds) = [];
hit_data(inds, :) = [];

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds, :) = [];

% Make sure all the data are numbers
inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {NaN};
hit_data = cell2mat(hit_data);

% Average data for identical ORFs that appear multiple times
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [20; 21; 22; 23; 24; 25];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
bleackley_macgillivray_2011.orfs = hit_strains;
bleackley_macgillivray_2011.ph = hit_data_names;
bleackley_macgillivray_2011.data = hit_data;
bleackley_macgillivray_2011.dataset_ids = hit_data_ids;

%% Save

save('./bleackley_macgillivray_2011.mat','bleackley_macgillivray_2011');

%% Print out

fid = fopen('./bleackley_macgillivray_2011.txt','w');
write_matrix_file(fid, bleackley_macgillivray_2011.orfs, bleackley_macgillivray_2011.ph, bleackley_macgillivray_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(bleackley_macgillivray_2011)
end

end
