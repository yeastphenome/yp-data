%% Breslow~Weissman, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
breslow_weissman_2008.pmid = 18622397;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(breslow_weissman_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/nmeth.1234-S5.xlsx', 'Deletion growth rates');

% Get indices of the data columns
ind_data = find(strcmp('Median Growth Rate', data.raw(1,:)));

hit_strains = data.raw(:,2);
hit_data = data.raw(:,ind_data);

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds)); 

hit_strains(inds) = [];
hit_data(inds,:) = [];

% Make sure all the data are numbers
inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {0};
hit_data = cell2mat(hit_data);

% Average data for identical ORFs that appear multiple times
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [26];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
breslow_weissman_2008.orfs = hit_strains;
breslow_weissman_2008.ph = hit_data_names;
breslow_weissman_2008.data = hit_data;
breslow_weissman_2008.dataset_ids = hit_data_ids;

%% Save

save('./breslow_weissman_2008.mat','breslow_weissman_2008');

%% Print out

fid = fopen('./breslow_weissman_2008.txt','w');
write_matrix_file(fid, breslow_weissman_2008.orfs, breslow_weissman_2008.ph, breslow_weissman_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(breslow_weissman_2008)
end

end
