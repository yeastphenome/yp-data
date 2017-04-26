%% Brett~Rao, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
brett_rao_2011.pmid = 21423800;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(brett_rao_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/journal.pone.0017619.s003.xlsx', 'Unsorted Data');

% Get the list of ORFs and the correponding data 
hit_strains = data.raw(:,2);

% Eliminate white spaces before/after ORF
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% Eliminate anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
data.raw(inds,:) = [];

% Make sure all the data are numbers
hit_data = data.raw(:, [4:6 8:10]);
inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {NaN};
hit_data = cell2mat(hit_data);

% Average data for identical ORFs that appear multiple times
[hit_strains,hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [59; 60; 61; 529; 530; 531];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
brett_rao_2011.orfs = hit_strains;
brett_rao_2011.ph = hit_data_names;
brett_rao_2011.data = hit_data;
brett_rao_2011.dataset_ids = hit_data_ids;

%% Save

save('./brett_rao_2011.mat','brett_rao_2011');

%% Print out

fid = fopen('./brett_rao_2011.txt','w');
write_matrix_file(fid, brett_rao_2011.orfs, brett_rao_2011.ph, brett_rao_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(brett_rao_2011)
end

end
