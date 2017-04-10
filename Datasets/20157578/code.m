%% Postma~Ralser, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
postma_ralser_2009.pmid = 20157578;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(postma_ralser_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat_a_obs_v2.0.xlsx', 'Mat_a_obs_v2.0.txt');

% Get the list of ORFs and the correponding data 
tested_orfs = tested.raw(2:end,1);

% Eliminate all white spaces & capitalize
tested_orfs = clean_orf(tested_orfs);

% Eliminate any empty indices
inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];

% If in gene name form, transform into ORF name
[tested_orfs, translated, ambiguous] = translate(tested_orfs);

% Find anything that doesn't look like an ORF
tested_orfs(ismember(tested_orfs, {'YLR287-A'})) = {'YLR287C-A'};
inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

% Make sure all ORFs are unique
tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES, hits_orfs] = read_data('textscan', './raw_data/hits_orfs.txt', '%s');

% Eliminate all white spaces & capitalize
hits_orfs = clean_orf(hits_orfs);

% If in gene name form, transform into ORF name
[hits_orfs, translated, ambiguous] = translate(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
hits_orfs(inds) = [];

% Make sure the ORFs are unique
hits_orfs = unique(hits_orfs);

% Get the data itself
hits_scores = ones(length(hits_orfs),1);

% Make sure the that all the hits are part of the tested set
[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustment
tested_orfs = [tested_orfs; missing];   % 1 ORF added

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [156];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
postma_ralser_2009.orfs = tested_orfs;
postma_ralser_2009.ph = hit_data_names;
postma_ralser_2009.data = zeros(length(postma_ralser_2009.orfs),length(postma_ralser_2009.ph));
postma_ralser_2009.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, postma_ralser_2009.orfs);
postma_ralser_2009.data(ind2) = hits_scores(ind1);

%% Save

save('./postma_ralser_2009.mat','postma_ralser_2009');

%% Print out

fid = fopen('./postma_ralser_2009.txt','w');
write_matrix_file(fid, postma_ralser_2009.orfs, postma_ralser_2009.ph, postma_ralser_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(postma_ralser_2009)
end

end
