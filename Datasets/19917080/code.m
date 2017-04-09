%% Arita~Costa, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
arita_costa_2009.pmid = 19917080;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(arita_costa_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, hits.raw] = read_data('xlsread','./raw_data/1471-2164-10-524-S1.XLS', 'Sheet1');

% Get the list of ORFs and the correponding data 
hits_sensitive_orfs = hits.raw(7:end,1);
hits_resistant_orfs = hits.raw(7:end,2);

% Eliminate all numeric indices
inds = cellfun(@isnumeric, hits_resistant_orfs);
hits_resistant_orfs(inds) = [];

% Eliminate all white spaces & capitalize
hits_resistant_orfs = clean_orf(hits_resistant_orfs);

% If in gene name form, transform into ORF name
[hits_resistant_orfs, translated, ambiguous] = translate(hits_resistant_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_resistant_orfs));
hits_resistant_orfs(inds) = [];

% Make sure the list is unique
hits_resistant_orfs = unique(upper(hits_resistant_orfs));

% Make the data
hits_resistant_scores = ones(length(hits_resistant_orfs),1);

% Do the same for the sensitive dataset
% Get rid of all the numeric data
inds = cellfun(@isnumeric, hits_sensitive_orfs);
hits_sensitive_orfs(inds) = [];

% Eliminate all white spaces & capitalize
hits_sensitive_orfs = clean_orf(hits_sensitive_orfs);

% If in gene name form, transform into ORF name
[hits_sensitive_orfs, translated, ambiguous] = translate(hits_sensitive_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_sensitive_orfs));
hits_sensitive_orfs(inds) = [];

% Make sure the list is unique
hits_sensitive_orfs = unique(upper(hits_sensitive_orfs));

% Make the data
hits_sensitive_scores = -ones(length(hits_sensitive_orfs),1);

% Combine to make tested
tested_orfs = [hits_resistant_orfs; hits_sensitive_orfs];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [163];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
arita_costa_2009.orfs = tested_orfs;
arita_costa_2009.ph = hit_data_names;
arita_costa_2009.data = zeros(length(arita_costa_2009.orfs),length(arita_costa_2009.ph));
arita_costa_2009.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_resistant_orfs, arita_costa_2009.orfs);
arita_costa_2009.data(ind2) = hits_resistant_scores(ind1);
[~,ind1,ind2] = intersect(hits_sensitive_orfs, arita_costa_2009.orfs);
arita_costa_2009.data(ind2) = hits_sensitive_scores(ind1);

%% Save

save('./arita_costa_2009.mat','arita_costa_2009');

%% Print out

fid = fopen('./arita_costa_2009.txt','w');
write_matrix_file(fid, arita_costa_2009.orfs, arita_costa_2009.ph, arita_costa_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(arita_costa_2009)
end

end
