%% Merz~Westermann, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
merz_westermann_2009.pmid = 19751518;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(merz_westermann_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/gb-2009-10-9-r95-s1.xlsx');

% Get the list of ORFs and the correponding data 
hits_orfs = data.raw(3:end,1);

% Eliminate all white spaces & capitalize
hits_orfs = clean_orf(hits_orfs);

% If in gene name form, transform into ORF name
[hits_orfs, translated, ambiguous] = translate(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
hits_orfs(inds) = [];

% Make unique
hits_orfs = unique(hits_orfs);

% Make the data
hits_data = -ones(size(hits_orfs));

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [158];

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/pet-Screen.xlsx', 'mat_alpha_obs');

% Load tested strains
tested_orfs = tested.raw(2:end,2);

% Remove numeric/empty indices
inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

% Eliminate all white spaces & capitalize
tested_orfs = clean_orf(tested_orfs);

% If in gene name form, transform into ORF name
[tested_orfs, translated, ambiguous] = translate(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

% Make unique
tested_orfs = unique(tested_orfs);

% Make sure the that all the hits are part of the tested set
[missing, ix] = setdiff(hits_orfs, tested_orfs);
tested_orfs = [tested_orfs; missing];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
merz_westermann_2009.orfs = tested_orfs;
merz_westermann_2009.ph = hit_data_names;
merz_westermann_2009.data = zeros(length(merz_westermann_2009.orfs),length(merz_westermann_2009.ph));
merz_westermann_2009.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, merz_westermann_2009.orfs);
merz_westermann_2009.data(ind2) = hits_data(ind1);

%% Save

save('./merz_westermann_2009.mat','merz_westermann_2009');

%% Print out

fid = fopen('./merz_westermann_2009.txt','w');
write_matrix_file(fid, merz_westermann_2009.orfs, merz_westermann_2009.ph, merz_westermann_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(merz_westermann_2009)
end

end
