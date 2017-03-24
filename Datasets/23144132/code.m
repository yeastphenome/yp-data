%% Smith~Bakalinsky, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
smith_bakalinsky_2013.pmid = 23144132;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(smith_bakalinsky_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested (Same as Ding~Bakalinksy, 2013)
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/YSC1054Y.copy.xlsx', 'mat_alpha_obs');
tested_orfs = tested.raw(2:end,2);

% Eliminate all white spaces & capitalize
tested_orfs = clean_orf(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

% Get a unique list
tested_orfs = unique(tested_orfs);

%% Load hit data (resistant)
[FILENAMES, hits_orfs] = read_data('textscan', './raw_data/hits_orfs.txt', '%s');

% Eliminate all white spaces & capitalize
hits_orfs = clean_orf(hits_orfs);

% If in gene name form, transform into ORF name
hits_orfs = translate(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
hits_orfs(inds) = [];

% If possible, fix the problem (typos, omissions etc.)
hits_orfs(ismember(hits_orfs, {'YAL001'})) = {'YAL001C'};

% Get a unique list
hits_orfs = unique(hits_orfs);

% Make data
hits_scores = ones(length(hits_orfs),1);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [427];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);


% If the dataset is discrete/binary and the tested strains were provided separately:
smith_bakalinsky_2013.orfs = tested_orfs;
smith_bakalinsky_2013.ph = hit_data_names;
smith_bakalinsky_2013.data = zeros(length(smith_bakalinsky_2013.orfs),length(smith_bakalinsky_2013.ph));
smith_bakalinsky_2013.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, smith_bakalinsky_2013.orfs);
smith_bakalinsky_2013.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./smith_bakalinsky_2013.mat','smith_bakalinsky_2013');

%% Print out

fid = fopen('./smith_bakalinsky_2013.txt','w');
write_matrix_file(fid, smith_bakalinsky_2013.orfs, smith_bakalinsky_2013.ph, smith_bakalinsky_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(smith_bakalinsky_2013)
end

end
