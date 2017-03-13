%% Gaupel~Tenniswood, 2014
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
gaupel_tenniswood_2004.pmid = 24968945;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(gaupel_tenniswood_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data
[FILENAMES{end+1}, data_sens.raw] = read_data('xlsread','./raw_data/1471-2164-15-528-s1.xlsx', 'Top sensitive strains');

hits_orfs1 = data_sens.raw(2:end,1);
hits_data1 = -cell2mat(data_sens.raw(2:end,3));

hits_orfs1 = strtrim(upper(hits_orfs1));

inds = find(~strncmp('Y', hits_orfs1,1));
hits_orfs1(inds) = [];
hits_data1(inds,:) = [];

[t,t2] = grpstats(hits_data1,hits_orfs1,{'gname','mean'});
hits_orfs1 = t;
hits_data1 = t2;

% Load data (2)
[FILENAMES{end+1}, data_res.raw] = read_data('xlsread','./raw_data/1471-2164-15-528-s1.xlsx', 'Top resistant strains');

hits_orfs2 = data_res.raw(2:end,1);
hits_data2 = -cell2mat(data_res.raw(2:end,3));

hits_orfs2 = strtrim(upper(hits_orfs2));

inds = find(~strncmp('Y', hits_orfs2,1));
hits_orfs2(inds) = [];
hits_data2(inds,:) = [];

[t,t2] = grpstats(hits_data2,hits_orfs2,{'gname','mean'});
hits_orfs2 = t;
hits_data2 = t2;

% Eliminate overlap between sensitive and resistant strains?
[~,ind1,ind2] = intersect(hits_orfs1, hits_orfs2);
hits_orfs1(ind1) = []; hits_data1(ind1,:) = [];
hits_orfs2(ind2) = []; hits_data2(ind2,:) = [];

% Combine part1 and part2
hits_orfs = [hits_orfs1; hits_orfs2];
hits_data = [hits_data1; hits_data2];

%% Tested strains

% Load tested strains
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/CompleteDeletionLibrary.xlsx');
tested_orfs = tested.raw(2:end,1);

% Remove all indices without a name
inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

% Take the unique set
tested_orfs = unique(strtrim(upper(tested_orfs)));

% Remove non-relevant indices
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Make sure that all the hits are part of the tested set
[missing, ix] = setdiff(hits_orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [569];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
gaupel_tenniswood_2004.orfs = tested_orfs;
gaupel_tenniswood_2004.ph = hit_data_names;
gaupel_tenniswood_2004.data = zeros(length(gaupel_tenniswood_2004.orfs),length(gaupel_tenniswood_2004.ph));
gaupel_tenniswood_2004.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
gaupel_tenniswood_2004.data(ind2,:) = hits_data(ind1,:);

%% Save

save('./gaupel_tenniswood_2004.mat','gaupel_tenniswood_2004');

%% Print out

fid = fopen('./gaupel_tenniswood_2004.txt','w');
write_matrix_file(fid, gaupel_tenniswood_2004.orfs, gaupel_tenniswood_2004.ph, gaupel_tenniswood_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(deutschbauer_giaevere_2005)
end

end
