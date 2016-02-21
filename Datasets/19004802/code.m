%% Ralser~Lehrach, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ralser_lehrach_2008.pmid = 19004802;

phenotypes = {'growth'};
treatments = {'2-DG'};

% Load data
[FILENAMES{end+1}, hits_orfs] = read_data('textread','./raw_data/hits_orfs.txt', '%s');
hits_orfs = unique(strtrim(upper(hits_orfs)));

hits_data = ones(size(hits_orfs));

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat_a_obs_v2.0.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

ralser_lehrach_2008.orfs = tested_orfs;
ralser_lehrach_2008.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
ralser_lehrach_2008.data(ind2,:) = hits_data(ind1,:);

ralser_lehrach_2008.ph = strcat(phenotypes, '; ', treatments);

save('./ralser_lehrach_2008.mat','ralser_lehrach_2008');

fid = fopen('./ralser_lehrach_2008.txt','w');
write_matrix_file(fid, ralser_lehrach_2008.orfs, ralser_lehrach_2008.ph, ralser_lehrach_2008.data);
fclose(fid);

end
