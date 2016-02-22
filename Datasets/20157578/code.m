%% Postma~Ralser, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

postma_ralser_2009.pmid = 20157578;

phenotypes = {'growth [spot assay]'};
treatments = {'hibernation'};

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat_a_obs_v2.0.xlsx', 'Mat_a_obs_v2.0.txt');
tested_orfs = tested.raw(2:end,1);


inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[FILENAMES, hits_orfs] = read_data('fopen', './raw_data/hits_orfs.txt', '%s');

hits_orfs = unique(hits_orfs{1});
hits_scores = zeros(length(hits_orfs),1)+1;

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustment
tested_orfs = [tested_orfs; missing];   % 1 ORF added

postma_ralser_2009.orfs = tested_orfs;
postma_ralser_2009.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
postma_ralser_2009.data(ind2,1) = hits_scores(ind1);

postma_ralser_2009.ph = [strcat(phenotypes, '; ', treatments)];

save('./postma_ralser_2009.mat','postma_ralser_2009');

fid = fopen('./postma_ralser_2009.txt','w');
write_matrix_file(fid, postma_ralser_2009.orfs, postma_ralser_2009.ph, postma_ralser_2009.data);
fclose(fid);

end
