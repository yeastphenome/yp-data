%% Mir~Cashikar, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

mir_cashikar_2009.pmid = 18936161;

phenotypes = {'growth [MIC]'};
treatments = {'heat stress (temperature [50ºC], duration [30 min])'};

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat_a.xlsx', 'mat_a_041902');
tested_orfs = tested.raw(3:end,2);


inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[FILENAMES, hits] = read_data('textscan', './raw_data/hits_orfs_scores.txt', '%s %d');
hits_orfs = upper(hits{1});
hits_scores = -hits{2};

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

mir_cashikar_2009.orfs = tested_orfs;
mir_cashikar_2009.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
mir_cashikar_2009.data(ind2,1) = hits_scores(ind1);

mir_cashikar_2009.ph = [strcat(phenotypes, '; ', treatments)];

save('./mir_cashikar_2009.mat','mir_cashikar_2009');

fid = fopen('./mir_cashikar_2009.txt','w');
write_matrix_file(fid, mir_cashikar_2009.orfs, mir_cashikar_2009.ph, mir_cashikar_2009.data);
fclose(fid);

end
