%% Orij~Smits, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

orij_smits_2012.pmid = 23021432;

phenotypes = {'cytosolic pH'};
treatments = {'standard'};

% Load hits
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/pH Screen raw.xlsx', 'initial screens');
hits_orfs = data.raw(2:end,1);
hits_scores = data.raw(2:end,3:4);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds,:) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_scores(inds,:) = [];

hits_orfs = upper(hits_orfs);
hits_scores = cell2mat(hits_scores);
hits_scores = nanmean(hits_scores, 2);

[hits_orfs, hits_scores] = grpstats(hits_scores,hits_orfs,{'gname','mean'});

orij_smits_2012.orfs = hits_orfs;
orij_smits_2012.data = hits_scores;

orij_smits_2012.ph = [strcat(phenotypes, '; ', treatments)];

save('./orij_smits_2012.mat','orij_smits_2012');

fid = fopen('./orij_smits_2012.txt','w');
write_matrix_file(fid, orij_smits_2012.orfs, orij_smits_2012.ph, orij_smits_2012.data);
fclose(fid);

end
