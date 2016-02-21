%% Ding~Bakalinsky, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

ding_bakalinsky_2013.pmid = 23828602;

phenotypes = {'growth [pooled CFU]'};
treatments = {'acetic acid [122.5 mM]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/YSC1054Y.copy.xlsx', 'mat_alpha_obs');
tested_orfs = tested.raw(2:end,2);

inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];

tested_orfs = unique(upper(clean_orf(tested_orfs)));

% Load data
fid = fopen('./raw_data/hits_genenames_R.txt');
hits_genenames_R = textscan(fid,'%s');
fclose(fid);

hits_genenames_R = hits_genenames_R{1};

hits_genenames_R = clean_genename(hits_genenames_R);
hits_orfs_R = translate(hits_genenames_R);
hits_orfs_R = unique(hits_orfs_R);
hits_scores_R = zeros(length(hits_orfs_R),1)+1;

[missing, ix] = setdiff(hits_orfs_R, tested_orfs);

ding_bakalinsky_2013.orfs = tested_orfs;
ding_bakalinsky_2013.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs_R, tested_orfs);
ding_bakalinsky_2013.data(ind2,1) = hits_scores_R(ind1);

ding_bakalinsky_2013.ph = [strcat(phenotypes, '; ', treatments)];

save('./ding_bakalinsky_2013.mat','ding_bakalinsky_2013');

fid = fopen('./ding_bakalinsky_2013.txt','w');
write_matrix_file(fid, ding_bakalinsky_2013.orfs, ding_bakalinsky_2013.ph, ding_bakalinsky_2013.data);
fclose(fid);

end
