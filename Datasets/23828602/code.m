%% Ding~Bakalinsky, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

ding_bakalinsky_2013.pmid = 23828602;

phenotypes = {'growth [pooled CFU]'};
treatments = {'acetic acid [122.5 mM]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/YSC1054Y.copy.xlsx', 'mat_alpha_obs');
tested_orfs = tested.raw(2:end,2);

inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];

tested_orfs = unique(upper(cleanOrf(tested_orfs)));

% Load data
fid = fopen('./raw_data/hits_genenames_R.txt');
hits_genenames_R = textscan(fid,'%s');
fclose(fid);

hits_genenames_R = hits_genenames_R{1};

hits_genenames_R = cleanGenename(hits_genenames_R);
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
return;

% Save data into database
dt = ding_bakalinsky_2013;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

