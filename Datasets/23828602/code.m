%% Ding~Bakalinsky, 2013
function FILENAMES = code()
FILENAMES = {};

ding_bakalinsky_2013.pmid = 23828602;

phenotypes = {'growth [pooled CFU]'};
treatments = {'acetic acid [122.5 mM]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/YSC1054Y.copy.xlsx', 'mat_alpha_obs');
tested_orfs = tested.raw(2:end,2);

inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
fid = fopen('raw_data/hits_genenames_R.txt');
hits_genenames_R = textscan(fid,'%s');
fclose(fid);

hits_genenames_R = hits_genenames_R{1};

hits_genenames_R = cellfun(@strtrim, hits_genenames_R,'UniformOutput',0);
hits_genenames_R = unique(upper(hits_genenames_R));
hits_orfs_R = unique(genename2orf(hits_genenames_R,'noannot'));
hits_scores_R = zeros(length(hits_orfs_R),1)+1;

inds = find(~strncmp('Y', hits_orfs_R,1));
hits_orfs_R(inds) = [];
hits_scores_R(inds) = [];

[missing, ix] = setdiff(hits_orfs_R, tested_orfs);

ding_bakalinsky_2013.orfs = tested_orfs;
ding_bakalinsky_2013.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs_R, tested_orfs);
ding_bakalinsky_2013.data(ind2,1) = hits_scores_R(ind1);

ding_bakalinsky_2013.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'ding_bakalinsky_2013.mat'],'ding_bakalinsky_2013');
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

