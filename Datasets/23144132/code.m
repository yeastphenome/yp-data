%% Smith~Bakalinsky, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
smith_bakalinsky_2013.pmid = 23144132;

phenotypes = {'growth [CFU]'};
treatments = {'AuNP [10-100 ug/ml]'};

% Load tested (Same as Ding~Bakalinksy, 2013)
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/YSC1054Y.copy.xlsx', 'mat_alpha_obs');
tested_orfs = tested.raw(2:end,2);

inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load resistant
fid = fopen('./raw_data/hits_orfs.txt');
hits = textscan(fid,'%s');
hits_orfs = hits{1};
fclose(fid);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_orfs = unique(upper(hits_orfs));

hits_scores = ones(length(hits_orfs),1);


smith_bakalinsky_2013.orfs = tested_orfs;
smith_bakalinsky_2013.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
smith_bakalinsky_2013.data(ind2,1) = hits_scores(ind1);

smith_bakalinsky_2013.ph = [strcat(phenotypes, '; ', treatments)];

save('./smith_bakalinsky_2013.mat','smith_bakalinsky_2013');
return;

% Save data into database
dt = smith_bakalinsky_2013;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


fid = fopen('./smith_bakalinsky_2013.txt','w');
write_matrix_file(fid, smith_bakalinsky_2013.orfs, smith_bakalinsky_2013.ph, smith_bakalinsky_2013.data);
fclose(fid);

end

