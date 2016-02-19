%% Garay~DeLuna, 2014
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

garay_deluna_2014.pmid = 24586198;

phenotypes = {'chronological life span'};
treatments = {'standard'};

% Load hits
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/journal.pgen.1004168.s013.xlsx', 'Genome-wide CLS screen');
hits_orfs = data.raw(2:end,1);
hits_scores = data.raw(2:end,3);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds,:) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_scores(inds,:) = [];

hits_orfs = upper(hits_orfs);
hits_scores = cell2mat(hits_scores);

[hits_orfs, hits_scores] = grpstats(hits_scores,hits_orfs,{'gname','mean'});

garay_deluna_2014.orfs = hits_orfs;
garay_deluna_2014.data = hits_scores;

garay_deluna_2014.ph = [strcat(phenotypes, '; ', treatments)];

save('./garay_deluna_2014.mat','garay_deluna_2014');
return;

% Save data into database
dt = garay_deluna_2014;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

fid = fopen('./garay_deluna_2014.txt','w');
write_matrix_file(fid, garay_deluna_2014.orfs, garay_deluna_2014.ph, garay_deluna_2014.data);
fclose(fid);

end

