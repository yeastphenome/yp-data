%% Teixeira~Sa-Correia, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
teixeira_sa_correia_2010.pmid = 20210661;

phenotypes = {'growth'};
treatments = {'Glu, 30%'};

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
[FILENAMES{end+1}, hits_genenames] = read_data('textread','./raw_data/hits_genenames.txt', '%s');
inds = find(cellfun(@isempty, hits_genenames) | cellfun(@isnumeric, hits_genenames));
hits_genenames(inds) = [];
hits_genenames = strtrim(hits_genenames);

hits_orfs = translate(hits_genenames);
hits_orfs(~is_orf(hits_orfs)) = [];

hits_orfs = unique(hits_orfs);

[missing, ix] = setdiff(hits_orfs, tested_orfs2);
hits_orfs(strcmp('YML009W-B',hits_orfs)) = {'YML010W-A'};

teixeira_sa_correia_2010.orfs = tested_orfs;
teixeira_sa_correia_2010.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
teixeira_sa_correia_2010.data(ind2) = -1;

teixeira_sa_correia_2010.ph = strcat(phenotypes, '; ', treatments);

save('./teixeira_sa_correia_2010.mat','teixeira_sa_correia_2010');
return;

% Save data into database
dt = teixeira_sa_correia_2010;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

fid = fopen('./teixeira_sa_correia_2010.txt','w');
write_matrix_file(fid, teixeira_sa_correia_2010.orfs, teixeira_sa_correia_2010.ph, teixeira_sa_correia_2010.data);
fclose(fid);

end

