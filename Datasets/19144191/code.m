%% Kemmer~Roberge, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));


addpath(genpath('../Yeast-Matlab-Utils/'));

FILENAMES = {};

kemmer_roberge_2009.pmid = 19144191;

phenotypes = {'growth (colony size)'};
treatments = {'dhMotC [60 uM]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/haploid set.xlsx', 'haploid set');
tested_orfs = tested.raw(6:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = clean_orf(tested_orfs);
tested_orfs(strcmp('YYKL138C', tested_orfs)) = {'YKL138C'};

tested_orfs = unique(tested_orfs);

% Load data
[FILENAMES{end+1}, hits_genenames] = read_data('textread','./raw_data/hits_genenames.txt', '%s');

hits_genenames = clean_genename(hits_genenames);
hits_orfs = translate(hits_genenames);
hits_orfs = unique(hits_orfs);

hits_scores = -ones(length(hits_orfs),1);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
hits_orfs(strcmp('YHR039C-A', hits_orfs)) = {'YHR039C-B'};

kemmer_roberge_2009.orfs = tested_orfs;
kemmer_roberge_2009.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
kemmer_roberge_2009.data(ind2,1) = hits_scores(ind1);

kemmer_roberge_2009.ph = [strcat(phenotypes, '; ', treatments)];

save('./kemmer_roberge_2009.mat','kemmer_roberge_2009');
return;

% Save data into database
dt = kemmer_roberge_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

fid = fopen('./kemmer_roberge_2009.txt','w');
write_matrix_file(fid, kemmer_roberge_2009.orfs, kemmer_roberge_2009.ph, kemmer_roberge_2009.data);
fclose(fid);

end

