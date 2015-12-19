%% Sambade~Kane, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

sambade_kane_2005.pmid = 15937126;

phenotypes = {'growth (colony size)'};
treatments = {'pH [7.5] CaCl2 [60 mM]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/ResGen 384 well set 14 plates.xlsx', 'ResGen MATa -384');
tested_orfs = tested.raw(4:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs(ismember(tested_orfs,{'YLR287-A'})) = {'YLR287C-A'};

tested_orfs = unique(upper(cleanOrf(tested_orfs)));

% Load data
[FILENAMES{end+1}, hits_genenames] = dataread('textread','./raw_data/hits_genenames.txt', '%s');
hits_orfs = translate(hits_genenames);

hits_scores = -ones(length(hits_orfs),1);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
hits_orfs(strcmpi('YHR039C-A', hits_orfs)) = {'YHR039C-B'};

[missing, ix] = setdiff(hits_orfs, tested_orfs);
tested_orfs = [tested_orfs; missing];   % 2 orfs to be added

sambade_kane_2005.orfs = tested_orfs;
sambade_kane_2005.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
sambade_kane_2005.data(ind2,1) = hits_scores(ind1);

sambade_kane_2005.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'sambade_kane_2005.mat'],'sambade_kane_2005');
return;

% Save data into database
dt = sambade_kane_2005;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

