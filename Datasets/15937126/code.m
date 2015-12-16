%% Sambade~Kane, 2005
function FILENAMES = code()
FILENAMES = {};

sambade_kane_2005.pmid = 15937126;

phenotypes = {'growth (colony size)'};
treatments = {'pH [7.5] CaCl2 [60 mM]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/ResGen 384 well set 14 plates.xlsx', 'ResGen MATa -384');
tested_orfs = tested.raw(4:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[FILENAMES{end+1}, hits_genenames] = dataread('textread','./raw_data/hits_genenames.txt', '%s');

hits_genenames = upper(hits_genenames);
hits_orfs = genename2orf(hits_genenames,'noannot');

hits_orfs(strcmpi('ada3', hits_orfs)) = {'YDR176W'};
hits_orfs(strcmpi('cak1', hits_orfs)) = {'YFL029C'};
hits_orfs(strcmpi('cwh36', hits_orfs)) = {'YCL007C'};
hits_orfs(strcmpi('rcs1', hits_orfs)) = {'YGL071W'};
hits_orfs(strcmpi('rmd7', hits_orfs)) = {'YER083C'};
hits_orfs(strcmpi('vma1', hits_orfs)) = {'YDL185W'};
hits_orfs(strcmpi('vma12', hits_orfs)) = {'YKL119C'};
hits_orfs(strcmpi('vma16', hits_orfs)) = {'YHR026W'};
hits_orfs(strcmpi('vma3', hits_orfs)) = {'YEL027W'};

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

hits_scores = -ones(length(hits_orfs),1);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = upper(hits_orfs);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
hits_orfs(strcmpi('YHR039C-A', hits_orfs)) = {'YHR039C-B'};

[missing, ix] = setdiff(hits_orfs, tested_orfs);
tested_orfs = [tested_orfs; missing];   % 1 orfs to be added

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

