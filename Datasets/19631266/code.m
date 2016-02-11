%% Zhou~Costa, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
% TESTED = not available
zhou_costa_2009.pmid = 19631266;

phenotypes = {'growth [spot assay]'};
treatments = {'NaAsO2 [0.075-1 mM]'};

% Load hits
[FILENAMES{end+1}, hits.raw] = readdata('xlsread','./raw_data/table.xlsx', 'table.csv');
hits_orfs = hits.raw(4:end,1);
hits_scores = hits.raw(4:end,4);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = upper(hits_orfs);
hits_scores(strcmp('Above 5', hits_scores)) = {5};
hits_scores = cell2mat(hits_scores);

[hits_orfs, hits_scores] = grpstats(hits_scores,hits_orfs,{'gname','mean'});

wt = 4.47;  %IC50 from paper
hits_scores = hits_scores - wt; % negative scores = more sensitive than WT, and viceversa.


zhou_costa_2009.orfs = hits_orfs;
zhou_costa_2009.data = hits_scores;

zhou_costa_2009.ph = [strcat(phenotypes, '; ', treatments)];

save('./zhou_costa_2009.mat','zhou_costa_2009');
return;

% Save data into database
dt = zhou_costa_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

