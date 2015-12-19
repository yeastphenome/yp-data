%% Huang~Paulovich, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

huang_paulovich_2013.pmid = 23382077;

phenotypes = {'growth'};
treatments = {'MMS'};

% Load plate maps
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/Mat_a_obs_v4 0.xls', 'DATA');
tested_orfs = tested.raw(2:end,2);
tested_orfs(cellfun(@isnumeric, tested_orfs)) = [];

tested_orfs = regexprep(tested_orfs, '\W','');
tested_orfs = unique(upper(tested_orfs));
tested_orfs(ismember(tested_orfs,{'YLR287A'})) = {'YLR287C-A'};
tested_orfs(~isorf(tested_orfs)) = [];


% Load data
[FILENAMES{end+1}, hits_genenames] = dataread('textread','./raw_data/huang_paulovich_2013_hits.txt', '%s');

hits_genenames = regexprep(hits_genenames, '\W','');
hits_orfs = translate(hits_genenames);
hits_orfs(~isorf(hits_orfs)) = [];

hits_scores = -ones(length(hits_orfs),1);

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
huang_paulovich_2013.orfs = tested_orfs;
huang_paulovich_2013.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(huang_paulovich_2013.orfs, hits_orfs);
huang_paulovich_2013.data(ind1,:) = hits_scores(ind2,:);

huang_paulovich_2013.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'huang_paulovich_2013.mat'],'huang_paulovich_2013');
return;

% Save data into database
dt = huang_paulovich_2013;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


end

