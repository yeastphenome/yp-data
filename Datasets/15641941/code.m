%% Outten~Culotta, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

outten_culotta_2005.source = {'http://www.biochemj.org/bj/388/bj3880093add.pdf'};
outten_culotta_2005.downloaddate = {'2014-03-10'};
outten_culotta_2005.pmid = 15641941;

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/outten_culotta_2005.xlsx', 'Sheet1');

hits_orfs = data.raw(:,1);

% Eliminate white spaces before/after ORF
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

hits_scores = cell2mat(data.raw(:,2));

% Adjust scores such that -1 = weakest phenotype, -4 = strongest phenotype.
hits_scores = hits_scores - 5;

% The 6 hits with no score (not sure why), set to 0
hits_scores(isnan(hits_scores)) = 0;

phenotypes = {'growth'};
treatments = {'hyperoxia'};


% Load tested genes
[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/Yeast Knockout -BY4741.xlsx', 'mat_a_060701.txt');

tested_orfs = data.raw(2:end,2);

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

tested_orfs = unique(upper(tested_orfs));

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% From supplement: MED2 (YDL005C) not available in Mat-a background, was used Mat-alpha
% instead. So it should be added it to the list of tested
tested_orfs = [tested_orfs; {'YDL005C'}];

% YHR039C-A is not in the list of tested set, but it has the alias
% YHR039C-B, which is in the tested set. So, these are likely to be the same gene, so YHR039C-A should be renamed.
hits_orfs(strcmp('YHR039C-A', hits_orfs)) = {'YHR039C-B'};

% Create dataset
outten_culotta_2005.orfs = tested_orfs;
outten_culotta_2005.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(outten_culotta_2005.orfs, hits_orfs);
outten_culotta_2005.data(ind1,:) = hits_scores(ind2,:);

outten_culotta_2005.ph = [strcat(phenotypes, '; ', treatments)];

save('./outten_culotta_2005.mat','outten_culotta_2005');
return;

% Save data into database
dt = outten_culotta_2005;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

