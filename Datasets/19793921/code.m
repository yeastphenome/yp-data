%% Kanki~Klionsky, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
kanki_klionsky_2009.pmid = 19793921;

phenotypes = {'mitophagy'};
treatments = {'YPL'};

% Load data
[FILENAMES{end+1}, data.raw] = readdata('xlsread','./raw_data/TableS2-2.xlsx', 'Sheet1');
hits_orfs = data.raw(4:end,1);
hits_scores = data.raw(4:end,2);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_scores(inds) = [];
hits_orfs = upper(hits_orfs);

inds = strcmp('x', hits_scores);
hits_scores(inds) = {NaN};
inds = strcmp('No', hits_scores);
hits_scores(inds) = {-1};
inds = strcmp('Yes', hits_scores);
hits_scores(inds) = {1};

hits_scores = cell2mat(hits_scores);

hits_orfs_u = unique(hits_orfs);
hits_scores_u = zeros(length(hits_orfs_u),1)+NaN;
for i = 1 : length(hits_orfs_u)
inds = find(strcmp(hits_orfs_u{i}, hits_orfs));
if length(inds) == 1
hits_scores_u(i) = hits_scores(inds);
elseif length(inds) == 2
tmp = prod(hits_scores(inds));
if isnan(tmp) || tmp < 0    % If at least one of the values is NaN or the values are of opposite sign, set them to NaN
hits_scores_u(i) = NaN;
else
hits_scores_u(i) = hits_scores(inds(1));    % If the values are of the same sign, pick the first number (they are the same)
end
else
fprintf('%s\t', hits_orfs_u{i});
for j = 1 : length(inds)
fprintf('%d,', hits_scores(inds(j)));
end
fprintf('\n');
end
end

kanki_klionsky_2009.orfs = hits_orfs_u;
kanki_klionsky_2009.data = hits_scores_u;

kanki_klionsky_2009.ph = strcat(phenotypes, '; ', treatments);

save('./kanki_klionsky_2009.mat','kanki_klionsky_2009');
return;

% Save data into database
dt = kanki_klionsky_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

