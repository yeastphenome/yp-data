%% Parsons~Boone, 2004
function FILENAMES = code()
FILENAMES = {};
% TESTED = not available

parsons_boone_2004.source = {'http://www.nature.com/nbt/journal/v22/n1/extref/nbt919-S2.xls'};
parsons_boone_2004.downloaddate = {'2014-02-20'};
parsons_boone_2004.pmid = 14661025;

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/nbt919-S2.xlsx', 'Sheet1');

inds = find(strcmp('ORF', data.raw(:,1)));

phenotypes = {'growth (colony size)'};
treatments = data.raw(inds,3:end)';

% Eliminate white spaces before/after ORF
data.raw(inds+1:end,1) = cellfun(@strtrim, data.raw(inds+1:end,1),'UniformOutput',0);

% Replace 'Inf' with Inf
t = data.raw(inds+1:end,3:end);
t(cellfun(@isnan, t)) = {0};

parsons_boone_2004.orfs = upper(data.raw(inds+1:end,1));
parsons_boone_2004.data = cell2mat(t);

% Flip the sign of the values, such that negative = sensitive, positive =
% resistant
parsons_boone_2004.data = -parsons_boone_2004.data;

parsons_boone_2004.ph = [strcat(phenotypes{1}, '; ', treatments)];


% Eliminate the essential genes
% load essential_genes_100908;
% [t,ind1,ind2] = intersect(parsons_boone_2004.orfs, essential_genes);
% parsons_boone_2004.orfs(ind1) = [];
% parsons_boone_2004.data(ind1,:) = [];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'parsons_boone_2004.mat'],'parsons_boone_2004');
return;

% Save data into database
dt = parsons_boone_2004;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
adj_ix = [4 6 8 1 2 3 5 7 9 10 11 12];
[~,adj_ix] = sort(adj_ix);

datasets.names(database_ix(adj_ix),:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix(adj_ix)));

end

