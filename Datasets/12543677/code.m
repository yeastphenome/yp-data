%% Blackburn~Avery, 2003
% DATA = blackburn_avery_2003
function FILENAMES = code()
FILENAMES = {};
% TESTED = not available

blackburn_avery_2003.source = {'manuscript PDF'};
blackburn_avery_2003.downloaddate = {'2014-02-12'};
blackburn_avery_2003.pmid = 12543677;

[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/blackburn_avery_2003_data.xlsx', 'data.txt');

phenotypes = {'growth (MIC)'};
treatments = data.raw(1,2:8)';

% Eliminate white spaces before/after ORF
data.raw(2:end,1) = cellfun(@strtrim, data.raw(2:end,1),'UniformOutput',0);

% Replace 'Inf' with Inf
t = data.raw(2:end,2:end);
t(~cellfun(@isnumeric, t)) = {Inf};

blackburn_avery_2003.orfs = upper(data.raw(2:end,1));
blackburn_avery_2003.data = cell2mat(t);

blackburn_avery_2003.ph = [strcat(phenotypes{1}, '; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(blackburn_avery_2003.orfs, essential_genes);
blackburn_avery_2003.orfs(ind1) = [];
blackburn_avery_2003.data(ind1,:) = [];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'blackburn_avery_2003.mat'],'blackburn_avery_2003');
return;

% Save data into database
dt = blackburn_avery_2003;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


end

