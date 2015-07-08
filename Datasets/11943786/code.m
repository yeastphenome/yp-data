%% Akache~Turcotte, 2002
% DATA = akache_turcotte_2002
function FILENAMES = code()
FILENAMES = {};

akache_turcotte_2002.source = {'manuscript PDF'};
akache_turcotte_2002.downloaddate = {'2014-02-10'};
akache_turcotte_2002.pmid = 11943786;

[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/akache_turcotte_2002_data.xlsx', 'Sheet1');

phenotypes = {'growth'};
treatments = data.raw(1,2:end)';
phenotype_severity = {'S';'SS';'R';'RR'};
phenotype_severity_num = [-1,-2,1,2];

% Eliminate white spaces before/after ORF
data.raw(2:end,1) = cellfun(@strtrim, data.raw(2:end,1),'UniformOutput',0);

% Replace string scores with numbers
t = data.raw(2:end,2:end);
for i = 1 : length(phenotype_severity)
    inds = find(strcmp(phenotype_severity{i}, t));
    t(inds) = {phenotype_severity_num(i)};
end
t(cellfun(@isnan,t))={0};

akache_turcotte_2002.orfs = upper(data.raw(2:end,1));
akache_turcotte_2002.data = cell2mat(t);
akache_turcotte_2002.ph = [strcat(phenotypes{1}, '; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(akache_turcotte_2002.orfs, essential_genes);
akache_turcotte_2002.orfs(ind1) = [];
akache_turcotte_2002.data(ind1,:) = [];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'akache_turcotte_2002.mat'],'akache_turcotte_2002');
return;

% Save data into database
dt = akache_turcotte_2002;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

