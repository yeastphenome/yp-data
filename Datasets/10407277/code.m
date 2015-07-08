%% Rieger~Maundrell, 1999
% DATA = rieger_maundrell_1999
function FILENAMES = code()
FILENAMES = {};

rieger_maundrell_1999.source = {'manuscript PDF'};
rieger_maundrell_1999.downloaddate = {'2014-01-31'};
rieger_maundrell_1999.pmid = 10407277;

[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/rieger_maundrell_1999_data.xlsx', 'Data');

phenotypes = {'growth'};
treatments = data.raw(2:end,1);
phenotype_severity = {'HS','S','R'};
phenotype_severity_num = [-2,-1,1];

% Eliminate white spaces before/after ORF
data.raw(1,2:end) = cellfun(@strtrim, data.raw(1,2:end),'UniformOutput',0);

% Replace string scores with numbers
t = data.raw(2:end,2:end);
for i = 1 : length(phenotype_severity)
    inds = find(strcmp(phenotype_severity{i}, t));
    t(inds) = {phenotype_severity_num(i)};
end

rieger_maundrell_1999.orfs = upper(data.raw(1,2:end))';
rieger_maundrell_1999.data = cell2mat(t');
rieger_maundrell_1999.data(isnan(rieger_maundrell_1999.data)) = 0;
rieger_maundrell_1999.ph = [strcat(phenotypes{1}, '; ', treatments)];

% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(rieger_maundrell_1999.orfs, essential_genes);
rieger_maundrell_1999.orfs(ind1) = [];
rieger_maundrell_1999.data(ind1,:) = [];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'rieger_maundrell_1999.mat'],'rieger_maundrell_1999');
return;

% Save data into database
dt = rieger_maundrell_1999;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
adj_ix = [1 3 4 7 8 9 15 21 19 30 34 35 37 39 40 10 45 43 46 55 59 67 11 12 13 14 16 17 18 20 22 23 24 25 26 27 28 29 31 32 33 36 38 41 42 44 47 48 49 50 51 52 53 54 56 57 58 2 5 6 60 61 62 63 64 65 66];
[~,adj_ix] = sort(adj_ix);

datasets.names(database_ix(adj_ix),:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix(adj_ix)));

end

