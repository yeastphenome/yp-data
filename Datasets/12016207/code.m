%% Desmoucelles~Daignan-Fornier, 2002
function FILENAMES = code()
FILENAMES = {};
% NOTE = 1) typo in data: YKR065W should be YKR065C 2) YCR002C is in the
% result, but not in tested.

desmoucelles_daignan_fornier_2002.source = {'manuscript PDF'};
desmoucelles_daignan_fornier_2002.downloaddate = {'2014-02-10'};
desmoucelles_daignan_fornier_2002.pmid = 12016207;

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/desmoucelles_daignan_fornier_2002_data.xlsx', 'data.txt');
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/euroscarf list.xlsx', '1_1.xlwb');

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, tested.raw(:,2)));
tested.raw(inds,:) = [];

inds = find(~strncmp('Y', tested.raw(:,2),1));
tested.raw(inds,:) = [];

tested.raw(:,2) = cellfun(@strtrim, tested.raw(:,2),'UniformOutput',0);

phenotypes = {'growth'};
treatments = {'MPA'};
phenotype_severity = {'S','SS','R','RR'};
phenotype_severity_num = [-1,-2,1,2];

% Eliminate white spaces before/after ORF
data.raw(1:end,1) = cellfun(@strtrim, data.raw(1:end,1),'UniformOutput',0);
inds = find(strcmp('YKR065W', data.raw(:,1)));
data.raw(inds,1) = {'YKR065C'};

% Replace string scores with numbers
t = data.raw(1:end,2);
for i = 1 : length(phenotype_severity)
inds = find(strcmp(phenotype_severity{i}, t));
t(inds) = {phenotype_severity_num(i)};
end
t(cellfun(@isnan,t))={0};

desmoucelles_daignan_fornier_2002.orfs = upper(tested.raw(:,2));
desmoucelles_daignan_fornier_2002.data = zeros(length(desmoucelles_daignan_fornier_2002.orfs),1);

[tmp,ind1,ind2] = intersect(desmoucelles_daignan_fornier_2002.orfs, upper(data.raw(1:end,1)));
desmoucelles_daignan_fornier_2002.data(ind1) = cell2mat(t(ind2));

desmoucelles_daignan_fornier_2002.ph = [strcat(phenotypes{1}, '; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(desmoucelles_daignan_fornier_2002.orfs, essential_genes);
desmoucelles_daignan_fornier_2002.orfs(ind1) = [];
desmoucelles_daignan_fornier_2002.data(ind1,:) = [];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'desmoucelles_daignan_fornier_2002.mat'],'desmoucelles_daignan_fornier_2002');
return;

% Save data into database
dt = desmoucelles_daignan_fornier_2002;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


end

