%% Xie~Huang, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};


xie_huang_2005.source = {'Jing Huang'};
xie_huang_2005.downloaddate = {'2013-03-21'};
xie_huang_2005.pmid = 15883373;

phenotypes = {'Growth, spot intensity on cell microarray'};
treatments = {'Rapamycin, 10 nM'; 'Rapamycin, 30 nM'};

[FILENAMES{end+1}, data.raw] = readdata('xlsread','./raw_data/rapa_cellarray_all.xls', 'Rapa10nm');
[FILENAMES{end+1}, data2.raw] = readdata('xlsread','./raw_data/rapa_cellarray_all.xls', 'Rapa30nm');

% Get indices of the data columns
ind_data = find(strcmp('Rapa/DMSO', data.raw(1,:)));
ind_data2 = find(strcmp('Rapa/DMSO', data2.raw(1,:)));

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,3)));
data.raw(inds,:) = [];
inds = find(cellfun(@isnumeric, data2.raw(:,3)));
data2.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,3)),strmatch('Y', data.raw(:,3)));
data.raw(inds,:) = [];
inds = setdiff(1:length(data2.raw(:,3)), strmatch('Y', data2.raw(:,3)));
data2.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,3) = cellfun(@strtrim, data.raw(:,3),'UniformOutput',0);
data2.raw(:,3) = cellfun(@strtrim, data2.raw(:,3),'UniformOutput',0);

data3.orfs = upper(data.raw(:,3));
data3.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data3.data)==0);
data3.data(inds) = {NaN};
data3.data = cell2mat(data3.data);

data4.orfs = upper(data2.raw(:,3));
data4.data = data2.raw(:,ind_data2);
inds = find(cellfun(@isnumeric, data4.data)==0);
data4.data(inds) = {NaN};
data4.data = cell2mat(data4.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data3.data, data3.orfs, {'gname','mean'});
xie_huang_2005.orfs = t;
xie_huang_2005.data = t2;
xie_huang_2005.ph = strcat(phenotypes, '; ', treatments);

[t,t2] = grpstats(data4.data, data4.orfs, {'gname','mean'});

[t3,ind1,ind2] = intersect(xie_huang_2005.orfs, t);
xie_huang_2005.data(ind1,end+size(t2,2)) = t2(ind2,:);

% The list contains both non-essential deletions and essential heterozygous
% mutants. Eliminate the essential genes

load essential_genes_100908;
[t,ind1,ind2] = intersect(xie_huang_2005.orfs, essential_genes);
xie_huang_2005.orfs(ind1) = [];
xie_huang_2005.data(ind1,:) = [];

save('./xie_huang_2005.mat','xie_huang_2005');
return;

% Save data into database
dt = xie_huang_2005;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
datasets_ids(i,1) = datasets(i).id;
datasets_names{i,1} = datasets(i).name;
if isempty(datasets(i).reporter)
datasets_names{i,2} = '';
else
datasets_names{i,2} = datasets(i).reporter;
end
datasets_names{i,3} = datasets(i).conditionset;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));


end

