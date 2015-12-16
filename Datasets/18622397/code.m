%% Breslow~Weissman, 2008
function FILENAMES = code()
FILENAMES = {};

breslow_weissman_2008.source = {'http://www.nature.com/nmeth/journal/v5/n8/extref/nmeth.1234-S4.xls'};
breslow_weissman_2008.downloaddate = {'2013-03-12'};
breslow_weissman_2008.pmid = 18622397;

phenotypes = {'Growth, flow cytometry'};
treatments = {''};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/nmeth.1234-S5.xlsx', 'Deletion growth rates');

% Get indices of the data columns
ind_data = find(strcmp('Median Growth Rate', data.raw(1,:)));

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,2)),strmatch('Y', data.raw(:,2)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,2) = cellfun(@strtrim, data.raw(:,2),'UniformOutput',0);

data2.orfs = upper(data.raw(:,2));
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
breslow_weissman_2008.orfs = t;
breslow_weissman_2008.data = t2;
breslow_weissman_2008.ph = phenotypes;

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'breslow_weissman_2008.mat'],'breslow_weissman_2008');
return;

% Save data into database
dt = breslow_weissman_2008;

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

