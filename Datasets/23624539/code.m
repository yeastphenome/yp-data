%% Galvan Marquez~Smith, 2013
% DATA = galvan_marquez_smith_2013
function FILENAMES = code()
FILENAMES = {};


galvan_marquez_smith_2013.source = {'Imelda Galvan Marquez'};
galvan_marquez_smith_2013.downloaddate = {'2013-05-19'};
galvan_marquez_smith_2013.pmid = 23624539;

phenotypes = {'Growth, colony size'};
treatments = {'Chitosan'};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/Chitosan Effect on GDA (raw data). Exp. 1 to 3. Imelda Galvan, 2013-1.xlsx', 'Sheet1');

% Get indices of the data columns
ind_data = find(strcmp('Ratio', data.raw(1,:)));

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,2)),strmatch('Y', data.raw(:,2)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,3) = cellfun(@strtrim, data.raw(:,2),'UniformOutput',0);

% Make sure all the data are numbers
t = data.raw(:,ind_data);
inds = find(cellfun(@isnumeric, t)==0);
t(inds) = {NaN};

data2.orfs = upper(data.raw(:,2));
data2.data = nanmean(cell2mat(t),2);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
galvan_marquez_smith_2013.orfs = t;
galvan_marquez_smith_2013.data = t2;
galvan_marquez_smith_2013.ph = strcat(phenotypes, '; ', treatments);


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(galvan_marquez_smith_2013.orfs, essential_genes);
galvan_marquez_smith_2013.orfs(ind1) = [];
galvan_marquez_smith_2013.data(ind1,:) = [];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'galvan_marquez_smith_2013.mat'],'galvan_marquez_smith_2013');
return;

% Save data into database
dt = galvan_marquez_smith_2013;

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

