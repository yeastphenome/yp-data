%% Svensson~Samson, 2011
function FILENAMES = code()
FILENAMES = {};


svensson_samson_2011.source = {'http://www.biomedcentral.com/content/supplementary/1752-0509-5-157-s1.xls'};
svensson_samson_2011.downloaddate = {'2013-03-08'};
svensson_samson_2011.pmid = 21978764;

phenotypes = {'Growth, GI50'};
treatments = {'MMS'};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/1752-0509-5-157-s1.xlsx', '2. Gi50 and R2 all strains');

% Get indices of the data columns
ind_data = 4;

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,2)),strmatch('Y', data.raw(:,2)));
data.raw(inds,:) = [];

% Separate deletions from DAMP strains. Personal communication from
% Peter Svensson: deletions are on plates 1-57, DAMPs are on plates
% 301-311

dels = 'ABCDEFGH';
plate = data.raw(:,1);
for i = 1 : length(dels)
plate = strtok(plate, dels(i));
end
plate = cellfun(@str2num, plate);

inds = find(plate > 57);
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
svensson_samson_2011.orfs = t;
svensson_samson_2011.data = t2;
svensson_samson_2011.ph = strcat(phenotypes, {'; '} , treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'svensson_samson_2011.mat'],'svensson_samson_2011');
return;

% Save data into database
dt = svensson_samson_2011;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
datasets_ids(i,1) = datasets(i).id;
datasets_names{i,1} = datasets(i).name;
datasets_names{i,2} = datasets(i).reporter;
datasets_names{i,3} = datasets(i).short_name;
datasets_names{i,4} = datasets(i).dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3 4]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));


end

