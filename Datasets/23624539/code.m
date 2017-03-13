%% Galvan Marquez~Smith, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
galvan_marquez_smith_2013.pmid = 23624539;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(galvan_marquez_smith_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Chitosan Effect on GDA (raw data). Exp. 1 to 3. Imelda Galvan, 2013-1.xlsx', 'Sheet1');

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

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [129];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
galvan_marquez_smith_2013.orfs = t;
galvan_marquez_smith_2013.ph = hit_data_names;
galvan_marquez_smith_2013.data = t2;
galvan_marquez_smith_2013.dataset_ids = hit_data_ids;

%% Save

save('./galvan_marquez_smith_2013.mat','galvan_marquez_smith_2013');

%% Print out

fid = fopen('./galvan_marquez_smith_2013.txt','w');
write_matrix_file(fid, galvan_marquez_smith_2013.orfs, galvan_marquez_smith_2013.ph, galvan_marquez_smith_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(galvan_marquez_smith_2013)
end

end
