%% Chesi~Gitler, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
chesi_gitler_2012.pmid = 22457822;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(chesi_gitler_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/yeast deletions Mn.xlsx', 'single deletion');

crr = zeros(size(data.raw,2)-1,3);  % Concentration, within-round replicate, round
for i = 2 : length(data.raw(1,:))
t = textscan(data.raw{1,i},'%d %s size-%d %d');
crr(i-1,:) = cell2mat(t([1 3 4]));
end

% Get indices of the data columns
ind_data = 2:19;

% Eliminate all white spaces & capitalize
data.raw(:,1) = clean_orf(data.raw(:,1));

% If in gene name form, transform into ORF name
data.raw(:,1) = translate(data.raw(:,1));

% Find anything that doesn't look like an ORF
temp = data.raw(:,1);
temp(ismember(data.raw(:,1), {'YLR287-A'})) = {'YLR287C-A'};
data.raw(:,1) = temp;
inds = find(~is_orf(data.raw(:,1)));
data.raw(inds, :) = [];

% Make sure all the data are numbers
t = data.raw(:,ind_data);
inds = find(~cellfun(@isnumeric, t));
t(inds) = {NaN};

ind_data_conc1 = [1 2 7 8 13 14];
ind_data_conc2 = [3 4 9 10 15 16];
ind_data_conc3 = [5 6 11 12 17 18];

data2.orfs = upper(data.raw(:,1));
data2.data(:,1) = nanmean(cell2mat(t(:,ind_data_conc1)),2);
data2.data(:,2) = nanmean(cell2mat(t(:,ind_data_conc2)),2);
data2.data(:,3) = nanmean(cell2mat(t(:,ind_data_conc3)),2);

% Normalize to no drug
data2.data(:,2) = data2.data(:,2) ./ data2.data(:,1);
data2.data(:,3) = data2.data(:,3) ./ data2.data(:,1);
data2.data(:,1) = [];

concentrations = unique(crr(:,1));
concentrations(1) = [];

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [247; 248];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
chesi_gitler_2012.orfs = t;
chesi_gitler_2012.ph = hit_data_names;
chesi_gitler_2012.data = t2;
chesi_gitler_2012.dataset_ids = hit_data_ids;

%% Save

save('./chesi_gitler_2012.mat','chesi_gitler_2012');

%% Print out

fid = fopen('./chesi_gitler_2012.txt','w');
write_matrix_file(fid, chesi_gitler_2012.orfs, chesi_gitler_2012.ph, chesi_gitler_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(chesi_gitler_2012)
end

end

