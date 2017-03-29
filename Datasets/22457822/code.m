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

hit_strains = data.raw(:,1);
hit_data = data.raw(:,ind_data);

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

hit_strains(ismember(hit_strains, {'YLR287-A'})) = {'YLR287C-A'};

inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% Make sure all the data are numbers
inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {NaN};

hit_data = cell2mat(hit_data);

ind_data_conc1 = [1 2 7 8 13 14];
ind_data_conc2 = [3 4 9 10 15 16];
ind_data_conc3 = [5 6 11 12 17 18];

hit_data2(:,1) = nanmean(hit_data(:,ind_data_conc1),2);
hit_data2(:,2) = nanmean(hit_data(:,ind_data_conc2),2);
hit_data2(:,3) = nanmean(hit_data(:,ind_data_conc3),2);

% Normalize to no drug
hit_data2(:,2) = hit_data2(:,2) ./ hit_data2(:,1);
hit_data2(:,3) = hit_data2(:,3) ./ hit_data2(:,1);
hit_data2(:,1) = [];

concentrations = unique(crr(:,1));
concentrations(1) = [];

% Average data for identical ORFs that appear multiple times
[hit_strains, hit_data] = grpstats(hit_data2, hit_strains, {'gname','mean'});

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
chesi_gitler_2012.orfs = hit_strains;
chesi_gitler_2012.ph = hit_data_names;
chesi_gitler_2012.data = hit_data;
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

