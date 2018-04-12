%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
choy_basrai_2013.pmid = 23825022;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(choy_basrai_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/Final HIS Data.xlsx', 'X YMB1004');
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/Final HIS Data.xlsx', 'X YMB1005');

% Get the list of ORFs and the correponding data 
hit_strains1 = data1(2:end,3);
hit_strains2 = data2(2:end,3);

% Get the data itself
hit_data1 = data1(2:end,4:6);
hit_data2 = data2(2:end,4:6);

% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);
hit_strains2 = clean_orf(hit_strains2);

inds = find(cellfun(@isnumeric, hit_strains1) | strcmp('NA', hit_strains1));
hit_strains1(inds) = [];
hit_data1(inds,:) = [];

inds = find(cellfun(@isnumeric, hit_strains2) | strcmp('NA', hit_strains2));
hit_strains2(inds) = [];
hit_data2(inds,:) = [];

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);

% If possible, fix the problem (typos, omissions etc.)
hit_strains1(ismember(hit_strains1, {'YCLO51W'})) = {'YCL051W'};
hit_strains1(ismember(hit_strains1, {'YGR122C-'})) = {'YGR122C-A'};
hit_strains1(ismember(hit_strains1, {'YHR139C-'})) = {'YHR139C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

t = hit_strains1(inds);
t = cellfun(@(S) S(1:end-1), t, 'Uniform', 0);

hit_strains1(inds) = t;


% If in gene name form, transform into ORF name
hit_strains2 = translate(hit_strains2);

% If possible, fix the problem (typos, omissions etc.)
hit_strains2(ismember(hit_strains2, {'YCLO51W'})) = {'YCL051W'};
hit_strains2(ismember(hit_strains2, {'YGR122C-'})) = {'YGR122C-A'};
hit_strains2(ismember(hit_strains2, {'YHR139C-'})) = {'YHR139C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

t = hit_strains2(inds);
t = cellfun(@(S) S(1:end-1), t, 'Uniform', 0);

hit_strains2(inds) = t;

hit_data1 = cell2mat(hit_data1);
hit_data2 = cell2mat(hit_data2);


% Sum all the values and take the average once
hit_data1 = nanmean(hit_data1, 2);
hit_data2 = nanmean(hit_data2, 2);

hit_strains = [hit_strains1; hit_strains2];
hit_data = [hit_data1; hit_data2];

[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16003];

%% TEST
[FILENAMES{end+1}, data3] = read_data('xlsread','./raw_data/table_s1.xlsx', 'Table 1');

hit_strains3 = data3(4:end,1);
hit_data3 = data3(4:end,2);

inds = find(cellfun(@isnumeric, hit_strains3));
hit_strains3(inds) = [];
hit_data3(inds) = [];

hit_data3  = cell2mat(hit_data3);

[~,ind1,ind2] = intersect(hit_strains3, hit_strains);
plot(hit_data3(ind1), hit_data(ind2),'o')

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
choy_basrai_2013.orfs = hit_strains;
choy_basrai_2013.ph = hit_data_names;
choy_basrai_2013.data = hit_data;
choy_basrai_2013.dataset_ids = hit_data_ids;

%% Save

save('./choy_basrai_2013.mat','choy_basrai_2013');

%% Print out

fid = fopen('./choy_basrai_2013.txt','w');
write_matrix_file(fid, choy_basrai_2013.orfs, choy_basrai_2013.ph, choy_basrai_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(choy_basrai_2013)
end

end

