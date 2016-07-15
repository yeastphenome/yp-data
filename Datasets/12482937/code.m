%% Chang~Brown, 2002
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
chang_brown_2002.pmid = 12482937;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(chang_brown_2002.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the hit_data

[FILENAMES{end+1}, data] = read_data('textscan','./raw_data/MMS Sensitive Strains.txt', '%s', 'Delimiter', '\n');

% Retrieve the strains and the data
hit_strains = {};
hit_data = [];
for i = 1:length(data)
    C = strsplit(data{i}, ' ');
    tmp = strsplit(C{1}, '/');
    hit_strains(i,1) = tmp(1);
    hit_data(i,1) = -length(C{end});
end

hit_strains = clean_genename(hit_strains);
hit_strains(strcmp('KIM3', hit_strains)) = {'YPR164W'}; % old name, no longer used in SGD, fixed manually

[hit_strains, translated] = translate(hit_strains);

hit_strains(~translated) = [];
hit_data(~translated) = [];

%% Tested strains

[FILENAMES{end+1}, tested_strains] = read_data('xlsread','./raw_data/Old array position-before Dec 8, 05.xlsx', 'Sheet1');
tested_strains = tested_strains(2:end, 4);

tested_strains = clean_orf(tested_strains);

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'YPL072WA'})) = {'YPL072W-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
tested_strains(inds) = [];

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

% 1 missing hit, adding to tested
tested_strains = [tested_strains; missing];

% Intersect the tested with the hits
all_data = zeros(length(tested_strains), 1);
[~, ind1, ind2] = intersect(tested_strains, hit_strains);
all_data(ind1) = hit_data(ind2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [65];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
chang_brown_2002.orfs = tested_strains;
chang_brown_2002.ph = hit_data_names;
chang_brown_2002.data = all_data;
chang_brown_2002.dataset_ids = hit_data_ids;

%% Save

save('./chang_brown_2002.mat','chang_brown_2002');

%% Print out

fid = fopen('./chang_brown_2002.txt','w');
write_matrix_file(fid, chang_brown_2002.orfs, chang_brown_2002.ph, chang_brown_2002.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(chang_brown_2002)
end

end

