%% Desmoucelles~Daignan-Fornier, 2002
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
desmoucelles_daignan_fornier_2002.pmid = 12016207;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(desmoucelles_daignan_fornier_2002.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested strains

% NOTE = 1) typo in data: YKR065W should be YKR065C 2) YCR002C is in the
% result, but not in tested.

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/euroscarf list.xlsx', '1_1.xlwb');

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested.raw(:,2));

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));

tested_strains(inds) = [];

tested_strains = unique(tested_strains);

%% Load the data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/desmoucelles_daignan_fornier_2002_data.xlsx', 'data.txt');

phenotype_severity = {'S','SS','R','RR'};
phenotype_severity_num = [-1,-2,1,2];

hit_strains = data.raw(:,1);
hit_data = data.raw(:,2);

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YKR065W'})) = {'YKR065C'};

% Replace string scores with numbers
for i = 1 : length(phenotype_severity)
    inds = find(strcmp(phenotype_severity{i}, hit_data));
    hit_data(inds) = {phenotype_severity_num(i)};
end
hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [66];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
desmoucelles_daignan_fornier_2002.orfs = tested_strains;
desmoucelles_daignan_fornier_2002.ph = hit_data_names;
desmoucelles_daignan_fornier_2002.data = zeros(length(desmoucelles_daignan_fornier_2002.orfs),length(desmoucelles_daignan_fornier_2002.ph));
desmoucelles_daignan_fornier_2002.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, desmoucelles_daignan_fornier_2002.orfs);
desmoucelles_daignan_fornier_2002.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./desmoucelles_daignan_fornier_2002.mat','desmoucelles_daignan_fornier_2002');

%% Print out

fid = fopen('./desmoucelles_daignan_fornier_2002.txt','w');
write_matrix_file(fid, desmoucelles_daignan_fornier_2002.orfs, desmoucelles_daignan_fornier_2002.ph, desmoucelles_daignan_fornier_2002.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(desmoucelles_daignan_fornier_2002)
end

end

