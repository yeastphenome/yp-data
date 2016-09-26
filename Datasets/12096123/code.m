%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
wilson_roach_2002.pmid = 12096123;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(wilson_roach_2002.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/suppdata.xlsx', 'SuppDataREV');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,3);

% Get the data itself
hit_data = data(:,2);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(strcmp('YORO36W', hit_strains)) = {'YOR036W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds) = [];

hit_data = strtrim(hit_data);

[dt,~,ic] = unique(hit_data);
dt_n = [1 1 2 -1 -2]';
hit_data = dt_n(ic);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [4949];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('xlsread','./raw_data/ResGen Diploids inventory.xlsx', 'Inventory');

tested_strains = tested_strains(:,1);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'TAL004W'})) = {'YAL004W'};
tested_strains(ismember(tested_strains, {'YELOO1C'})) = {'YEL001C'};
tested_strains(ismember(tested_strains, {'KL187C'})) = {'YKL187C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

tested_strains(inds) = [];

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
wilson_roach_2002.orfs = tested_strains;
wilson_roach_2002.ph = hit_data_names;
wilson_roach_2002.data = zeros(length(wilson_roach_2002.orfs),length(wilson_roach_2002.ph));
wilson_roach_2002.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, wilson_roach_2002.orfs);
wilson_roach_2002.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./wilson_roach_2002.mat','wilson_roach_2002');

%% Print out

fid = fopen('./wilson_roach_2002.txt','w');
write_matrix_file(fid, wilson_roach_2002.orfs, wilson_roach_2002.ph, wilson_roach_2002.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(wilson_roach_2002)
end

end

