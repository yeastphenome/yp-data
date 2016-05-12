%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
firstauthor_lastauthor_YYYY.pmid = 12345678;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', './extras/YeastPhenome_<PMID>_datasets_list.txt','%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/data.xlsx', 'Spreadsheet name');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);

% Get the data itself
hit_data = ones(size(hit_strains)); % if the dataset is binary
hit_data = data(:,2:4); % if the dataset is discrete or binary
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YAL001'})) = {'YAL001C'};

% If not possible, eliminate the entry
hit_strains(ismember(hit_strains, {'BLANK'})) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [1 2 3 4 5 6];
[hit_data_ids, hit_data] = grpstats(hit_data', hit_data_ids, {'gname','mean'});
hit_data = hit_data';
hit_data_ids = cellfun(@str2num, hit_data_ids);

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('xlsread','./raw_data/tested_strains.xlsx', 'Spreadsheet name');

% Eliminate all white spaces & capitalize
tested_strains = clean_genename(tested_strains);

% If in gene name form, transform into ORF name
tested_strains = translate(tested_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'YAL001'})) = {'YAL001C'};

% If not possible, eliminate the entry
tested_strains(ismember(tested_strains, {'BLANK'})) = [];

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

% If it seems reasonable, add the missing hits to the list of tested strains
tested_strains = [tested_strains; missing];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
firstauthor_lastauthor_YYYY.orfs = hit_strains;
firstauthor_lastauthor_YYYY.ph = hit_data_names;
firstauthor_lastauthor_YYYY.data = hit_data;
firstauthor_lastauthor_YYYY.dataset_ids = hit_data_ids;

% If the dataset is discrete/binary and the tested strains were provided separately:
firstauthor_lastauthor_YYYY.orfs = tested_strains;
firstauthor_lastauthor_YYYY.ph = hit_data_names;
firstauthor_lastauthor_YYYY.data = zeros(length(firstauthor_lastauthor_YYYY.orfs),length(firstauthor_lastauthor_YYYY.ph));
firstauthor_lastauthor_YYYY.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, firstauthor_lastauthor_YYYY.orfs);
firstauthor_lastauthor.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./firstauthor_lastauthor_YYYY.mat','firstauthor_lastauthor_YYYY');

%% Print out

fid = fopen('./firstauthor_lastauthor_YYYY.txt','w');
write_matrix_file(fid, firstauthor_lastauthor_YYYY.orfs, firstauthor_lastauthor_YYYY.ph, firstauthor_lastauthor_YYYY.data);
fclose(fid);

end

