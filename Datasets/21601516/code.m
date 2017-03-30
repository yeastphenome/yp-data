%% Suzuki~Yoshida, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
suzuki_yoshida_2011.pmid = 21601516;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(suzuki_yoshida_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/hits_data.xlsx');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hits_orfs = data.raw(1:end,1);

% Get the data itself
hits_data = cell2mat(data.raw(1:end,2));

% Eliminate all white spaces & capitalize
hits_orfs = clean_genename(hits_orfs);

% If in gene name form, transform into ORF name
hits_orfs = translate(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

% If the same strain is present more than once, average its values
[hits_orfs,hits_data] = grpstats(hits_data, hits_orfs,{'gname','mean'});

%% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/YKOmatalpha_GSH_list070508.xlsx');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
tested_orfs = tested.raw(3:end,3);

% Eliminate all white spaces & capitalize
tested_orfs = clean_orf(tested_orfs);

% Find empty matrices and remove them
inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

% If in gene name form, transform into ORF name
tested_orfs = translate(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

% Get a unique list of ORFs
tested_orfs = unique(tested_orfs);

% Isolate difference between tested and orfs
[missing, ix] = setdiff(hits_orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [125];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
suzuki_yoshida_2011.orfs = tested_orfs;
suzuki_yoshida_2011.ph = hit_data_names;
suzuki_yoshida_2011.data = zeros(length(suzuki_yoshida_2011.orfs),length(suzuki_yoshida_2011.ph));
suzuki_yoshida_2011.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, suzuki_yoshida_2011.orfs);
suzuki_yoshida_2011.data(ind2,:) = hits_data(ind1,:);

%% Save

save('./suzuki_yoshida_2011.mat','suzuki_yoshida_2011');

%% Print out

fid = fopen('./suzuki_yoshida_2011.txt','w');
write_matrix_file(fid, suzuki_yoshida_2011.orfs, suzuki_yoshida_2011.ph, suzuki_yoshida_2011.data);
fclose(fid);


%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(suzuki_yoshida_2011)
end

end
