%% Uluisik~Koc, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
uluisik_koc_2011.pmid = 21035538;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(uluisik_koc_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested
[FILENAMES{end+1}, tested_orfs] = read_data('textread','./raw_data/tested_strains.txt', '%s');

% Eliminate all white spaces & capitalize
tested_orfs = clean_orf(tested_orfs);

% If in gene name form, transform into ORF name
tested_orfs = translate(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

% Get unique list of ORFs
tested_orfs = unique(tested_orfs);

%% Load resistant data 1
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/mmc2.xlsx');

% Get the list of ORFs and the correponding data 
hits_orfs_1 = data.raw(4:end,1);
hits_data_1 = data.raw(4:end,2:5);

% Translate sign into value
hits_data_1(strcmp('+',hits_data_1)) = {1};
hits_data_1(strcmp('-',hits_data_1)) = {0};

% Turn into numeric matrix
hits_data_1 = cell2mat(hits_data_1);

% Do numeric transformations of data
hits_data_1 = hits_data_1 - repmat(hits_data_1(1,:),size(hits_data_1,1),1);    % Subtract the 1st row (WT)
hits_data_1 = hits_data_1 - repmat(hits_data_1(:,1),1,size(hits_data_1,2));     % Subtract the 1st col (0 mM);
hits_data_1(1,:) = [];
hits_data_1(:,1) = [];
hits_orfs_1(1) = [];

% Eliminate all white spaces & capitalize
hits_orfs_1 = clean_orf(hits_orfs_1);

% If in gene name form, transform into ORF name
hits_orfs_1 = translate(hits_orfs_1);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs_1));
hits_orfs_1(inds) = [];
hits_data_1(inds,:) = [];

%% Load sensitive data 2
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/mmc3.xlsx');

% Get the list of ORFs and the correponding data 
hits_orfs_2 = data.raw(4:end,1);
hits_data_2 = data.raw(4:end,2:8);

% Translate sign into value
hits_data_2(strcmp('+',hits_data_2)) = {-1};
hits_data_2(strcmp('-',hits_data_2)) = {0};

% Turn into numeric matrix
hits_data_2 = cell2mat(hits_data_2);

% Eliminate empty ORFs
inds = find(cellfun(@isempty, hits_orfs_2) | cellfun(@isnumeric, hits_orfs_2));
hits_orfs_2(inds) = [];
hits_data_2(inds,:) = [];

% Do numeric transformations
hits_data_2 = hits_data_2 - repmat(hits_data_2(1,:),size(hits_data_2,1),1);    % Subtract the 1st row (WT)
hits_data_2 = hits_data_2 - repmat(hits_data_2(:,1),1,size(hits_data_2,2));     % Subtract the 1st col (0 mM);
hits_data_2(1,:) = [];
hits_data_2(:,1) = [];
hits_orfs_2(1) = [];

% Eliminate all white spaces & capitalize
hits_orfs_2 = clean_orf(hits_orfs_2);

% If in gene name form, transform into ORF name
hits_orfs_2 = translate(hits_orfs_2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs_2));
hits_orfs_2(inds) = [];
hits_data_2(inds,:) = [];

% Check for any hits missing from tested
[missing, ix] = setdiff(hits_orfs_2, tested_orfs);    % 3 ORFs missing (eliminated)
hits_orfs_2(ix) = [];
hits_data_2(ix,:) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [438; 439; 440; 441; 442; 443; 147; 444; 445];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
uluisik_koc_2011.orfs = tested_orfs;
uluisik_koc_2011.ph = hit_data_names;
uluisik_koc_2011.data = zeros(length(uluisik_koc_2011.orfs),length(uluisik_koc_2011.ph));
uluisik_koc_2011.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs_1, uluisik_koc_2011.orfs);
uluisik_koc_2011.data(ind2,7:9) = hits_data_1(ind1,:);

[~,ind1,ind2] = intersect(hits_orfs_2, uluisik_koc_2011.orfs);
uluisik_koc_2011.data(ind2,1:6) = hits_data_2(ind1,:);

%% Save
save('./uluisik_koc_2011.mat','uluisik_koc_2011');

%% Print out
fid = fopen('./uluisik_koc_2011.txt','w');
write_matrix_file(fid, uluisik_koc_2011.orfs, uluisik_koc_2011.ph, uluisik_koc_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(uluisik_koc_2011)
end

end
