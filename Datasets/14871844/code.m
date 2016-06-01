%% Aouida~Ramotar, 2004
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
aouida_ramotar_2004.pmid = 14871844;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(aouida_ramotar_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/HU haploid.xlsx');
tested_orfs = tested.raw(4:end,5);

% Eliminate all white spaces & capitalize
tested_orfs = clean_orf(tested_orfs);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));  

tested_orfs = unique(tested_orfs);


%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/CAN_2-1-04_Aouida.xlsx');
hits_orfs = data.raw(:,1);

% Eliminate all white spaces & capitalize
hits_orfs = clean_orf(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));

hits_data = strtrim(data.raw(:,2));

% Data conversion
hits_data(strcmp('R+++', hits_data)) = {9};     % (>500-fold more resistant -> log2(500) ~ 9)
hits_data(strcmp('R++++', hits_data)) = {10};
hits_data(strcmp('S+', hits_data)) = {-1};
hits_data(strcmp('S++', hits_data)) = {-2};
hits_data(strcmp('S+++', hits_data)) = {-3};
hits_data(strcmp('S++++', hits_data)) = {-4};

hits_data = cell2mat(hits_data);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

%%

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [96];

%%

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
aouida_ramotar_2004.orfs = tested_orfs;
aouida_ramotar_2004.ph = hit_data_names;
aouida_ramotar_2004.data = zeros(length(aouida_ramotar_2004.orfs),length(aouida_ramotar_2004.ph));
aouida_ramotar_2004.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, aouida_ramotar_2004.orfs);
aouida_ramotar_2004.data(ind2,:) = hits_data(ind1,:);


save('./aouida_ramotar_2004.mat','aouida_ramotar_2004');

fid = fopen('./aouida_ramotar_2004.txt','w');
write_matrix_file(fid, aouida_ramotar_2004.orfs, aouida_ramotar_2004.ph, aouida_ramotar_2004.data);
fclose(fid);

end

