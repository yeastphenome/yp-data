%% Stauffer~Powers, 2015
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
stauffer_powers_2015.pmid = 26466677;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(stauffer_powers_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the hit strains

[FILENAMES{end+1}, hit_strains] = read_data('textread','./raw_data/Table1.txt', '%s', 'delimiter', '\n');
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% Translate the hit strains
hit_strains = translate(hit_strains);

hit_strains(ismember(hit_strains, {'RPL19-A'})) = {'YBR084C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];

hit_strains = clean_orf(hit_strains);
hit_strains = unique(hit_strains);

% Get the data itself
hit_data = zeros(size(hit_strains))-1;

%% Load tested strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Deletion Collection List for Anastasia.xlsx');

% Collect tested strains
tested_strains = data(5:end, 1);

% Clean the orfs
tested_strains = clean_orf(tested_strains);

% Find anything that doesn't look like an ORF
tested_strains(cellfun(@isnumeric, tested_strains)) = [];
tested_strains(ismember(tested_strains, {'YLR287-A'})) = {'YLR287C-A'};
inds = find(~is_orf(tested_strains));
tested_strains(inds) = [];

tested_strains = translate(tested_strains);

tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [766];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
stauffer_powers_2015.orfs = tested_strains;
stauffer_powers_2015.ph = hit_data_names;
stauffer_powers_2015.data = zeros(length(stauffer_powers_2015.orfs),length(stauffer_powers_2015.ph));
stauffer_powers_2015.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, stauffer_powers_2015.orfs);
stauffer_powers_2015.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./stauffer_powers_2015.mat','stauffer_powers_2015');

%% Print out

fid = fopen('./stauffer_powers_2015.txt','w');
write_matrix_file(fid, stauffer_powers_2015.orfs, stauffer_powers_2015.ph, stauffer_powers_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(stauffer_powers_2015)
end

end

