%% Luban~Schmidt, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
luban_schmidt_2005.pmid = 15908144;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(luban_schmidt_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested_orfs] = read_data('textscan', './raw_data/list_of_used_knockouts_PhD_Thesis_Luban.txt', '%s\n');
tested_orfs = clean_orf(tested_orfs);

tested_orfs(strcmp('YBL098V', tested_orfs)) = {'YBL098W'};
tested_orfs(strcmp('YDR07SW', tested_orfs)) = {'YDR075W'};
tested_orfs(strcmp('YDR27SW', tested_orfs)) = {'YDR275W'};
tested_orfs(strcmp('YDR51SW', tested_orfs)) = {'YDR515W'};
tested_orfs(strcmp('YDRS41C', tested_orfs)) = {'YDR541C'};
tested_orfs(strcmp('YEL0I6C', tested_orfs)) = {'YEL016C'};
tested_orfs(strcmp('YIIL016C', tested_orfs)) = {'YHL016C'};
tested_orfs(strcmp('YIIL017W', tested_orfs)) = {'YHL017W'};
tested_orfs(strcmp('YHL0I9C', tested_orfs)) = {'YHL019C'};
tested_orfs(strcmp('YJR09JC', tested_orfs)) = {'YJR091C'};
tested_orfs(strcmp('YNL09SC', tested_orfs)) = {'YNL095C'};
tested_orfs(strcmp('YPLOI8W', tested_orfs)) = {'YPL018W'};
tested_orfs(strcmp('YPL07LC', tested_orfs)) = {'YPL071C'};

inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));

tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, hit_orfs] = read_data('textscan', './raw_data/list_of_pet_mutants.txt', '%s');
hit_orfs = clean_orf(hit_orfs);

inds = find(~is_orf(hit_orfs));
disp(hit_orfs(inds));

hit_orfs = unique(hit_orfs);
hit_data = -ones(size(hit_orfs));

missing = setdiff(hit_orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [417];

%% Prepare the final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
luban_schmidt_2005.orfs = tested_orfs;
luban_schmidt_2005.ph = hit_data_names;
luban_schmidt_2005.data = zeros(length(luban_schmidt_2005.orfs),length(luban_schmidt_2005.ph));
luban_schmidt_2005.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_orfs, luban_schmidt_2005.orfs);
luban_schmidt_2005.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./luban_schmidt_2005.mat','luban_schmidt_2005');

%% Print out

fid = fopen('./luban_schmidt_2005.txt','w');
write_matrix_file(fid, luban_schmidt_2005.orfs, luban_schmidt_2005.ph, luban_schmidt_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(luban_schmidt_2005)
end

end
