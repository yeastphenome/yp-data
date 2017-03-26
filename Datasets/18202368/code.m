%% Nyswaner~Garfinkel,2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
nyswaner_garfinkel_2008.pmid = 18202368;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(nyswaner_garfinkel_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Matalphakos counted.xlsx', 'Sheet2');
tested_orfs = tested.raw(6:end,2);

tested_orfs = clean_orf(tested_orfs);

% Fix typos
tested_orfs(find(strcmp('YLR228', tested_orfs))) = {'YLR228C'};
tested_orfs(find(strcmp('YMR062', tested_orfs))) = {'YMR062C'};
tested_orfs(find(strcmp('YYKL138C', tested_orfs))) = {'YKL138C'};

inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));

tested_orfs(inds) = [];
tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/genetics.107.082602-9.xlsx', 'Sheet1');
hits_orfs = data.raw(:,1);

hits_orfs = clean_orf(hits_orfs);

inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));

hits_orfs(inds) = [];

missing = setdiff(hits_orfs, tested_orfs);

% Adding 3 ORFs to the list of tested
tested_orfs = [tested_orfs; missing];

hits_data = ones(size(hits_orfs));

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [164];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
nyswaner_garfinkel_2008.orfs = tested_orfs;
nyswaner_garfinkel_2008.ph = hit_data_names;
nyswaner_garfinkel_2008.data = zeros(length(nyswaner_garfinkel_2008.orfs),length(nyswaner_garfinkel_2008.ph));
nyswaner_garfinkel_2008.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, nyswaner_garfinkel_2008.orfs);
nyswaner_garfinkel_2008.data(ind2,:) = hits_data(ind1,:);

%% Save

save('./nyswaner_garfinkel_2008.mat','nyswaner_garfinkel_2008');

%% Print out

fid = fopen('./nyswaner_garfinkel_2008.txt','w');
write_matrix_file(fid, nyswaner_garfinkel_2008.orfs, nyswaner_garfinkel_2008.ph, nyswaner_garfinkel_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(nyswaner_garfinkel_2008)
end

end
