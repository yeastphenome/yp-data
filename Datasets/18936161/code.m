%% Mir~Cashikar, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
mir_cashikar_2009.pmid = 18936161;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(mir_cashikar_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat_a.xlsx', 'mat_a_041902');
tested_orfs = tested.raw(3:end,2);

tested_orfs = clean_orf(tested_orfs);

tested_orfs(strcmp('YLR287-A', tested_orfs)) = {'YLR287C-A'};

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = translate(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));  

tested_orfs(inds) = [];
tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES, hits] = read_data('textscan', './raw_data/hits_orfs_scores.txt', '%s %d');
hits_orfs = hits{1};
hits_scores = -hits{2};

hits_orfs = clean_orf(hits_orfs);

hits_orfs = translate(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds)); 

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [100];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
mir_cashikar_2009.orfs = tested_orfs;
mir_cashikar_2009.ph = hit_data_names;
mir_cashikar_2009.data = zeros(length(mir_cashikar_2009.orfs),length(mir_cashikar_2009.ph));
mir_cashikar_2009.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, mir_cashikar_2009.orfs);
mir_cashikar_2009.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./mir_cashikar_2009.mat','mir_cashikar_2009');

%% Print out

fid = fopen('./mir_cashikar_2009.txt','w');
write_matrix_file(fid, mir_cashikar_2009.orfs, mir_cashikar_2009.ph, mir_cashikar_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(mir_cashikar_2009)
end

end
