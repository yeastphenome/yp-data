%% Chamilos~Kontoyiannis, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
chamilos_kontoyiannis_2008.pmid = 18212113;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(chamilos_kontoyiannis_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested_orfs] = read_data('textread','./raw_data/tested_orfs.txt', '%s');

tested_orfs = clean_orf(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds)); 

tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, C] = read_data('textscan', './raw_data/data_genenames.txt', '%s %f','delimiter','\t');

genenames = C{1};
raw_data = C{2};

% Normalize to WT
inds = find(strcmp('WT', genenames));
raw_data = raw_data./raw_data(inds)-1;

genenames(inds) = [];
raw_data(inds) = [];

orfs = translate(genenames);

missing_orfs = setdiff(orfs, tested_orfs);
tested_orfs = [tested_orfs; missing_orfs];  % Adding 2 missing strains to the list of tested strains.

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [80];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
chamilos_kontoyiannis_2008.orfs = tested_orfs;
chamilos_kontoyiannis_2008.ph = hit_data_names;
chamilos_kontoyiannis_2008.data = zeros(length(chamilos_kontoyiannis_2008.orfs),length(chamilos_kontoyiannis_2008.ph));
chamilos_kontoyiannis_2008.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(orfs, chamilos_kontoyiannis_2008.orfs);
chamilos_kontoyiannis_2008.data(ind2,:) = raw_data(ind1,:);

%% Save

save('./chamilos_kontoyiannis_2008.mat','chamilos_kontoyiannis_2008');

%% Print out

fid = fopen('./chamilos_kontoyiannis_2008.txt','w');
write_matrix_file(fid, chamilos_kontoyiannis_2008.orfs, chamilos_kontoyiannis_2008.ph, chamilos_kontoyiannis_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(chamilos_kontoyiannis_2008)
end

end
