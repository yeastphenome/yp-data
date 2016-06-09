%% Vasquez-Soto~Norambuena, 2015
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
vasquez_soto_norambuena_2015.pmid = 26209329;

hit_data_ids = [767];

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(vasquez_soto_norambuena_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% All Strains
[FILENAMES{end+1}, all_orfs] = read_data('xlsread','./raw_data/YSC1054.YKO.mat_alpha.v1.0.xlsx', 'mat_alpha_obs');

% Get the list of orfs
tested_orfs = all_orfs(:,2);

% Clean them up
tested_orfs = clean_orf(tested_orfs);

% Find anything that doesn't look like an ORF and remove it
inds = find(~is_orf(tested_orfs));
tested_orfs(inds, :) = [];

%% Hit Strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/40659_2015_32_MOESM1_ESM.xlsx', 'Sheet1');

% Get the list of ORFs
hit_orfs = data(3:end, 1);
hit_orfs = clean_orf(hit_orfs);

% Find anything that doesn't look like an ORF and remove it
inds = cellfun(@isnumeric, hit_orfs);
hit_orfs(inds, :) = [];
data(inds, :) = [];

%% Data
% Make an array of zeros
final_data = zeros(size(tested_orfs));

% Separate sensitive vs resistant in hit_data
hit_data = data(3:end, 4);
hit_data(strcmp('R', hit_data)) = {-2};
indx = ~cellfun(@isnumeric, hit_data);
hit_data(indx) = {-1};
hit_data = cell2mat(hit_data);

% Match the hit names with all the names
[~,ind1,ind2] = intersect(tested_orfs, hit_orfs);
final_data(ind1) = hit_data(ind2);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

vasquez_soto_norambuena_2015.orfs = tested_orfs;
vasquez_soto_norambuena_2015.ph = hit_data_names;
vasquez_soto_norambuena_2015.data = final_data;
vasquez_soto_norambuena_2015.dataset_ids = hit_data_ids;

%% Save

save('./vasquez_soto_norambuena_2015.mat','vasquez_soto_norambuena_2015');

fid = fopen('./vasquez_soto_norambuena_2015.txt','w');
write_matrix_file(fid, vasquez_soto_norambuena_2015.orfs, vasquez_soto_norambuena_2015.ph, vasquez_soto_norambuena_2015.data);
fclose(fid);

end


