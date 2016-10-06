%% Prescott~Simmonds, 2014

function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
prescott_simmonds_2014.pmid = 24569176;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(prescott_simmonds_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the hit strains
% Get the list of ORFs
hit_strains = {'ARP7'; 'RSC58'; 'RPB7'};

% Translate genenames
hit_strains = translate(hit_strains);

% Create data
hit_data = zeros(size(hit_strains))-1;

%% Load the Tested strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Copy of YSC1057_Essential_het_diploid_obs_v5.0.xlsx', 'Essential Heterozygous Diploid');
tested_strains = data(2:end, 2);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
tested_strains(inds) = [];

% If the same strain is present more than once, remove repeats
tested_strains = unique(tested_strains);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [520];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
prescott_simmonds_2014.orfs = tested_strains;
prescott_simmonds_2014.ph = hit_data_names;
prescott_simmonds_2014.data = zeros(length(prescott_simmonds_2014.orfs),length(prescott_simmonds_2014.ph));
prescott_simmonds_2014.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, prescott_simmonds_2014.orfs);
prescott_simmonds_2014.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./prescott_simmonds_2014.mat','prescott_simmonds_2014');

%% Print out

fid = fopen('./prescott_simmonds_2014.txt','w');
write_matrix_file(fid, prescott_simmonds_2014.orfs, prescott_simmonds_2014.ph, prescott_simmonds_2014.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(prescott_simmonds_2014)
end

end