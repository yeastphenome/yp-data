%% Junne~Hoepfner, 2015
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
junne_hoepfner_2015.pmid = 25616894;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(junne_hoepfner_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data
% First two datasets

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Junnet_et_al_HIPHOPrawdata.xlsx', 'Decatransin Cpd 1 Exp 1');

% Get the list of ORFs and the correponding data 
all_strains = data(2:end,5);

% separate hom from het & retrieve data
hom_indx = find(~cellfun(@isempty, (strfind(data(2:end,1), 'HOP'))));
het_indx = find(~cellfun(@isempty, (strfind(data(2:end,1), 'HIP'))));
hom_strains_1 = all_strains(hom_indx);
het_strains_1 = all_strains(het_indx);
hom_data_1 = data(2:end, 7);
hom_data_1 = cell2mat(hom_data_1(hom_indx));
het_data_1 = data(2:end, 7);
het_data_1 = cell2mat(het_data_1(het_indx));
   
% Eliminate all white spaces & capitalize
hom_strains_1 = clean_orf(hom_strains_1);
het_strains_1 = clean_orf(het_strains_1);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hom_strains_1));
hom_strains_1(inds) = [];
hom_data_1(inds) = [];

inds = find(~is_orf(het_strains_1));
het_strains_1(inds) = [];
het_data_1(inds) = [];

% Second set of data
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Junnet_et_al_HIPHOPrawdata.xlsx', 'Decatransin Cpd 2 Exp 2');

% Get the list of ORFs and the correponding data 
all_strains = data(2:end,5);

% separate hom from het & retrieve data
hom_indx = find(~cellfun(@isempty, (strfind(data(2:end,1), 'HOP'))));
het_indx = find(~cellfun(@isempty, (strfind(data(2:end,1), 'HIP'))));
hom_strains_2 = all_strains(hom_indx);
het_strains_2 = all_strains(het_indx);
hom_data_2 = data(2:end, 7);
hom_data_2 = cell2mat(hom_data_2(hom_indx));
het_data_2 = data(2:end, 7);
het_data_2 = cell2mat(het_data_2(het_indx));
   
% Eliminate all white spaces & capitalize
hom_strains_2 = clean_orf(hom_strains_2);
het_strains_2 = clean_orf(het_strains_2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hom_strains_2));
hom_strains_2(inds) = [];
hom_data_2(inds) = [];

inds = find(~is_orf(het_strains_2));
het_strains_2(inds) = [];
het_data_2(inds) = [];

% Third set of data
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Junnet_et_al_HIPHOPrawdata.xlsx', 'Cotransin Cpd 2');

% Get the list of ORFs and the correponding data 
all_strains = data(2:end,5);

% separate hom from het & retrieve data
hom_indx = find(~cellfun(@isempty, (strfind(data(2:end,1), 'HOP'))));
het_indx = find(~cellfun(@isempty, (strfind(data(2:end,1), 'HIP'))));
hom_strains_3 = all_strains(hom_indx);
het_strains_3 = all_strains(het_indx);
hom_data_3 = data(2:end, 7);
hom_data_3 = cell2mat(hom_data_3(hom_indx));
het_data_3 = data(2:end, 7);
het_data_3 = cell2mat(het_data_3(het_indx));
   
% Eliminate all white spaces & capitalize
hom_strains_3 = clean_orf(hom_strains_3);
het_strains_3 = clean_orf(het_strains_3);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hom_strains_3));
hom_strains_3(inds) = [];
hom_data_3(inds) = [];

inds = find(~is_orf(het_strains_3));
het_strains_3(inds) = [];
het_data_3(inds) = [];

% Fourth set of data
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Junnet_et_al_HIPHOPrawdata.xlsx', 'Cotansin Cpd 3 (HUN-7293)');

% Get the list of ORFs and the correponding data 
all_strains = data(2:end,5);

% separate hom from het & retrieve data
hom_indx = find(~cellfun(@isempty, (strfind(data(2:end,1), 'HOP'))));
het_indx = find(~cellfun(@isempty, (strfind(data(2:end,1), 'HIP'))));
hom_strains_4 = all_strains(hom_indx);
het_strains_4 = all_strains(het_indx);
hom_data_4 = data(2:end, 7);
hom_data_4 = cell2mat(hom_data_4(hom_indx));
het_data_4 = data(2:end, 7);
het_data_4 = cell2mat(het_data_4(het_indx));
   
% Eliminate all white spaces & capitalize
hom_strains_4 = clean_orf(hom_strains_4);
het_strains_4 = clean_orf(het_strains_4);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hom_strains_4));
hom_strains_4(inds) = [];
hom_data_4(inds) = [];

inds = find(~is_orf(het_strains_4));
het_strains_4(inds) = [];
het_data_4(inds) = [];

% Combine the data into one matrix
all_strains = [hom_strains_1; het_strains_1; hom_strains_2; het_strains_2; hom_strains_3; het_strains_3; hom_strains_4; het_strains_4];
all_strains = unique(all_strains);
all_data = nan(length(all_strains), 8);

[~, ind1, ind2] = intersect(all_strains, hom_strains_1);
all_data(ind1, 1) = hom_data_1(ind2);

[~, ind1, ind2] = intersect(all_strains, het_strains_1);
all_data(ind1, 2) = het_data_1(ind2);

[~, ind1, ind2] = intersect(all_strains, hom_strains_2);
all_data(ind1, 3) = hom_data_2(ind2);

[~, ind1, ind2] = intersect(all_strains, het_strains_2);
all_data(ind1, 4) = het_data_2(ind2);

[~, ind1, ind2] = intersect(all_strains, hom_strains_3);
all_data(ind1, 5) = hom_data_3(ind2);

[~, ind1, ind2] = intersect(all_strains, het_strains_3);
all_data(ind1, 6) = het_data_3(ind2);

[~, ind1, ind2] = intersect(all_strains, hom_strains_4);
all_data(ind1, 7) = hom_data_4(ind2);

[~, ind1, ind2] = intersect(all_strains, het_strains_4);
all_data(ind1, 8) = het_data_4(ind2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [772; 773; 5262; 5263; 5264; 5265; 5266; 5267];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
junne_hoepfner_2015.orfs = all_strains;
junne_hoepfner_2015.ph = hit_data_names;
junne_hoepfner_2015.data = all_data;
junne_hoepfner_2015.dataset_ids = hit_data_ids;

%% Save

save('./junne_hoepfner_2015.mat','junne_hoepfner_2015');

%% Print out

fid = fopen('./junne_hoepfner_2015.txt','w');
write_matrix_file(fid, junne_hoepfner_2015.orfs, junne_hoepfner_2015.ph, junne_hoepfner_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(junne_hoepfner_2015)
end

end

