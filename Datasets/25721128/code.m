%% Hendry~Brown, 2015
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
hendry_brown_2015.pmid = 25721128;

phenotypes = {'Rnr3 abundance'};
treatments = {'Untreated conditions'; 'YPD, 0.03% MMS'};

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hendry_brown_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread', './raw_data/TableS1.xlsx', 'S1');
[FILENAMES{end+1}, data2] = read_data('xlsread', './raw_data/TableS2.xlsx', 'S2');

% Get the list of ORFs
hit_strains = data(11:end, 5);
hit_strains2 = data2(11:end, 5);

% Clean up ORFs
hit_strains = clean_orf(hit_strains);
hit_strains2 = clean_orf(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds)); 
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));

% Retrieve data from files
hit_data = cell2mat(data(11:end, 8));
hit_data2 = cell2mat(data2(11:end, 8));

% Average any repeated value
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

% Combined the 2 datasets
hit_strains3 = unique([hit_strains; hit_strains2]);
hit_data3 = nan(length(hit_strains3),2);

[~,ind1,ind2] = intersect(hit_strains3, hit_strains);
hit_data3(ind1,1) = hit_data(ind2);

[~,ind1,ind2] = intersect(hit_strains3, hit_strains2);
hit_data3(ind1,2) = hit_data2(ind2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [693;606];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

hendry_brown_2015.orfs = hit_strains3;
hendry_brown_2015.ph = hit_data_names;
hendry_brown_2015.data = hit_data3;
hendry_brown_2015.dataset_ids = hit_data_ids;


%% Save

save('./hendry_brown_2015.mat','hendry_brown_2015');

fid = fopen('./hendry_brown_2015.txt','w');
write_matrix_file(fid, hendry_brown_2015.orfs, hendry_brown_2015.ph, hendry_brown_2015.data);
fclose(fid);

end
