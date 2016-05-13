%% Choi~Basrai, 2015
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
choy_basrai_2015.pmid = 26646153;

hit_data_ids = [768];

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(choy_basrai_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load file from first test
all_data = [];
for p = 1 : 14
    [FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', ['S1 P' num2str(p)]);
    all_data = [all_data; data];
end

% Get the list of ORFs
strains = all_data(:,1);

% Clean up ORFs
strains = clean_orf(strains);

% Find anything that doesn't look like an ORF and remove it
inds = find(~is_orf(strains));
all_data(inds, :) = [];
strains(inds, :) = [];

% Get data from hits
hit_data = all_data(:,6);
indx = ~cellfun(@isnumeric, hit_data);
hit_data(indx) = {NaN};
hit_data = cell2mat(hit_data);

% Load file from second test
all_data2 = [];
for p = 1 : 14
    [FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', ['S2 P' num2str(p)]);
    all_data2 = [all_data2; data];
end

% Get the list of ORFs
strains2 = all_data2(:,1);

% Clean up ORFs
strains2 = clean_orf(strains2);

% Find anything that doesn't look like an ORF and remove it
inds = find(~is_orf(strains2));
all_data2(inds, :) = [];
strains2(inds, :) = [];

% Get data from hits
hit_data2 = all_data2(:,6);
indx = ~cellfun(@isnumeric, hit_data2);
hit_data2(indx) = {NaN};
hit_data2 = cell2mat(hit_data2);

% Combine the data and find the average
[strains, hit_data] = grpstats(hit_data, strains, {'gname','mean'});
[strains2, hit_data2] = grpstats(hit_data2, strains2, {'gname','mean'});

final_strains = unique([strains; strains2]);
final_data = nan(length(final_strains),2);

[~,ind1,ind2] = intersect(strains, final_strains);
final_data(ind2,1) = hit_data(ind1);

[~,ind1,ind2] = intersect(strains2, final_strains);
final_data(ind2,2) = hit_data2(ind1);

final_data = nanmean(final_data, 2);

final_data = 1 ./ final_data;

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);


choy_basrai_2015.orfs = final_strains;
choy_basrai_2015.ph = hit_data_names;
choy_basrai_2015.data = final_data;
choy_basrai_2015.dataset_ids = hit_dataset_ids;

%% Save

save('./choy_basrai_2015.mat','choy_basrai_2015');

fid = fopen('./choy_basrai_2015.txt','w');
write_matrix_file(fid, choy_basrai_2015.orfs, choy_basrai_2015.ph, choy_basrai_2015.data);
fclose(fid);

end
