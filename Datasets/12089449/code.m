%% Jorgensen~Tyers, 2002
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
jorgensen_tyers_2002.pmid = 12089449;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(jorgensen_tyers_2002.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};


%% Hap strains

% Load hit strains
[FILENAMES{end+1}, data_hap] = read_data('xlsread','./raw_data/data.xlsx', 'haps');
[FILENAMES{end+1}, data_het] = read_data('xlsread','./raw_data/data.xlsx', 'hets');

% Get the list of ORFs and the correponding data 
hap_strains = data_hap(:,1);
hap_data = data_hap(:,3); 

het_strains = data_het(:,1);
het_data = data_het(:,2); 
    
% Eliminate all white spaces & capitalize
hap_strains = clean_orf(hap_strains);
het_strains = clean_orf(het_strains);

inds = find(cellfun(@isnumeric, het_strains));
het_strains(inds) = [];
het_data(inds,:) = [];

hap_strains = translate(hap_strains);
het_strains = translate(het_strains);

% Normalize data to WT
hap_data(~cellfun(@isnumeric, hap_data)) = {NaN};
hap_data = cell2mat(hap_data);

het_data(~cellfun(@isnumeric, het_data)) = {NaN};
het_data = cell2mat(het_data);

ind_wt = find(strcmp('WT', hap_strains));
hap_data = hap_data ./ repmat(hap_data(ind_wt,:),length(hap_data),1) - 1;

ind_wt = find(strcmp('NA', het_strains));  % Seems to be the WT based on the description of the data in the paper
het_data = het_data ./ repmat(het_data(ind_wt,:),length(het_data),1) - 1;

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hap_strains));
disp(hap_strains(inds));  

hap_strains(inds) = [];
hap_data(inds,:) = [];

inds = find(~is_orf(het_strains));
disp(het_strains(inds));

het_strains(inds) = [];
het_data(inds,:) = [];

% Merge the 2 datasets

hit_strains = unique([hap_strains; het_strains]);
hit_data = nan(length(hit_strains),2);

[~,ind1,ind2] = intersect(hap_strains, hit_strains);
hit_data(ind2,1) = hap_data(ind1);
[~,ind1,ind2] = intersect(het_strains, hit_strains);
hit_data(ind2,2) = het_data(ind1);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

hit_data_ids = [68; 709];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
jorgensen_tyers_2002.orfs = hit_strains;
jorgensen_tyers_2002.ph = hit_data_names;
jorgensen_tyers_2002.data = hit_data;
jorgensen_tyers_2002.dataset_ids = hit_data_ids;


%% Save

save('./jorgensen_tyers_2002.mat','jorgensen_tyers_2002');

%% Print out

fid = fopen('./jorgensen_tyers_2002.txt','w');
write_matrix_file(fid, jorgensen_tyers_2002.orfs, jorgensen_tyers_2002.ph, jorgensen_tyers_2002.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(jorgensen_tyers_2002)
end

end

