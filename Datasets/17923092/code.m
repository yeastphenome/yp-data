%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
mcclellan_frydman_2007.pmid = 17923092;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(mcclellan_frydman_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/mmc2.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data(4:end,1);
hit_strains2 = data(4:end,4);

% Get the data itself
hit_data1 = cell2mat(data(4:end,2));
hit_data2 = cell2mat(data(4:end,5));

inds = find(cellfun(@isnumeric, hit_strains1));
hit_strains1(inds) = [];
hit_data1(inds) = [];

inds = find(cellfun(@isnumeric, hit_strains2));
hit_strains2(inds) = [];
hit_data2(inds) = [];

% Eliminate all white spaces & capitalize
hit_strains1 = clean_genename(hit_strains1);
hit_strains2 = clean_genename(hit_strains2);

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);
hit_strains2 = translate(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  
hit_strains1(inds) = [];
hit_data1(inds) = [];

inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  
hit_strains2(inds) = [];
hit_data2(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = nan(length(hit_strains),2);

[~,ind1,ind2] = intersect(hit_strains, hit_strains1);
hit_data(ind1,1) = hit_data1(ind2,1);

[~,ind1,ind2] = intersect(hit_strains, hit_strains2);
hit_data(ind1,2) = hit_data2(ind2,1);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11785 11786]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
mcclellan_frydman_2007.orfs = hit_strains;
mcclellan_frydman_2007.ph = hit_data_names;
mcclellan_frydman_2007.data = hit_data;
mcclellan_frydman_2007.dataset_ids = hit_data_ids;

%% Save

save('./mcclellan_frydman_2007.mat','mcclellan_frydman_2007');

%% Print out

fid = fopen('./mcclellan_frydman_2007.txt','w');
write_matrix_file(fid, mcclellan_frydman_2007.orfs, mcclellan_frydman_2007.ph, mcclellan_frydman_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(mcclellan_frydman_2007)
end

end

