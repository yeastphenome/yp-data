%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
jo_vulpe_2007.pmid = 17785683;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(jo_vulpe_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/toxsci-07-0195-File012_kfm226.xlsx', 'Table 1');
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/toxsci-07-0195-File013_kfm226.xlsx', 'Table 1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data1(:,1);
hit_strains2 = data2(:,1);

% Get the data itself
hit_data1 = data1(:,4);
hit_data2 = data2(:,4);
   
% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);
hit_strains2 = clean_orf(hit_strains2);

inds = find(cellfun(@isnumeric, hit_strains1));
hit_strains1(inds) = [];
hit_data1(inds,:) = [];

inds = find(cellfun(@isnumeric, hit_strains2));
hit_strains2(inds) = [];
hit_data2(inds,:) = [];

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);
hit_strains2 = translate(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

hit_strains1(inds) = [];
hit_data1(inds,:) = [];

inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

hit_strains2(inds) = [];
hit_data2(inds,:) = [];

hit_data1 = cell2mat(hit_data1);
hit_data2 = cell2mat(hit_data2);

% Merge the 2 datasets
hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = nan(length(hit_strains),2);

[~,ind1,ind2] = intersect(hit_strains1, hit_strains);
hit_data(ind2,1) = hit_data1(ind1,1);
[~,ind1,ind2] = intersect(hit_strains2, hit_strains);
hit_data(ind2,2) = hit_data2(ind1,1);


% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [501; 1346];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
jo_vulpe_2007.orfs = hit_strains;
jo_vulpe_2007.ph = hit_data_names;
jo_vulpe_2007.data = hit_data;
jo_vulpe_2007.dataset_ids = hit_data_ids;

%% Save

save('./jo_vulpe_2007.mat','jo_vulpe_2007');

%% Print out

fid = fopen('./jo_vulpe_2007.txt','w');
write_matrix_file(fid, jo_vulpe_2007.orfs, jo_vulpe_2007.ph, jo_vulpe_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(jo_vulpe_2007)
end

end

