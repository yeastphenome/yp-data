%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
north_vulpe_2011.pmid = 21912624;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(north_vulpe_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/s007.xlsx', 'Sheet1');
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/s008.xlsx', 'Sheet1');
[FILENAMES{end+1}, data3] = read_data('xlsread','./raw_data/s009.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data1(4:end,1);
hit_strains2 = data2(4:end,1);
hit_strains3 = data3(4:end,1);

% Get the data itself
hit_data1 = data1(4:end,3:8);
hit_data2 = data2(4:end,3:8);
hit_data3 = data3(4:end,3:8);

% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);
hit_strains2 = clean_orf(hit_strains2);
hit_strains3 = clean_orf(hit_strains3);

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);
hit_strains2 = translate(hit_strains2);
hit_strains3 = translate(hit_strains3);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  
inds = find(~is_orf(hit_strains3));
disp(hit_strains3(inds));  

hit_data1 = cell2mat(hit_data1);
hit_data2 = cell2mat(hit_data2);
hit_data3 = cell2mat(hit_data3);

hit_data1(isnan(hit_data1)) = 0;
hit_data2(isnan(hit_data2)) = 0;
hit_data3(isnan(hit_data3)) = 0;

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids1 = [513 1330 1331 775 1332 1333]';
hit_data_ids2 = [1334 1335 1336 1337 1338 1339]';
hit_data_ids3 = [1340 1342 1341 1343 1345 1344]';
hit_data_ids = [hit_data_ids1; hit_data_ids2; hit_data_ids3];

hit_strains = unique([hit_strains1; hit_strains2; hit_strains3]);

[~,ind1,ind2] = intersect(hit_strains1, hit_strains);
t1 = zeros(length(hit_strains), length(hit_data_ids1));
t1(ind2,:) = hit_data1(ind1,:);

[~,ind1,ind2] = intersect(hit_strains2, hit_strains);
t2 = zeros(length(hit_strains), length(hit_data_ids2));
t2(ind2,:) = hit_data2(ind1,:);

[~,ind1,ind2] = intersect(hit_strains3, hit_strains);
t3 = zeros(length(hit_strains), length(hit_data_ids3));
t3(ind2,:) = hit_data3(ind1,:);

hit_data = [t1 t2 t3];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
north_vulpe_2011.orfs = hit_strains;
north_vulpe_2011.ph = hit_data_names;
north_vulpe_2011.data = hit_data;
north_vulpe_2011.dataset_ids = hit_data_ids;

%% Save

save('./north_vulpe_2011.mat','north_vulpe_2011');

%% Print out

fid = fopen('./north_vulpe_2011.txt','w');
write_matrix_file(fid, north_vulpe_2011.orfs, north_vulpe_2011.ph, north_vulpe_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(north_vulpe_2011)
end

end

