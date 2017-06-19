%% Ohnuki~Ohya, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ohnuki_ohya_2007.pmid = 17351076;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ohnuki_ohya_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/mmc2.xlsx');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(3:end,1);

% Get the data itself
hit_data = data(3:end, 3:4); 

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YOLO57W'})) = {'YOL057W'};
hit_strains(ismember(hit_strains, {'YBRF182C-A'})) = {'YBR182C-A'};
hit_strains(ismember(hit_strains, {'YKLO72W'})) = {'YKL072W'};
hit_strains(ismember(hit_strains, {'YOLO62C'})) = {'YOL062C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];  
hit_data(inds, :) = [];

% Transform the data
hit_data = cell2mat(hit_data(:,2)) ./ cell2mat(hit_data(:,1));
hit_data = log(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11826];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
ohnuki_ohya_2007.orfs = hit_strains;
ohnuki_ohya_2007.ph = hit_data_names;
ohnuki_ohya_2007.data = hit_data;
ohnuki_ohya_2007.dataset_ids = hit_data_ids;

%% Save

save('./ohnuki_ohya_2007.mat','ohnuki_ohya_2007');

%% Print out

fid = fopen('./ohnuki_ohya_2007.txt','w');
write_matrix_file(fid, ohnuki_ohya_2007.orfs, ohnuki_ohya_2007.ph, ohnuki_ohya_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ohnuki_ohya_2007)
end

end