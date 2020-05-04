%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
kuroda_avalos_2019.pmid = 31734159;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(kuroda_avalos_2019.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/1-s2.0-S2405471219303825-mmc2.xlsx', '1st screen');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = [data(6:1030,1); data(6:2689,6); data(6:522,11); data(6:229,20)];

% Get the data itself
hit_data1 = cell2mat([data(6:1030,[2,4]); data(6:2689,[7,9]); data(6:522,[12,14])]); 
hit_data2 = horzcat(cell2mat(data(6:229,21)), zeros(229-6+1,1)+nan);
hit_data = [hit_data1; hit_data2];
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'VER093C-A'})) = {'YER093C-A'};
hit_strains(ismember(hit_strains, {'YMROB4W'})) = {'YMR084W'};
hit_strains(ismember(hit_strains, {'YARD02C-A'})) = {'YAR002C-A'};
hit_strains(ismember(hit_strains, {'YFLOOIW'})) = {'YFL001W'};
hit_strains(ismember(hit_strains, {'YMLOIOC-B'})) = {'YML010C-B'};
hit_strains(ismember(hit_strains, {'VER091C'})) = {'YER091C'};
hit_strains(ismember(hit_strains, {'VCR086W'})) = {'YCR086W'};
hit_strains(ismember(hit_strains, {'YNLO15W'})) = {'YNL015W'};
hit_strains(ismember(hit_strains, {'VER064C'})) = {'YER064C'};
hit_strains(ismember(hit_strains, {'VBR285W'})) = {'YBR285W'};
hit_strains(ismember(hit_strains, {'YMR08IC'})) = {'YMR081C'};
hit_strains(ismember(hit_strains, {'YJRIOOC'})) = {'YJR100C'};
hit_strains(ismember(hit_strains, {'YGR2B8W'})) = {'YGR288W'};
hit_strains(ismember(hit_strains, {'YPROT4C'})) = {'YPR014C'};
hit_strains(ismember(hit_strains, {'VCR087W'})) = {'YCR087W'};
hit_strains(ismember(hit_strains, {'YARD50W'})) = {'YAR050W'};
hit_strains(ismember(hit_strains, {'VCL075W'})) = {'YCL075W'};
hit_strains(ismember(hit_strains, {'YJR09IC'})) = {'YJR091C'};
hit_strains(ismember(hit_strains, {'YJR044C6'})) = {'YJR044C'};
hit_strains(ismember(hit_strains, {'VCR031C'})) = {'YCR031C'};

% Couldn't fix one remaining gene unequivocably: YLR25TW (T = 1 or 7?)

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds)); 


% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16414; 16411];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
kuroda_avalos_2019.orfs = hit_strains;
kuroda_avalos_2019.ph = hit_data_names;
kuroda_avalos_2019.data = hit_data;
kuroda_avalos_2019.dataset_ids = hit_data_ids;

%% Save

save('./kuroda_avalos_2019.mat','kuroda_avalos_2019');

%% Print out

fid = fopen('./kuroda_avalos_2019.txt','w');
write_matrix_file(fid, kuroda_avalos_2019.orfs, kuroda_avalos_2019.ph, kuroda_avalos_2019.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(kuroda_avalos_2019)
end

end

