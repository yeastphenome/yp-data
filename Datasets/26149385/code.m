%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
sauerwald_rapaport_2015.pmid = 26149385;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(sauerwald_rapaport_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/zmb999100949sd2.xlsx', 'Supplementary Table 1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data1(2:end,2);

% Get the data itself
hit_data1 = -ones(size(hit_strains1));
   
% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);

inds = find(cellfun(@isnumeric, hit_strains1));
hit_strains1(inds) = [];
hit_data1(inds,:) = [];

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

% If the same strain is present more than once, average its values
[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});

%% Dataset #2

[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/zmb999100949sd3.xlsx', 'Supplemental Table 2');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains2 = data2(2:end,2);

% Get the data itself
hit_data2 = cell2mat(data2(2:end,1))-1;
   
% Eliminate all white spaces & capitalize
hit_strains2 = clean_orf(hit_strains2);

% Fix typos
hit_strains2(ismember(hit_strains2, {'YOLO62C'})) = {'YOL062C'};
hit_strains2(ismember(hit_strains2, {'YOLO57W'})) = {'YOL057W'};

% If in gene name form, transform into ORF name
hit_strains2 = translate(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

% If the same strain is present more than once, average its values
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

%% Merge

hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = zeros(length(hit_strains),2);

[~,ind1,ind2] = intersect(hit_strains1, hit_strains);
hit_data(ind2,1) = hit_data1(ind1,1);
[~,ind1,ind2] = intersect(hit_strains2, hit_strains);
hit_data(ind2,2) = hit_data2(ind1,1);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11865; 11866];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
sauerwald_rapaport_2015.orfs = hit_strains;
sauerwald_rapaport_2015.ph = hit_data_names;
sauerwald_rapaport_2015.data = hit_data;
sauerwald_rapaport_2015.dataset_ids = hit_data_ids;

%% Save

save('./sauerwald_rapaport_2015.mat','sauerwald_rapaport_2015');

%% Print out

fid = fopen('./sauerwald_rapaport_2015.txt','w');
write_matrix_file(fid, sauerwald_rapaport_2015.orfs, sauerwald_rapaport_2015.ph, sauerwald_rapaport_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(sauerwald_rapaport_2015)
end

end

