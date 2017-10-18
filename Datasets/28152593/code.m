%% Lesbesgue~Lemeer, 2017
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
lebesgue_lemeer_2017.pmid = 28152593;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(lebesgue_lemeer_2017.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/pr6b00691_si_001.xlsx', 'S5-Deletion-Overexpression ');

% Get the list of ORFs and the correponding data 
hit_strains = data(3:end,1);

% Get the data itself
hit_data = cell2mat(data(3:end,3));

% Remove numeric indices
inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YBRF182C-A'})) = {'YBR182C-A'};
hit_strains(ismember(hit_strains, {'YJL206-A'})) = {'YJL206C-A'};
hit_strains(ismember(hit_strains, {'YKLO72W'})) = {'YKL072W'};
hit_strains(ismember(hit_strains, {'YLR287-A'})) = {'YLR287C-A'};
hit_strains(ismember(hit_strains, {'YOLO57W'})) = {'YOL057W'};
hit_strains(ismember(hit_strains, {'YOLO62C'})) = {'YOL062C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11864];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
lebesgue_lemeer_2017.orfs = hit_strains;
lebesgue_lemeer_2017.ph = hit_data_names;
lebesgue_lemeer_2017.data = hit_data;
lebesgue_lemeer_2017.dataset_ids = hit_data_ids;

%% Save

save('./lebesgue_lemeer_2017.mat','lebesgue_lemeer_2017');

%% Print out

fid = fopen('./lebesgue_lemeer_2017.txt','w');
write_matrix_file(fid, lebesgue_lemeer_2017.orfs, lebesgue_lemeer_2017.ph, lebesgue_lemeer_2017.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(lebesgue_lemeer_2017)
end

end