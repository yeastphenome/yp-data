%% Zatorska~Strahl, 2017
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
zatorska_strahl_2017.pmid = 28598353;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(zatorska_strahl_2017.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/ijms-18-01226-s001.xlsx', 'Table S1');

% Get the list of ORFs and the correponding data 
hit_strains = data(4:end,1);

% Get the data itself
hit_data = cell2mat(data(4:end, 5)); 
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Remove the DaMP trials
types = data(4:end, 6);
inds = find(ismember(types, {'DAmP'}));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% Find anything that doesn't look like an ORF
hit_strains(ismember(hit_strains, {'YLR287-A'})) = {'YLR287C-A'};
hit_strains(ismember(hit_strains, {'YBRF182C-A'})) = {'YBR182C-A'};
hit_strains(ismember(hit_strains, {'YOLO57W'})) = {'YOL057W'};
hit_strains(ismember(hit_strains, {'YOLO62C'})) = {'YOL062C'};
hit_strains(ismember(hit_strains, {'YBRF182C-A'})) = {'YBR182C-A'};
hit_strains(ismember(hit_strains, {'YKLO72W'})) = {'YKL072W'};
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];


% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11862];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
zatorska_strahl_2017.orfs = hit_strains;
zatorska_strahl_2017.ph = hit_data_names;
zatorska_strahl_2017.data = hit_data;
zatorska_strahl_2017.dataset_ids = hit_data_ids;

%% Save

save('./zatorska_strahl_2017.mat','zatorska_strahl_2017');

%% Print out

fid = fopen('./zatorska_strahl_2017.txt','w');
write_matrix_file(fid, zatorska_strahl_2017.orfs, zatorska_strahl_2017.ph, zatorska_strahl_2017.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(zatorska_strahl_2017)
end

end