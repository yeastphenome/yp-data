%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
segura_wang_korbel_2017.pmid = 28818866;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(segura_wang_korbel_2017.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS5.xlsx', 'Tabelle1');

drugs = {'Campt','Doxo','HU','MMS'};

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);
hit_data = data(:,3:4);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {NaN};
hit_data = cell2mat(hit_data);
hit_data = (hit_data(:,1) - hit_data(:,2))./hit_data(:,2);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds)); 

% Get the data itself
hit_strains_u = unique(hit_strains);
hit_data_u = zeros(length(hit_strains_u), length(drugs));

drugs_inds = [];
for i = 1 : length(drugs)
    drugs_inds(i) = find(strcmp(upper(drugs{i}), hit_strains));
end

for i = 1 : length(drugs)
    inds_start = drugs_inds(i)+1;
    inds_end = length(hit_strains);
    if i < length(drugs)
        inds_end = drugs_inds(i+1);
    end
    hit_strains_this = hit_strains(inds_start:inds_end);
    hit_data_this = hit_data(inds_start:inds_end,1);
    
    [~,ind1,ind2] = intersect(hit_strains_this, hit_strains_u);
    hit_data_u(ind2,i) = hit_data_this(ind1,1);
end

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_u));
disp(hit_strains_u(inds)); 

hit_strains_u(inds) = [];
hit_data_u(inds,:) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16150; 16149; 16147; 16148];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
segura_wang_korbel_2017.orfs = hit_strains_u;
segura_wang_korbel_2017.ph = hit_data_names;
segura_wang_korbel_2017.data = hit_data_u;
segura_wang_korbel_2017.dataset_ids = hit_data_ids;

%% Save

save('./segura_wang_korbel_2017.mat','segura_wang_korbel_2017');

%% Print out

fid = fopen('./segura_wang_korbel_2017.txt','w');
write_matrix_file(fid, segura_wang_korbel_2017.orfs, segura_wang_korbel_2017.ph, segura_wang_korbel_2017.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(segura_wang_korbel_2017)
end

end

