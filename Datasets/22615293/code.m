%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hoepfner_parker_2012.pmid = 22615293;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hoepfner_parker_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

t = {'Cpd1','Cpd2','Cpd3','Cpd4','Vori'};
compound = cell(length(t),1);

for i = 1:5
    filename = ['./raw_data/HIP-HOP Scores ', t{i}, '.txt'];
    [FILENAMES{end+1}, data] = read_data('readtable', filename, 'Delimiter', '\t', 'HeaderLines', 0);
    
    data.SYSTEMATIC_NAME = clean_orf(data.SYSTEMATIC_NAME);
    data.SYSTEMATIC_NAME = translate(data.SYSTEMATIC_NAME);
    
    if i == 1
        data_all = data;
    else
        data_all = [data_all; data];
    end
    
    compound{i} = num2str(unique(data.COMPOUND));
end

compound_concentrations = unique(data_all.COMPOUND_CONCENTRATION);

hit_strains = unique(data_all.SYSTEMATIC_NAME);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% Get the data itself
hit_data = nan(length(hit_strains), length(compound_concentrations)*2);
compound_concentrations_all = [compound_concentrations; compound_concentrations];
hip_hop = cell(length(compound_concentrations)*2,1);
hip_hop(1:length(compound_concentrations)) = {'HOP'};
hip_hop(length(compound_concentrations)+1:end) = {'HIP'};

for i = 1:length(compound_concentrations)
    inds = find(strcmp(compound_concentrations{i}, data_all.COMPOUND_CONCENTRATION) & strcmp('HOP', data_all.EXPERIMENT_TYPE));
    [~,ind1,ind2] = intersect(data_all.SYSTEMATIC_NAME(inds), hit_strains);
    hit_data(ind2,i) = data_all.SCORE(inds(ind1));
    
    inds = find(strcmp(compound_concentrations{i}, data_all.COMPOUND_CONCENTRATION) & strcmp('HIP', data_all.EXPERIMENT_TYPE));
    [~,ind1,ind2] = intersect(data_all.SYSTEMATIC_NAME(inds), hit_strains);
    hit_data(ind2,length(compound_concentrations)+i) = data_all.SCORE(inds(ind1));
end

missing_vals = find(sum(~isnan(hit_data),1)==0);
compound_concentrations_all(missing_vals) = [];
hip_hop(missing_vals) = [];
hit_data(:, missing_vals) = [];

tmp = split(compound_concentrations_all, '_');
compounds_all = tmp(:,1);
doses_all = tmp(:,2);
doses_all = cellfun(@str2num, doses_all);

compound_all_names = cell(length(compounds_all),1);
for i = 1 : length(compounds_all)
    ind = find(strcmp(compounds_all{i}, compound));
    compound_all_names{i} = t{ind};
end
   
[num, txt, raw] = xlsread('./raw_data/doses_datasetids.xlsx','Sheet1');

hit_data_ids = nan(length(compound_all_names),1);
for i = 1 : length(compound_all_names)
    ind = find(strcmp(compound_all_names{i}, raw(:,1)) & strcmp(hip_hop{i}, raw(:,2)) & (cell2mat(raw(:,3)) == doses_all(i)));
    hit_data_ids(i) = raw{ind,4};
end

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
hoepfner_parker_2012.orfs = hit_strains;
hoepfner_parker_2012.ph = hit_data_names;
hoepfner_parker_2012.data = hit_data;
hoepfner_parker_2012.dataset_ids = hit_data_ids;

%% Save

save('./hoepfner_parker_2012.mat','hoepfner_parker_2012');

%% Print out

fid = fopen('./hoepfner_parker_2012.txt','w');
write_matrix_file(fid, hoepfner_parker_2012.orfs, hoepfner_parker_2012.ph, hoepfner_parker_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(hoepfner_parker_2012)
end

end

