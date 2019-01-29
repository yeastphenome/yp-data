%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
kwon_kwak_2016.pmid = 27660105;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(kwon_kwak_2016.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table1.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);

% Get the data itself
hit_data = data(2:end,4);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_data = cell2mat(hit_data);
data(:,6) = cellfun(@strtrim, data(:,6), 'UniformOutput',0);

% Split into datasets
hit_strains_u = unique(hit_strains);
hit_data_u = zeros(length(hit_strains_u), 2);

inds = find(strcmp('HOP', data(2:end,6)));
[~,ind1,ind2] = intersect(hit_strains(inds), hit_strains_u);
hit_data_u(ind2,1) = hit_data(inds(ind1),1);

inds = find(strcmp('HIP', data(2:end,6)));
[~,ind1,ind2] = intersect(hit_strains(inds), hit_strains_u);
hit_data_u(ind2,2) = hit_data(inds(ind1),1);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16140; 16143];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
kwon_kwak_2016.orfs = hit_strains_u;
kwon_kwak_2016.ph = hit_data_names;
kwon_kwak_2016.data = hit_data_u;
kwon_kwak_2016.dataset_ids = hit_data_ids;

%% Save

save('./kwon_kwak_2016.mat','kwon_kwak_2016');

%% Print out

fid = fopen('./kwon_kwak_2016.txt','w');
write_matrix_file(fid, kwon_kwak_2016.orfs, kwon_kwak_2016.ph, kwon_kwak_2016.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(kwon_kwak_2016)
end

end

