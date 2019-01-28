%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
krol_skoneczna_2015.pmid = 25756177;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(krol_skoneczna_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('textscan','./raw_data/hits.txt', '%s %d %s');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data{1};

% Get the data itself
hit_data = data{2};

hit_collection = data{3};
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains_u = unique(hit_strains);
hit_data_u = zeros(length(hit_strains_u), 2);

inds = find(strcmp('ess', hit_collection));
[~,ind1,ind2] = intersect(hit_strains(inds), hit_strains_u);
hit_data_u(ind2,2) = hit_data(inds(ind1),1);

inds = find(~strcmp('ess', hit_collection));
[~,ind1,ind2] = intersect(hit_strains(inds), hit_strains_u);
hit_data_u(ind2,1) = hit_data(inds(ind1),1);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16151; 16152];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
krol_skoneczna_2015.orfs = hit_strains_u;
krol_skoneczna_2015.ph = hit_data_names;
krol_skoneczna_2015.data = hit_data_u;
krol_skoneczna_2015.dataset_ids = hit_data_ids;

%% Save

save('./krol_skoneczna_2015.mat','krol_skoneczna_2015');

%% Print out

fid = fopen('./krol_skoneczna_2015.txt','w');
write_matrix_file(fid, krol_skoneczna_2015.orfs, krol_skoneczna_2015.ph, krol_skoneczna_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(krol_skoneczna_2015)
end

end

