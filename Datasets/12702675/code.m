%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
anderson_kohn_2003.pmid = 12702675;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(anderson_kohn_2003.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('textread','./raw_data/hits.txt', '%s %d','delimiter','\t');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data{1};

% Get the data itself
hit_data = data{2};
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = 71;

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
anderson_kohn_2003.orfs = hit_strains;
anderson_kohn_2003.ph = hit_data_names;
anderson_kohn_2003.data = hit_data;
anderson_kohn_2003.dataset_ids = hit_data_ids;

%% Save

save('./anderson_kohn_2003.mat','anderson_kohn_2003');

%% Print out

fid = fopen('./anderson_kohn_2003.txt','w');
write_matrix_file(fid, anderson_kohn_2003.orfs, anderson_kohn_2003.ph, anderson_kohn_2003.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(anderson_kohn_2003)
end

end

