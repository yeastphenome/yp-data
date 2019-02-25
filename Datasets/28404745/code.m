%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
dalton_conibear_2017.pmid = 28404745;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(dalton_conibear_2017.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/mc-E16-11-0772-s08.xlsx', 'master a alpha results');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);

% Get the data itself
hit_data = data(2:end,4:5);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_data(find(~cellfun(@isnumeric, hit_data))) = {NaN};
hit_data = cell2mat(hit_data);

% Original z-scores: "A positive Z-score indicates mutants have reduced
% surface levels of the chimera."
hit_data = -hit_data;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16263; 16262];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
dalton_conibear_2017.orfs = hit_strains;
dalton_conibear_2017.ph = hit_data_names;
dalton_conibear_2017.data = hit_data;
dalton_conibear_2017.dataset_ids = hit_data_ids;

%% Save

save('./dalton_conibear_2017.mat','dalton_conibear_2017');

%% Print out

fid = fopen('./dalton_conibear_2017.txt','w');
write_matrix_file(fid, dalton_conibear_2017.orfs, dalton_conibear_2017.ph, dalton_conibear_2017.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(dalton_conibear_2017)
end

end

