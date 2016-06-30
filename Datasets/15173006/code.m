%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
wu_brown_2004.pmid = 15173006;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(wu_brown_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('read_matrix_file','./raw_data/inline-supplementary-material-1.txt', 4, 1);

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data.labels_row(:,1);
hit_data = data.data;
  
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
[FILENAMES{end+1}, t] = read_data('readtable','./extras/dataset_id_name.txt','delimiter','\t','ReadVariableNames',false);

[~,ind1,ind2] = intersect(t.Var2, upper(data.labels_col));
hit_data_ids = t.Var1(ind1);
hit_data = hit_data(:,ind1);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
wu_brown_2004.orfs = hit_strains;
wu_brown_2004.ph = hit_data_names;
wu_brown_2004.data = hit_data;
wu_brown_2004.dataset_ids = hit_data_ids;

%% Save

save('./wu_brown_2004.mat','wu_brown_2004');

%% Print out

fid = fopen('./wu_brown_2004.txt','w');
write_matrix_file(fid, wu_brown_2004.orfs, wu_brown_2004.ph, wu_brown_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(wu_brown_2004)
end

end

