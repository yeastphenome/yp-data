%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
brown_brown_2006.pmid = 16738548;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(brown_brown_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
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

[~,ind1,ind2] = intersect(t.Var2, data.labels_col);
hit_data_ids = t.Var1(ind1);
hit_data = hit_data(:,ind2);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
brown_brown_2006.orfs = hit_strains;
brown_brown_2006.ph = hit_data_names;
brown_brown_2006.data = hit_data;
brown_brown_2006.dataset_ids = hit_data_ids;

%% Save

save('./brown_brown_2006.mat','brown_brown_2006');

%% Print out

fid = fopen('./brown_brown_2006.txt','w');
write_matrix_file(fid, brown_brown_2006.orfs, brown_brown_2006.ph, brown_brown_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(brown_brown_2006)
end

end

