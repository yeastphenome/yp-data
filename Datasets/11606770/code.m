%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
birrell_brown_2001.pmid = 11606770;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(birrell_brown_2001.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('read_matrix_file','./raw_data/inline-supplementary-material-1.txt', 4, 1);

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data.labels_row(:,1);

ind = find(ismember(data.labels_col, {'UVB','UVC'}));
data.labels_col = data.labels_col(ind);
hit_data = data.data(:,ind);
  
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
hit_data_ids = [745; 466];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
birrell_brown_2001.orfs = hit_strains;
birrell_brown_2001.ph = hit_data_names;
birrell_brown_2001.data = hit_data;
birrell_brown_2001.dataset_ids = hit_data_ids;

%% Save

save('./birrell_brown_2001.mat','birrell_brown_2001');

%% Print out

fid = fopen('./birrell_brown_2001.txt','w');
write_matrix_file(fid, birrell_brown_2001.orfs, birrell_brown_2001.ph, birrell_brown_2001.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(birrell_brown_2001)
end

end

