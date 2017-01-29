%% Ando~Shima, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ando_shima_2007.pmid = 16989656;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ando_shima_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Freeze-thaw stress.xlsx');

% Get the list of ORFs and the correponding data 
hit_strains = data(4:end,1);

% Get the data itself
temp_data = data(4:end, 4);
hit_data = nan(length(hit_strains),1);
indx = find(~strcmp(temp_data, 'ND'));
hit_data(indx) = cell2mat(temp_data(indx));
hit_data = log(hit_data);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];  
hit_data(inds) = [];  

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [502];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
ando_shima_2007.orfs = hit_strains;
ando_shima_2007.ph = hit_data_names;
ando_shima_2007.data = hit_data;
ando_shima_2007.dataset_ids = hit_data_ids;

%% Save

save('./ando_shima_2007.mat','ando_shima_2007');

%% Print out

fid = fopen('./ando_shima_2007.txt','w');
write_matrix_file(fid, ando_shima_2007.orfs, ando_shima_2007.ph, ando_shima_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ando_shima_2007)
end

end

