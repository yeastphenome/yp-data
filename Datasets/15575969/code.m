%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
haugen_vanhouten_2004.pmid = 15575969;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(haugen_vanhouten_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/gb-2004-5-12-r95-s5.xlsx', 'Query2');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);

% Get the data itself
hit_data = data(:,26);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));

hit_strains(inds) = [];
hit_data(inds) = [];

hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [485];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
haugen_vanhouten_2004.orfs = hit_strains;
haugen_vanhouten_2004.ph = hit_data_names;
haugen_vanhouten_2004.data = hit_data;
haugen_vanhouten_2004.dataset_ids = hit_data_ids;

%% Save

save('./haugen_vanhouten_2004.mat','haugen_vanhouten_2004');

%% Print out

fid = fopen('./haugen_vanhouten_2004.txt','w');
write_matrix_file(fid, haugen_vanhouten_2004.orfs, haugen_vanhouten_2004.ph, haugen_vanhouten_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(haugen_vanhouten_2004)
end

end
