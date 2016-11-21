%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
lee_giaever_2005.pmid = 16121259;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(lee_giaever_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/pgen.0010024.sd001.xlsx', 'OrfGeneData');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);
hit_screens = data(1,8:3:end)';

% Get the data itself
hit_data = cell2mat(data(2:end,8:3:end));
   
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
[FILENAMES{end+1}, screens] = read_data('readtable','./extras/screen_datasetids.txt', 'delimiter','\t','ReadVariableNames',0);
[~,ind1,ind2] = intersect(screens.Var1, hit_screens);

hit_data_ids(ind2) = screens.Var2(ind1);
hit_data_ids = hit_data_ids';

[hit_data_ids, hit_data] = grpstats(hit_data', hit_data_ids, {'gname','mean'});
hit_data = hit_data';
hit_data_ids = cellfun(@str2num, hit_data_ids);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
lee_giaever_2005.orfs = hit_strains;
lee_giaever_2005.ph = hit_data_names;
lee_giaever_2005.data = hit_data;
lee_giaever_2005.dataset_ids = hit_data_ids;

%% Save

save('./lee_giaever_2005.mat','lee_giaever_2005');

%% Print out

fid = fopen('./lee_giaever_2005.txt','w');
write_matrix_file(fid, lee_giaever_2005.orfs, lee_giaever_2005.ph, lee_giaever_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(lee_giaever_2005)
end

end

