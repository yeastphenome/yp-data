%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
schluter_conibear_2008.pmid = 18216282;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(schluter_conibear_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/DataSet1.xlsx', 'Table S1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);

% Get the data itself
hit_data = data(:,3:8);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YLR287-A'})) = {'YLR287C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16030 16030 16031 16031 16032 16032]';
[hit_data_ids, hit_data] = grpstats(hit_data', hit_data_ids, {'gname','mean'});
hit_data = hit_data';
hit_data_ids = cellfun(@str2num, hit_data_ids);


%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
schluter_conibear_2008.orfs = hit_strains;
schluter_conibear_2008.ph = hit_data_names;
schluter_conibear_2008.data = hit_data;
schluter_conibear_2008.dataset_ids = hit_data_ids;

%% Save

save('./schluter_conibear_2008.mat','schluter_conibear_2008');

%% Print out

fid = fopen('./schluter_conibear_2008.txt','w');
write_matrix_file(fid, schluter_conibear_2008.orfs, schluter_conibear_2008.ph, schluter_conibear_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(schluter_conibear_2008)
end

end

