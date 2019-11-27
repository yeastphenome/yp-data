%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
north_eide_2012.pmid = 22685415;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(north_eide_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/journal.pgen.1002699.s002.xlsx', 'Table 1');
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/journal.pgen.1002699.s003.xlsx', 'Table 1');
[FILENAMES{end+1}, data3] = read_data('xlsread','./raw_data/journal.pgen.1002699.s004.xlsx', 'Table 1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = [data1(3:end,1); data2(3:end,1); data3(3:end,1)];

% Get the data itself
hit_data = [cell2mat(data1(3:end,3:4)); cell2mat(data2(3:end,3:4)); zeros(length(data3(3:end,1)),2)];
hit_data(isnan(hit_data)) = 0;

inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16282; 16283];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
north_eide_2012.orfs = hit_strains;
north_eide_2012.ph = hit_data_names;
north_eide_2012.data = hit_data;
north_eide_2012.dataset_ids = hit_data_ids;

%% Save

save('./north_eide_2012.mat','north_eide_2012');

%% Print out

fid = fopen('./north_eide_2012.txt','w');
write_matrix_file(fid, north_eide_2012.orfs, north_eide_2012.ph, north_eide_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(north_eide_2012)
end

end

