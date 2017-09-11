%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
thorpe_dawes_2004.pmid = 15087496;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(thorpe_dawes_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/05888Table2.xlsx', 'Total');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,2);

% Get the data itself
hit_data = data(2:end,4:2:12);
hit_data2 = data(2:end,5:2:13);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_data(cellfun(@isnumeric, hit_data)) = {0};
hit_data(find(strcmpi('S', hit_data))) = {-1};
hit_data(find(strcmpi('R', hit_data))) = {1};

hit_data(~cellfun(@isnumeric, hit_data)) = {NaN};
hit_data = cell2mat(hit_data);

% Coefficients
hit_data2 = cell2mat(hit_data2);
hit_data2 = hit_data2 + 1;
hit_data2(isnan(hit_data2)) = 1;

hit_data = hit_data .* hit_data2;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [4959 4957 488 4954 4953]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
thorpe_dawes_2004.orfs = hit_strains;
thorpe_dawes_2004.ph = hit_data_names;
thorpe_dawes_2004.data = hit_data;
thorpe_dawes_2004.dataset_ids = hit_data_ids;

%% Save

save('./thorpe_dawes_2004.mat','thorpe_dawes_2004');

%% Print out

fid = fopen('./thorpe_dawes_2004.txt','w');
write_matrix_file(fid, thorpe_dawes_2004.orfs, thorpe_dawes_2004.ph, thorpe_dawes_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(thorpe_dawes_2004)
end

end

