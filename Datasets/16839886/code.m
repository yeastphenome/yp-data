%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
slessareva_dohlman_2006.pmid = 16839886;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(slessareva_dohlman_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/hits.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);

% Get the data itself
hit_data = data(2:end,4);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_data2 = regexp(hit_data, '[0-9\.]*(?=% )','match');
hit_data2 = vertcat(hit_data2{:});
hit_data2 = cellfun(@str2num, hit_data2);
hit_data2 = hit_data2/100;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data2, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5236];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
slessareva_dohlman_2006.orfs = hit_strains;
slessareva_dohlman_2006.ph = hit_data_names;
slessareva_dohlman_2006.data = hit_data;
slessareva_dohlman_2006.dataset_ids = hit_data_ids;

%% Save

save('./slessareva_dohlman_2006.mat','slessareva_dohlman_2006');

%% Print out

fid = fopen('./slessareva_dohlman_2006.txt','w');
write_matrix_file(fid, slessareva_dohlman_2006.orfs, slessareva_dohlman_2006.ph, slessareva_dohlman_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(slessareva_dohlman_2006)
end

end

