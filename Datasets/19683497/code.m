%% Fillingham~Andrews, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
fillingham_andrews_2009.pmid = 19683497;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(fillingham_andrews_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/mmc3.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
hit_strains = data.raw(:,1);

% Eliminate white spaces before/after ORF
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
data.raw(inds,:) = [];

% Get the data
hit_data = data.raw(:,2);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, hit_data)==0);
hit_data(inds) = {NaN};
hit_data = cell2mat(hit_data);

% Average data for identical ORFs that appear multiple times
[hit_strains,hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [14];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
fillingham_andrews_2009.orfs = hit_strains;
fillingham_andrews_2009.ph = hit_data_names;
fillingham_andrews_2009.data = hit_data;
fillingham_andrews_2009.dataset_ids = hit_data_ids;

%% Save

save('./fillingham_andrews_2009.mat','fillingham_andrews_2009');

%% Print out

fid = fopen('./fillingham_andrews_2009.txt','w');
write_matrix_file(fid, fillingham_andrews_2009.orfs, fillingham_andrews_2009.ph, fillingham_andrews_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(fillingham_andrews_2009)
end

end
