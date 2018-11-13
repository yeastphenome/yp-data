%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
zhu_li_2015.pmid = 25823586;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(zhu_li_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Zhu el al G3_qCTF_HI Summary.xlsx', 'Sheet1');

% Eliminate the missing values
data(cellfun(@isnumeric, data(:,3)),:) = [];

data(1,:) = [];

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,3);

% Get the data itself
inds = find(~cellfun(@isnumeric, data(:,4)));
data(inds,4) = {NaN};
hit_data = cell2mat(data(:,4));

% Average all the controls
inds = find(strncmp(hit_strains, 'Control', 7));
ctrl_avg = nanmean(hit_data(inds));

% Normalize by the controls
hit_data = hit_data / ctrl_avg;
 
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YLR287-A'})) = {'YLR287C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If not possible, eliminate the entry
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16006];


%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
zhu_li_2015.orfs = hit_strains;
zhu_li_2015.ph = hit_data_names;
zhu_li_2015.data = hit_data;
zhu_li_2015.dataset_ids = hit_data_ids;

%% Save

save('./zhu_li_2015.mat','zhu_li_2015');

%% Print out

fid = fopen('./zhu_li_2015.txt','w');
write_matrix_file(fid, zhu_li_2015.orfs, zhu_li_2015.ph, zhu_li_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(zhu_li_2015)
end

end

