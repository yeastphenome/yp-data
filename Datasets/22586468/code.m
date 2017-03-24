%% Peyroche~Plateau, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
peyroche_plateau_2012.pmid = 22586468;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(peyroche_plateau_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data 

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/journal.pone.0036343.s004.xlsx', 'data');

% Get indices of the data columns
ind_data = 5:8;

% Eliminate anything that doesn't look like an ORF
inds = find(~is_orf(data.raw(:,1)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,1) = clean_orf(data.raw(:,1));

% Retrieve the data from the data set
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% If in gene name form, transform into ORF name
data.raw(:,1) = translate(data.raw(:,1));

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data.raw(:,1), {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [57; 58; 525; 526];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
peyroche_plateau_2012.orfs = t;
peyroche_plateau_2012.ph = hit_data_names;
peyroche_plateau_2012.data = t2;
peyroche_plateau_2012.dataset_ids = hit_data_ids;

%% Save

save('./peyroche_plateau_2012.mat','peyroche_plateau_2012');

%% Print out

fid = fopen('./peyroche_plateau_2012.txt','w');
write_matrix_file(fid, peyroche_plateau_2012.orfs, peyroche_plateau_2012.ph, peyroche_plateau_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(peyroche_plateau_2012)
end

end
