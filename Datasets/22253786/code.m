%% Blackman~Nislow, 2012
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
blackman_nislow_2012.pmid = 22253786;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(blackman_nislow_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/pone.0029798.s003.xlsx', 'elescomol');

% Get the list of ORFs
strains = data(2:end, 1);

% Clean up ORFs
strains = cellfun(@(x) strtok(x, ':'), strains, 'UniformOutput', false); 
strains = clean_orf(strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(strains));
strains(inds) = [];

%% Get data from hits
% Find the rows with het
ind_het = find(~cellfun(@isnumeric, data(2:end, 15)));
% Make a data matrix
hit_data = zeros(length(strains), 2);
hit_data(:, 1) = cell2mat(data(2:end, 3)); 
hit_data(ind_het, 2) = hit_data(ind_het, 1);
hit_data(ind_het, 1) = NaN;

% Average any repeated value
[strains, hit_data] = grpstats(hit_data, strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [86; 87];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

blackman_nislow_2012.orfs = strains;
blackman_nislow_2012.data = hit_data;
blackman_nislow_2012.ph = hit_data_names;
blackman_nislow_2012.dataset_ids = hit_data_ids;

%% Save

save('./blackman_nislow_2012.mat','blackman_nislow_2012');

%% Print out

fid = fopen('./blackman_nislow_2012.txt','w');
write_matrix_file(fid, blackman_nislow_2012.orfs, blackman_nislow_2012.ph, blackman_nislow_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(blackman_nislow_2012)
end

end
