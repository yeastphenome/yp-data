%% Parsons~Boone, 2004
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
% TESTED = not available

parsons_boone_2004.pmid = 14661025;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(parsons_boone_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/nbt919-S2.xlsx', 'Sheet1');

inds = find(strcmp('ORF', data.raw(:,1)));

hit_strains = data.raw(inds+1:end,1);

hit_data = data.raw(inds+1:end,3:end);
hit_data(cellfun(@isnan, hit_data)) = {0};
hit_data = cell2mat(hit_data);

treatments = data.raw(inds,3:end)';

% Flip the sign of the values, such that negative = sensitive, positive =
% resistant
hit_data = -hit_data;

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% Load the mapping of treatments to dataset ids
[FILENAMES{end+1}, T] = read_data('textread','./extras/datasets.txt', '%d %s','delimiter','\t');
treatments_ids = T{1};
treatments_names = T{2};

hit_data_ids = nan(size(hit_data,2),1);
[~,ind1,ind2] = intersect(treatments, treatments_names);
hit_data_ids(ind1) = treatments_ids(ind2);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
parsons_boone_2004.orfs = hit_strains;
parsons_boone_2004.ph = hit_data_names;
parsons_boone_2004.data = hit_data;
parsons_boone_2004.dataset_ids = hit_data_ids;

%% Save

save('./parsons_boone_2004.mat','parsons_boone_2004');

%% Print out

fid = fopen('./parsons_boone_2004.txt','w');
write_matrix_file(fid, parsons_boone_2004.orfs, parsons_boone_2004.ph, parsons_boone_2004.data);
fclose(fid);

end

