%% Ooi~Boeke, 2001
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ooi_boeke_2001.pmid = 11701889;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ooi_boeke_2001.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/ooi_boeke_2001.xlsx');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);

% Get the data itself
hit_data = zeros(length(hit_strains), 2);
indx = find(~cellfun(@isempty, strfind(data(:,3), '+')));
hit_data(indx, 1) = -1;
indx = find(~cellfun(@isempty, strfind(data(:,2), '+')));
hit_data(indx, 2) = -1;
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
hit_strains(ismember(hit_strains, {'GPE2'})) = {'YAL056W'};
inds = find(~is_orf(hit_strains));

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [2; 315];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
ooi_boeke_2001.orfs = hit_strains;
ooi_boeke_2001.ph = hit_data_names;
ooi_boeke_2001.data = hit_data;
ooi_boeke_2001.dataset_ids = hit_data_ids;

%% Save

save('./ooi_boeke_2001.mat','ooi_boeke_2001');

%% Print out

fid = fopen('./ooi_boeke_2001.txt','w');
write_matrix_file(fid, ooi_boeke_2001.orfs, ooi_boeke_2001.ph, ooi_boeke_2001.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ooi_boeke_2001)
end

end

