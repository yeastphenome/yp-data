%% Lee~Giaever, 2014
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
lee_giaever_2014.pmid = 24723613;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(lee_giaever_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load text file
[FILENAMES{end+1}, data] = read_data('read_matrix_file', './raw_data/fitness_defect_matrix_hom.txt',1,1);
data.labels_col = strtrim(data.labels_col);

% Get compound names
[FILENAMES{end+1}, compounds] = read_data('readtable','./raw_data/1250217s1.xlsx', 'Sheet','compound library');
mu = native2unicode([206 188],'UTF-8');

inds = find(strcmp([mu 'M'], compounds.unit));
compounds.unit(inds) = {'uM'};

inds = find(strcmp('-', compounds.PCID));
compounds.PCID(inds) = {'NaN'};

compounds.PCID = cellfun(@str2num, compounds.PCID);

% Match screenIDs with dataset ids
[FILENAMES{end+1}, screen2dataset] = read_data('readtable', './extras/screenid_datasetid.txt','delimiter','\t','ReadVariableNames', false);
[~, ind1, ind2] = intersect(compounds.ScreenID, screen2dataset.Var1);
compounds.datasetid = nan(size(compounds.ScreenID));
compounds.datasetid(ind1) = screen2dataset.Var2(ind2);

% Match data with dataset ids
[~,ind1,ind2] = intersect(compounds.ScreenID, data.labels_col);
data.datasetid = nan(size(data.labels_col));
data.datasetid(ind2) = compounds.datasetid(ind1);

% Clean up the ORFs
data.labels_row = clean_orf(data.labels_row);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(data.labels_row));
disp(data.labels_row(inds));

% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
[data.datasetid, data.data] = grpstats(data.data', data.datasetid, {'gname','mean'});
data.data = data.data';
data.datasetid = cellfun(@str2num, data.datasetid);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, data.datasetid);
data.dataname = cell(size(data.datasetid));
data.dataname(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
lee_giaever_2014.orfs = data.labels_row;
lee_giaever_2014.ph = data.dataname;
lee_giaever_2014.data = data.data;
lee_giaever_2014.dataset_ids = data.datasetid;


%% Save
save('./lee_giaever_2014.mat','lee_giaever_2014');

fid = fopen('./lee_giaever_2014.txt','w');
write_matrix_file(fid, lee_giaever_2014.orfs, lee_giaever_2014.ph, lee_giaever_2014.data);
fclose(fid);
end
