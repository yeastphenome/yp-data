%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
eide_harper_2005.pmid = 16168084;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(eide_harper_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/13059_2005_1112_MOESM2_ESM.xlsx', 'Normalized concentrations');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);
hit_screens = data(1,2:14)';

% Get the data itself
hit_data = cell2mat(data(2:end,2:14));
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YILO16W'})) = {'YIL016W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
[FILENAMES{end+1}, data] = read_data('xlsread','./extras/datasets.xlsx', 'Sheet1');
[~,inds] = ismember(hit_screens, data(:,1));
hit_data_ids = cell2mat(data(inds,2));

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('xlsread','./raw_data/13059_2005_1112_MOESM1_ESM.xlsx', 'Genes Analyzed');
tested_strains = tested_strains(2:end,1);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'YELOO1C'})) = {'YEL001C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

tested_strains(inds) = [];

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
eide_harper_2005.orfs = tested_strains;
eide_harper_2005.ph = hit_data_names;
eide_harper_2005.data = zeros(length(eide_harper_2005.orfs),length(eide_harper_2005.ph));
eide_harper_2005.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, eide_harper_2005.orfs);
eide_harper_2005.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./eide_harper_2005.mat','eide_harper_2005');

%% Print out

fid = fopen('./eide_harper_2005.txt','w');
write_matrix_file(fid, eide_harper_2005.orfs, eide_harper_2005.ph, eide_harper_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(eide_harper_2005)
end

end

