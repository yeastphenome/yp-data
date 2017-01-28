%% Yu~Bellaoui, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
yu_bellaoui_2008.pmid = 19043571;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(yu_bellaoui_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data: diploid

[FILENAMES{end+1}, diploid] = read_data('xlsread','./raw_data/13_15_data.xlsx', '13&amp;15 diploid (Figure 2)');

% Get the list of ORFs and the correponding data 
dip_strains = diploid(2:end,1);
C = regexp(dip_strains, '::', 'split');
dip_strains = {};
for i = 1:length(C)
    dip_strains = [dip_strains; C{i}{1}];
end

% Get the type
dip_type = diploid(2:end, 11);

% Get the data itself
dip_data = diploid(2:end,3:4);
dip_data = cell2mat(dip_data);
   
% Eliminate all white spaces & capitalize
dip_strains = clean_orf(dip_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(dip_strains));
dip_strains(inds) = [];
dip_data(inds) = [];
dip_type(inds) = [];

% Find which are het and which are hom
data = nan(length(dip_strains), 4);

indx = find(strcmp(dip_type, 'het'));
data(indx,1) = dip_data(indx,1);
data(indx,2) = dip_data(indx,2);

indx = find(strcmp(dip_type, 'hom'));
data(indx,3) = dip_data(indx,1);
data(indx,4) = dip_data(indx,2);

% If the same strain is present more than once, average its values
[dip_strains, data] = grpstats(data, dip_strains, {'gname','mean'});

%% Load the data: haploid - 13

[FILENAMES{end+1}, haploid13] = read_data('xlsread','./raw_data/13_15_data.xlsx', '13 haploid (Figure 3)');

% Get the list of ORFs and the correponding data 
hap13_strains = haploid13(2:end,1);

% Get the data itself
hap13_data = haploid13(2:end,4);
hap13_data = cell2mat(hap13_data);
   
% Eliminate all white spaces & capitalize
hap13_strains = clean_orf(hap13_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hap13_strains));
hap13_strains(inds) = [];
hap13_data(inds) = [];

%% Load the data: haploid - 15

[FILENAMES{end+1}, haploid15] = read_data('xlsread','./raw_data/13_15_data.xlsx', '15 haploid (Figure 3)');

% Get the list of ORFs and the correponding data 
hap15_strains = haploid15(2:end,1);

% Get the data itself
hap15_data = haploid15(2:end,4);
hap15_data = cell2mat(hap15_data);
   
% Eliminate all white spaces & capitalize
hap15_strains = clean_orf(hap15_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hap15_strains));
hap15_strains(inds) = [];
hap15_data(inds) = [];

%% Combine all the data

hit_strains = unique([dip_strains; hap13_strains; hap15_strains]);
hit_data = nan(length(hit_strains), 6);

[~, ind1, ind2] = intersect(hit_strains, dip_strains);
hit_data(ind1, 1:4) = data(ind2, :);

[~, ind1, ind2] = intersect(hit_strains, hap13_strains);
hit_data(ind1, 5) = hap13_data(ind2);

[~, ind1, ind2] = intersect(hit_strains, hap15_strains);
hit_data(ind1, 6) = hap15_data(ind2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5007; 5008; 508; 4982; 511; 4983];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
yu_bellaoui_2008.orfs = hit_strains;
yu_bellaoui_2008.ph = hit_data_names;
yu_bellaoui_2008.data = hit_data;
yu_bellaoui_2008.dataset_ids = hit_data_ids;

%% Save

save('./yu_bellaoui_2008.mat','yu_bellaoui_2008');

%% Print out

fid = fopen('./yu_bellaoui_2008.txt','w');
write_matrix_file(fid, yu_bellaoui_2008.orfs, yu_bellaoui_2008.ph, yu_bellaoui_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(yu_bellaoui_2008)
end

end

