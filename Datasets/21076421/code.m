%% Baryshnikova~Myers, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
baryshnikova_myers_2010.pmid = 21076421;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(baryshnikova_myers_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Supplementary_data_1_SMF_standard_100209.xlsx');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hits_orfs = data.raw(:,1);
hits_data = data.raw(:,2);

% Eliminate the data for non-deletion strains
t = regexp(hits_orfs,'_','split');
inds = find(cellfun(@length, t) > 1);
hits_orfs(inds) = [];
hits_data(inds) = [];

% Eliminate all white spaces & capitalize
hits_orfs = clean_orf(hits_orfs);

% If in gene name form, transform into ORF name
hits_orfs = translate(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

% Find and remove any non-numeric data values
hits_data = cell2mat(hits_data);

% If the same strain is present more than once, average its values
[t,t2] = grpstats(hits_data, hits_orfs,{'mean','gname'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [540];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

baryshnikova_myers_2010.orfs = t2;
baryshnikova_myers_2010.ph = hit_data_names;
baryshnikova_myers_2010.data = t;
baryshnikova_myers_2010.dataset_ids = hit_data_ids;

%% Save
save('./baryshnikova_myers_2010.mat','baryshnikova_myers_2010');

%% Print out
fid = fopen('./baryshnikova_myers_2010.txt','w');
write_matrix_file(fid, baryshnikova_myers_2010.orfs, baryshnikova_myers_2010.ph, baryshnikova_myers_2010.data);
fclose(fid);

%% Save to DB (admin)
addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(baryshnikova_myers_2010)
end

end
