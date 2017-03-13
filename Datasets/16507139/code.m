%% Narayanaswamy~Marcotte, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
narayanaswamy_marcotte_2006.pmid = 16507139;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(narayanaswamy_marcotte_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};


%% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/gb-2006-7-1-r6-s2.xlsx');

orfs = data.raw(41:end,1);
raw_data = data.raw(41:end,2:47);
phenotypes = data.raw(40,2:47)';

orfs = clean_orf(orfs);

orfs(strcmp('YLR287-A', orfs)) = {'YLR287C-A'};

inds = find(~is_orf(orfs));
disp(orfs(inds));

raw_data = cell2mat(raw_data);

% Multiply intensity by penetrance
int_inds = find(~cellfun(@isempty, regexp(phenotypes,'_INT')));
pen_inds = find(~cellfun(@isempty, regexp(phenotypes,'_PEN')));

raw_data(:,int_inds) = raw_data(:,int_inds) .* raw_data(:,pen_inds)/4;

% Average between observers
obs1_inds = find(~cellfun(@isempty, regexp(phenotypes,'22_')));
obs2_inds = find(~cellfun(@isempty, regexp(phenotypes,'24_')));

raw_data(:,obs1_inds) = (raw_data(:,obs1_inds) + raw_data(:,obs2_inds))/2;

shmoo_inds = find(~cellfun(@isempty, regexp(phenotypes,'SHMOOS')));

final_inds = intersect([int_inds; shmoo_inds], obs1_inds);
raw_data2 = raw_data(:,final_inds);
phenotypes2 = phenotypes(final_inds);

% Remove OTHER and NORMAL as it's not clear what phenotype that is
inds_other = find(~cellfun(@isempty, regexp(phenotypes2, 'OTHER')));
inds_normal = find(~cellfun(@isempty, regexp(phenotypes2, 'NORMAL')));
phenotypes2([inds_other; inds_normal]) = [];
raw_data2(:,[inds_other; inds_normal]) = [];

% Map the phenotypes to a common standard
[FILENAMES{end+1}, tmp] = read_data('xlsread', './raw_data/phenotype_mapping.xlsx','Sheet1');
phmap.ph_old = tmp(3:end,1);
phmap.ph_new = tmp(1,2:end)';
hit_data_ids = cell2mat(tmp(2,2:end)');
phmap.mat = cell2mat(tmp(3:end,2:end));

raw_data3 = zeros(size(raw_data2,1), length(phmap.ph_new));
for i = 1 : length(phenotypes2)
    ind = find(strcmp(phenotypes2{i}, phmap.ph_old));
    ind2 = find(abs(phmap.mat(ind,:))>0);
    raw_data3(:,ind2) = raw_data3(:,ind2) + raw_data2(:,ind)*phmap.mat(ind,ind2);
end

% Note: some of the phenotypes cancel each other out (e.g., small and large cells) because they were both annotated to the same ORFs.

[hit_data,hit_orfs] = grpstats(raw_data3, orfs,{'mean','gname'});


%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
narayanaswamy_marcotte_2006.orfs = hit_orfs;
narayanaswamy_marcotte_2006.ph = hit_data_names;
narayanaswamy_marcotte_2006.data = hit_data;
narayanaswamy_marcotte_2006.dataset_ids = hit_data_ids;

%% Save

save('./narayanaswamy_marcotte_2006.mat','narayanaswamy_marcotte_2006');

%% Print out

fid = fopen('./narayanaswamy_marcotte_2006.txt','w');
write_matrix_file(fid, narayanaswamy_marcotte_2006.orfs, narayanaswamy_marcotte_2006.ph, narayanaswamy_marcotte_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(narayanaswamy_marcotte_2006)
end

end

