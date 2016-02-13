%% Narayanaswamy~Marcotte, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
narayanaswamy_marcotte_2006.pmid = 16507139;

treatments = {'standard'};

% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/gb-2006-7-1-r6-s2.xlsx');

orfs = data.raw(41:end,1);
raw_data = data.raw(41:end,2:47);
phenotypes = data.raw(40,2:47)';

orfs = clean_orf(orfs);

inds = find(~is_orf(orfs));
disp(orfs(inds));
orfs(strcmp('YLR287-A', orfs)) = {'YLR287C-A'};

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
phmap.ph_old = tmp(2:end,1);
phmap.ph_new = tmp(1,2:end)';
phmap.mat = cell2mat(tmp(2:end,2:end));

raw_data3 = zeros(size(raw_data2,1), length(phmap.ph_new));
for i = 1 : length(phenotypes2)
    ind = find(strcmp(phenotypes2{i}, phmap.ph_old));
    ind2 = find(abs(phmap.mat(ind,:))>0);
    raw_data3(:,ind2) = raw_data3(:,ind2) + raw_data2(:,ind)*phmap.mat(ind,ind2);
end

% Note: some of the phenotypes cancel each other out (e.g., small and large cells) because they were both annotated to the same ORFs.

[t,t2] = grpstats(raw_data3, orfs,{'mean','gname'});

narayanaswamy_marcotte_2006.orfs = t2;
narayanaswamy_marcotte_2006.data = t;
narayanaswamy_marcotte_2006.ph = strcat(phmap.ph_new', '; ', treatments);
narayanaswamy_marcotte_2006.ph = narayanaswamy_marcotte_2006.ph';

save('./narayanaswamy_marcotte_2006.mat','narayanaswamy_marcotte_2006');

end

