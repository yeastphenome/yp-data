%% Narayanaswamy~Marcotte, 2006
% DATA = narayanaswamy_marcotte_2006
function FILENAMES = code()
FILENAMES = {};
narayanaswamy_marcotte_2006.pmid = 16507139;

phenotypes = {'expression of PIS1'};
treatments = {'standard'};

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/gb-2006-7-1-r6-s2.xlsx');

orfs = data.raw(41:end,1);
raw_data = data.raw(41:end,2:47);
phenotypes = data.raw(40,2:47)';

inds = find(cellfun(@isempty, orfs) | cellfun(@isnumeric, orfs));
orfs(inds) = [];
raw_data(inds,:) = [];

inds = find(~strncmp('Y', orfs,1));
orfs(inds) = [];
raw_data(inds,:) = [];

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

phenotypes3 = {'slight shmoos', 'normal shmoos', 'other shmoos', 'large cells', 'small cells','round cells','pointed cells','elongated cells','pseudohyphal cells',...
    'clumpy cells','budding cells','polarized bud growth cells'};

% Remore OTHER as it's not clear what phenotype that is
phenotypes2(end) = [];
raw_data2(:,end) = [];

orfs = upper(strtrim(orfs));
[t,t2] = grpstats(raw_data2, orfs,{'mean','gname'});

narayanaswamy_marcotte_2006.orfs = t2;
narayanaswamy_marcotte_2006.data = t;
narayanaswamy_marcotte_2006.ph = strcat(phenotypes3, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'narayanaswamy_marcotte_2006.mat'],'narayanaswamy_marcotte_2006');
return;

% Save data into database
dt = narayanaswamy_marcotte_2006;

% datasets = get_datasets_for_paper(dt);
%
% [~,database_ix] = sortrows(datasets.names);
% [~,ph_ix] = sort(dt.ph);
ph_ix = 1:length(dt.ph);
%
% % Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
% datasets.names(database_ix,:)
% dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, [548 549 592:601]);


end

