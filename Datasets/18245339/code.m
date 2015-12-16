%% Abe~Minegishi, 2008
function FILENAMES = code()
FILENAMES = {};
abe_minegishi_2008.pmid = 18245339;

phenotypes = {'growth (culture turbidity)'};
treatments = {'high pressure';'low temperature'};

% Load tested strains
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/mat_alpha_041902.xlsx', 'mat_alpha_041902.txt');

tested_orfs = tested.raw(4:end,3);

inds = find(cellfun(@isempty, tested_orfs));
tested_orfs(inds) = [];

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(upper(cellfun(@strtrim, tested_orfs,'UniformOutput',0)));

% Eliminate all spaces from within the ORF name
tested_orfs = regexprep(tested_orfs, '\W','');

inds = find(cellfun(@isempty, regexp(tested_orfs,'Y[A-P][RL][0-9]{3}[CW](-[ABC])*')));
tested_orfs(inds) = [];

<<<<<<< Updated upstream
% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/Table1_Abe_Genetics.xlsx', 'Table 1 (2)');
=======
[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/Table1_Abe_Genetics.xlsx', 'Table 1 (2)');
>>>>>>> Stashed changes

hits_orfs = data.raw(7:end,3);
hits_data = data.raw(7:end,[26 31]);

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_orfs = upper(regexprep(hits_orfs, '\W',''));

inds = find(cellfun(@isempty, regexp(hits_orfs,'Y[A-P][RL][0-9]{3}[CW](-[ABC])*')));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_data = cell2mat(hits_data);
hits_data = hits_data/100;  % transform percent into fractions

[missing, ix] = setdiff(hits_orfs, tested_orfs);    % nothing missing

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(hits_data, hits_orfs, {'gname','mean'});

abe_minegishi_2008.orfs = tested_orfs;
abe_minegishi_2008.data = zeros(length(tested_orfs),2);

[~,ind1,ind2] = intersect(t, tested_orfs);
abe_minegishi_2008.data(ind2,:) = t2(ind1,:);

abe_minegishi_2008.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'abe_minegishi_2008.mat'],'abe_minegishi_2008');
return;

% Save data into database
dt = abe_minegishi_2008;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

