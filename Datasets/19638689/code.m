%% Auesukaree~Harashima, 2009
% DATA = auesukaree_harashima_2009
function FILENAMES = code()
FILENAMES = {};
auesukaree_harashima_2009.pmid = 19638689;

phenotypes = {'growth'};
treatments = {'EtOH','MeOH','propanol','NaCl','H2O2','37C'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','raw_data/Mat alpha_KOset list.xlsx');
tested_orfs = tested.raw(4:end,2);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
[FILENAMES{end+1}, data_hits{1}] = dataread('textread','raw_data/ethanol_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{2}] = dataread('textread','raw_data/methanol_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{3}] = dataread('textread','raw_data/propanol_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{4}] = dataread('textread','raw_data/nacl_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{5}] = dataread('textread','raw_data/h2o2_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{6}] = dataread('textread','raw_data/heat_sensitivity_hits.txt', '%s');

data_hits_orfs = cell(size(data_hits));

data_hits_orfs{1} = genename2orf(data_hits{1});
data_hits_orfs{1}(strcmp('rlr1', data_hits_orfs{1})) = {'YNL139C'};
data_hits_orfs{1}(strcmp('sur4', data_hits_orfs{1})) = {'YLR372W'};
data_hits_orfs{1} = unique(upper(data_hits_orfs{1}));

[missing, ix] = setdiff(data_hits_orfs{1}, tested_orfs);
data_hits_orfs{1}(strcmp('YHR039C-A', data_hits_orfs{1})) = {'YHR039C-B'};


data_hits_orfs{2} = genename2orf(data_hits{2});
data_hits_orfs{2}(strcmp('fen1', data_hits_orfs{2})) = [];  % ambiguous gene name
data_hits_orfs{2}(strcmp('fmp13', data_hits_orfs{2})) = {'YKR016W'};
data_hits_orfs{2} = unique(upper(data_hits_orfs{2}));

[missing, ix] = setdiff(data_hits_orfs{2}, tested_orfs);


data_hits_orfs{3} = genename2orf(data_hits{3});
data_hits_orfs{3}(strcmp('caf17', data_hits_orfs{3})) = {'YJR122W'};
data_hits_orfs{3}(strcmp('fmp13', data_hits_orfs{3})) = {'YKR016W'};
data_hits_orfs{3}(strcmp('kem1', data_hits_orfs{3})) = {'YGL173C'};
data_hits_orfs{3}(strcmp('ppa1', data_hits_orfs{3})) = [];  % ambiguous gene name
data_hits_orfs{3}(strcmp('tfp1', data_hits_orfs{3})) = {'YDL185W'};
data_hits_orfs{3} = unique(upper(data_hits_orfs{3}));

[missing, ix] = setdiff(data_hits_orfs{3}, tested_orfs);


data_hits_orfs{4} = genename2orf(data_hits{4});
data_hits_orfs{4}(strcmp('tfp1', data_hits_orfs{4})) = {'YDL185W'};
data_hits_orfs{4} = unique(upper(data_hits_orfs{4}));

[missing, ix] = setdiff(data_hits_orfs{4}, tested_orfs);


data_hits_orfs{5} = genename2orf(data_hits{5});
data_hits_orfs{5}(strcmp('cup5', data_hits_orfs{5})) = {'YEL027W'};
data_hits_orfs{5} = unique(upper(data_hits_orfs{5}));

[missing, ix] = setdiff(data_hits_orfs{5}, tested_orfs);


data_hits_orfs{6} = genename2orf(data_hits{6});
data_hits_orfs{6}(strcmp('cup5', data_hits_orfs{6})) = {'YEL027W'};
data_hits_orfs{6}(strcmp('ppa1', data_hits_orfs{6})) = [];  % ambiguous gene name
data_hits_orfs{6}(strcmp('sur4', data_hits_orfs{6})) = {'YLR372W'};
data_hits_orfs{6}(strcmp('vid21', data_hits_orfs{6})) = {'YDR359C'};
data_hits_orfs{6}(strcmp('vps66', data_hits_orfs{6})) = {'YPR139C'};
data_hits_orfs{6}(strcmp('vps69', data_hits_orfs{6})) = {'YPR087W'};
data_hits_orfs{6} = unique(upper(data_hits_orfs{6}));

[missing, ix] = setdiff(data_hits_orfs{6}, tested_orfs);


auesukaree_harashima_2009.orfs = tested_orfs;
auesukaree_harashima_2009.data = zeros(length(tested_orfs),length(treatments));

for i = 1 : length(treatments)
    [~,ind1,ind2] = intersect(data_hits_orfs{i}, tested_orfs);
    auesukaree_harashima_2009.data(ind2,i) = -1;
end

auesukaree_harashima_2009.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'auesukaree_harashima_2009.mat'],'auesukaree_harashima_2009');
return;

% Save data into database
dt = auesukaree_harashima_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
ph_ix = ph_ix([2 3 4 5 6 1]);
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

