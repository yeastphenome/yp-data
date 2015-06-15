%% GENERAL COMMENTS
%
% - Sometimes the downloaded XLS files had to be resaved as XLSX to facilitate import into Matlab.
% - The various cells in this file, relative to the various datasets, are independent from each other and can be run individually.
% - If the same ORF appears in a dataset more than once, its values are averaged.
% - If an ORF appears in the hit list but not in the set of tested ORFs, it
% gets added to the tested list.
% - If there is overlap between positive and negative hits in a given
% screen (e.g., sensitive and resistant mutants), they get removed from
% both lists for inconsistency.
% - "NaN" always means "data not available". E.g., 1) the mutant wasn't tested; 2) the mutant was tested
% but its phenotype is unknown for technical reasons.
% - 0 = absence of phenotype (mutant behaves like expected/wild-type)

%%

addpath(genpath('Datasets/'))
addpath(genpath('Utils/Matlab/'))


% PART II


%% Auesukaree~Harashima, 2009
% DATA = auesukaree_harashima_2009
auesukaree_harashima_2009.pmid = 19638689;

phenotypes = {'growth'};
treatments = {'EtOH','MeOH','propanol','NaCl','H2O2','37C'};

% Load tested
[tested.txt, tested.num, tested.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Auesukaree~Harashima/Mat alpha_KOset list.xlsx');
tested_orfs = tested.raw(4:end,2);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
data_hits{1} = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Auesukaree~Harashima/ethanol_sensitivity_hits.txt','%s');
data_hits{2} = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Auesukaree~Harashima/methanol_sensitivity_hits.txt','%s');
data_hits{3} = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Auesukaree~Harashima/propanol_sensitivity_hits.txt','%s');
data_hits{4} = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Auesukaree~Harashima/nacl_sensitivity_hits.txt','%s');
data_hits{5} = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Auesukaree~Harashima/h2o2_sensitivity_hits.txt','%s');
data_hits{6} = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Auesukaree~Harashima/heat_sensitivity_hits.txt','%s');

data_hits_orfs = cell(size(data_hits));

data_hits_orfs{1} = genename2orf_sgd(data_hits{1});
data_hits_orfs{1}(strcmp('rlr1', data_hits_orfs{1})) = {'YNL139C'};
data_hits_orfs{1}(strcmp('sur4', data_hits_orfs{1})) = {'YLR372W'};
data_hits_orfs{1} = unique(upper(data_hits_orfs{1}));

[missing, ix] = setdiff(data_hits_orfs{1}, tested_orfs);
data_hits_orfs{1}(strcmp('YHR039C-A', data_hits_orfs{1})) = {'YHR039C-B'};


data_hits_orfs{2} = genename2orf_sgd(data_hits{2});
data_hits_orfs{2}(strcmp('fen1', data_hits_orfs{2})) = [];  % ambiguous gene name
data_hits_orfs{2}(strcmp('fmp13', data_hits_orfs{2})) = {'YKR016W'};
data_hits_orfs{2} = unique(upper(data_hits_orfs{2}));

[missing, ix] = setdiff(data_hits_orfs{2}, tested_orfs);


data_hits_orfs{3} = genename2orf_sgd(data_hits{3});
data_hits_orfs{3}(strcmp('caf17', data_hits_orfs{3})) = {'YJR122W'};
data_hits_orfs{3}(strcmp('fmp13', data_hits_orfs{3})) = {'YKR016W'};
data_hits_orfs{3}(strcmp('kem1', data_hits_orfs{3})) = {'YGL173C'};
data_hits_orfs{3}(strcmp('ppa1', data_hits_orfs{3})) = [];  % ambiguous gene name
data_hits_orfs{3}(strcmp('tfp1', data_hits_orfs{3})) = {'YDL185W'};
data_hits_orfs{3} = unique(upper(data_hits_orfs{3}));

[missing, ix] = setdiff(data_hits_orfs{3}, tested_orfs);


data_hits_orfs{4} = genename2orf_sgd(data_hits{4});
data_hits_orfs{4}(strcmp('tfp1', data_hits_orfs{4})) = {'YDL185W'};
data_hits_orfs{4} = unique(upper(data_hits_orfs{4}));

[missing, ix] = setdiff(data_hits_orfs{4}, tested_orfs);


data_hits_orfs{5} = genename2orf_sgd(data_hits{5});
data_hits_orfs{5}(strcmp('cup5', data_hits_orfs{5})) = {'YEL027W'};
data_hits_orfs{5} = unique(upper(data_hits_orfs{5}));

[missing, ix] = setdiff(data_hits_orfs{5}, tested_orfs);


data_hits_orfs{6} = genename2orf_sgd(data_hits{6});
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

save('Datasets/Phenotypes/2009_Auesukaree~Harashima/auesukaree_harashima_2009.mat','auesukaree_harashima_2009'); 

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

%% Aouida~Ramotar, 2004
% DATA = aouida_ramotar_2004
aouida_ramotar_2004.pmid = 14871844;

phenotypes = {'growth'};
treatments = {'bleomycin'};

% Load tested
[tested.txt, tested.num, tested.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2004_Aouida~Ramotar/HU haploid.xlsx');
tested_orfs = tested.raw(4:end,5);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
[data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2004_Aouida~Ramotar/CAN_2-1-04_Aouida.xlsx');
hits_orfs = data.raw(:,1);
hits_orfs = strtrim(upper(hits_orfs));



hits_data = strtrim(data.raw(:,2));

% Data conversion
hits_data(strcmp('R+++', hits_data)) = {9};     % (>500-fold more resistant -> log2(500) ~ 9)
hits_data(strcmp('R++++', hits_data)) = {10};
hits_data(strcmp('S+', hits_data)) = {-1};
hits_data(strcmp('S++', hits_data)) = {-2};
hits_data(strcmp('S+++', hits_data)) = {-3};
hits_data(strcmp('S++++', hits_data)) = {-4};

hits_data = cell2mat(hits_data);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];


aouida_ramotar_2004.orfs = tested_orfs;
aouida_ramotar_2004.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
aouida_ramotar_2004.data(ind2) = hits_data(ind1);

aouida_ramotar_2004.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2004_Aouida~Ramotar/aouida_ramotar_2004.mat','aouida_ramotar_2004'); 

% Save data into database
dt = aouida_ramotar_2004;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Mira~Sa-Correia, 2009
% DATA = mira_sa_correia_2009
mira_sa_correia_2009.pmid = 19220866;

phenotypes = {'growth'};
treatments = {'propionic acid'};

% Load tested
[tested.txt, tested.num, tested.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Mira~Sa-Correia/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
hits_genenames = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Mira~Sa-Correia/hits_genenames.txt','%s');
hits_genenames_moderate = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Mira~Sa-Correia/hits_genenames_moderate.txt','%s');

hits_genenames = strtrim(lower(hits_genenames));
hits_genenames_moderate = strtrim(lower(hits_genenames_moderate));

hits_genenames_strong = setdiff(hits_genenames, hits_genenames_moderate);


hits_orfs_moderate = genename2orf_sgd(hits_genenames_moderate);
hits_orfs_moderate(strcmp('set7', hits_orfs_moderate)) = {'YDR257C'};
hits_orfs_moderate = unique(upper(hits_orfs_moderate));

hits_orfs_strong = genename2orf_sgd(hits_genenames_strong);
hits_orfs_strong(strcmp('mrp16', hits_orfs_strong)) = [];    % typo? genename doesn't exist
hits_orfs_strong(strcmp('rhr2', hits_orfs_strong)) = {'YIL053W'};
hits_orfs_strong(strcmp('tfp1', hits_orfs_strong)) = {'YDL185W'};
hits_orfs_strong(strcmp('tfp3', hits_orfs_strong)) = {'YPL234C'};
hits_orfs_strong = unique(upper(hits_orfs_strong));

inds = find(~strncmp('Y', hits_orfs_moderate,1));
hits_orfs_moderate(inds) = [];

inds = find(~strncmp('Y', hits_orfs_strong,1));
hits_orfs_strong(inds) = [];

[missing, ix] = setdiff(hits_orfs_strong, tested_orfs);
hits_orfs_strong(ix) = [];      % 3 ORFs eliminated from the hit list

[missing, ix] = setdiff(hits_orfs_moderate, tested_orfs);

hits_data_strong = zeros(length(hits_orfs_strong),1)-2;
hits_data_moderate = zeros(length(hits_orfs_moderate),1)-1;


mira_sa_correia_2009.orfs = tested_orfs;
mira_sa_correia_2009.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs_strong, tested_orfs);
mira_sa_correia_2009.data(ind2) = hits_data_strong(ind1);

[~,ind1,ind2] = intersect(hits_orfs_moderate, tested_orfs);
mira_sa_correia_2009.data(ind2) = hits_data_moderate(ind1);


mira_sa_correia_2009.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2009_Mira~Sa-Correia/mira_sa_correia_2009.mat','mira_sa_correia_2009'); 

% Save data into database
dt = mira_sa_correia_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Teixeira~Sa-Correia, 2009
% DATA = teixeira_sa_correia_2009
teixeira_sa_correia_2009.pmid = 19633105;

phenotypes = {'growth'};
treatments = {'EtOH'};

% Load tested
[tested.txt, tested.num, tested.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Teixeira~Sa-Correia/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
[data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Teixeira~Sa-Correia/TableS1_suplementary_material.xlsx');
hits_genenames = data.raw(7:end,1);
inds = find(cellfun(@isempty, hits_genenames) | cellfun(@isnumeric, hits_genenames));
hits_genenames(inds) = [];

hits_orfs = genename2orf_sgd(hits_genenames);
hits_orfs(strcmp('opi8', hits_orfs))= {'YKR035C'};
hits_orfs(strcmp('opi9', hits_orfs)) = {'YLR338W'};
hits_orfs(strcmp('ppa1', hits_orfs)) = [];  % ambiguous genename
hits_orfs(strcmp('rlr1', hits_orfs)) = {'YNL139C'};
hits_orfs(strcmp('soy1', hits_orfs)) = {'YBR194W'};
hits_orfs(strcmp('tfp3', hits_orfs)) = {'YPL234C'};
hits_orfs(strcmp('vps66', hits_orfs)) = {'YPR139C'};

hits_orfs = upper(hits_orfs);
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];


[missing, ix] = setdiff(hits_orfs, tested_orfs);
hits_orfs(ix) = [];      % 24 ORFs eliminated from the hit list (the list contains hits from both the hap and the het collection, but the list of tested strains is only available for the hap collection)

teixeira_sa_correia_2009.orfs = tested_orfs;
teixeira_sa_correia_2009.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
teixeira_sa_correia_2009.data(ind2) = -1;

teixeira_sa_correia_2009.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2009_Teixeira~Sa-Correia/teixeira_sa_correia_2009.mat','teixeira_sa_correia_2009'); 

% Save data into database
dt = teixeira_sa_correia_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
database_ix = database_ix(2);
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Teixeira~Sa-Correia, 2010
% DATA = teixeira_sa_correia_2010
teixeira_sa_correia_2010.pmid = 20210661;

phenotypes = {'growth'};
treatments = {'Glu, 30%'};

% Load tested
[tested.txt, tested.num, tested.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2010_Teixeira~Sa-Correia/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
hits_genenames = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2010_Teixeira~Sa-Correia/hits_genenames.txt','%s');
inds = find(cellfun(@isempty, hits_genenames) | cellfun(@isnumeric, hits_genenames));
hits_genenames(inds) = [];
hits_genenames = strtrim(hits_genenames);

hits_orfs = genename2orf_sgd(hits_genenames);
hits_orfs(strcmp('opi9', hits_orfs)) = {'YLR338W'};
hits_orfs(strcmp('ppa1', hits_orfs)) = [];  % ambiguous genename
hits_orfs(strcmp('soy1', hits_orfs)) = {'YBR194W'};

hits_orfs = upper(hits_orfs);
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];

hits_orfs = unique(hits_orfs);


[missing, ix] = setdiff(hits_orfs, tested_orfs);
hits_orfs(strcmp('YML009W-B',hits_orfs)) = {'YML010W-A'};

teixeira_sa_correia_2010.orfs = tested_orfs;
teixeira_sa_correia_2010.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
teixeira_sa_correia_2010.data(ind2) = -1;

teixeira_sa_correia_2010.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2010_Teixeira~Sa-Correia/teixeira_sa_correia_2010.mat','teixeira_sa_correia_2010'); 

% Save data into database
dt = teixeira_sa_correia_2010;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Mira~Sa-Correia, 2010
% DATA = mira_sa_correia_2010
mira_sa_correia_2010.pmid = 20973990;

phenotypes = {'growth'};
treatments = {'acetic acid'};

% Load tested
[tested.txt, tested.num, tested.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2010_Mira~Sa-Correia/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
[data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2010_Mira~Sa-Correia/1475-2859-9-79-s1.xlsx');
hits_genenames = data.raw(9:end,1);
hits_data = data.raw(9:end,3);

inds = find(cellfun(@isempty, hits_genenames) | cellfun(@isnumeric, hits_genenames));
hits_genenames(inds) = [];
hits_data(inds,:) = [];

hits_genenames = strtrim(hits_genenames);

hits_orfs = genename2orf_sgd(hits_genenames);
hits_orfs(strcmp('ace1', hits_orfs)) = {'YGL166W'};
hits_orfs(strcmp('api2', hits_orfs)) = {'YDR525W'};
hits_orfs(strcmp('brp1', hits_orfs)) = {'YGL007W'};
hits_orfs(strcmp('bud19', hits_orfs)) = {'YJL188C'};
hits_orfs(strcmp('bud26', hits_orfs)) = {'YDR241W'};
hits_orfs(strcmp('cos16', hits_orfs)) = {'YCR044C'};
hits_orfs(strcmp('cup5', hits_orfs)) = {'YEL027W'};
hits_orfs(strcmp('kem1', hits_orfs)) = {'YGL173C'};
hits_orfs(strcmp('opi8', hits_orfs)) = {'YKR035C'};
hits_orfs(strcmp('opi9', hits_orfs)) = {'YLR338W'};
hits_orfs(strcmp('rhr2', hits_orfs)) = {'YIL053W'};
hits_orfs(strcmp('see1', hits_orfs)) = {'YIL064W'};
hits_orfs(strcmp('sur4', hits_orfs)) = {'YLR372W'};
hits_orfs(strcmp('tfp1', hits_orfs)) = {'YDL185W'};
hits_orfs(strcmp('vps61', hits_orfs)) = {'YDR136C'};
hits_orfs(strcmp('vps66', hits_orfs)) = {'YPR139C'};

inds = find(strcmp('ppa1', hits_orfs));  % ambiguous genename
hits_orfs(inds) = [];
hits_data(inds,:) = [];


hits_orfs = upper(hits_orfs);
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_data(strcmp('++', hits_data)) = {-2};
hits_data(strcmp('+', hits_data)) = {-1};

hits_data = cell2mat(hits_data);
hits_data(isnan(hits_data)) = 0;

[missing, ix] = setdiff(hits_orfs, tested_orfs);    % 8 ORFs deleted
hits_orfs(ix) = [];
hits_data(ix,:) = [];


[t,t2] = grpstats(hits_data, hits_orfs, {'mean','gname'});
hits_orfs = t2;
hits_data = t;

mira_sa_correia_2010.orfs = tested_orfs;
mira_sa_correia_2010.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
mira_sa_correia_2010.data(ind2) = hits_data(ind1);

mira_sa_correia_2010.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2010_Mira~Sa-Correia/mira_sa_correia_2010.mat','mira_sa_correia_2010'); 

% Save data into database
dt = mira_sa_correia_2010;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Uluisik~Koc, 2011
% DATA = uluisik_koc_2011
uluisik_koc_2011.pmid = 21035538;

phenotypes = {'growth'};
treatments = {'BA, 20 mM','BA, 30 mM','BA, 40 mM','BA, 50 mM', 'BA, 60 mM', 'BA, 70 mM','BA, 100 mM','BA, 125 mM','BA, 150 mM'};

% Load tested
tested_orfs = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Uluisik~Koc/tested_strains.txt','%s');

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

tested_orfs = unique(upper(strtrim(tested_orfs)));

uluisik_koc_2011.orfs = tested_orfs;
uluisik_koc_2011.data = zeros(length(tested_orfs),length(treatments));
uluisik_koc_2011.ph = strcat(phenotypes, '; ', treatments);

% Load data 1
[data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Uluisik~Koc/mmc2.xlsx');
hits_orfs = data.raw(4:end,1);
hits_data = data.raw(4:end,2:5);
hits_treatments = {'BA, 100 mM','BA, 125 mM','BA, 150 mM'};

hits_data(strcmp('+',hits_data)) = {1};
hits_data(strcmp('-',hits_data)) = {0};

hits_data = cell2mat(hits_data);

inds = find(cellfun(@isempty, hits_orfs) | cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_data = hits_data - repmat(hits_data(1,:),size(hits_data,1),1);    % Subtract the 1st row (WT)
hits_data = hits_data - repmat(hits_data(:,1),1,size(hits_data,2));     % Subtract the 1st col (0 mM);
hits_data(1,:) = [];
hits_data(:,1) = [];
hits_orfs(1) = [];

hits_orfs = strtrim(upper(hits_orfs));
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];



[missing, ix] = setdiff(hits_orfs, tested_orfs);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
[~,ind3,ind4] = intersect(hits_treatments, treatments);
uluisik_koc_2011.data(ind2,ind4) = hits_data(ind1,ind3);

% Load data 2
[data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Uluisik~Koc/mmc3.xlsx');
hits_orfs = data.raw(4:end,1);
hits_data = data.raw(4:end,2:8);
hits_treatments = {'BA, 20 mM','BA, 30 mM','BA, 40 mM','BA, 50 mM','BA, 60 mM','BA, 70 mM'};

hits_data(strcmp('+',hits_data)) = {1};
hits_data(strcmp('-',hits_data)) = {0};

hits_data = cell2mat(hits_data);

inds = find(cellfun(@isempty, hits_orfs) | cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_data = hits_data - repmat(hits_data(1,:),size(hits_data,1),1);    % Subtract the 1st row (WT)
hits_data = hits_data - repmat(hits_data(:,1),1,size(hits_data,2));     % Subtract the 1st col (0 mM);
hits_data(1,:) = [];
hits_data(:,1) = [];
hits_orfs(1) = [];

hits_orfs = strtrim(upper(hits_orfs));
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);    % 3 ORFs missing (eliminated)
hits_orfs(ix) = [];     
hits_data(ix,:) = [];

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
[~,ind3,ind4] = intersect(hits_treatments, treatments);
uluisik_koc_2011.data(ind2,ind4) = hits_data(ind1,ind3);


save('~/Laboratory/Datasets/Phenotypes/2011_Uluisik~Koc/uluisik_koc_2011.mat','uluisik_koc_2011'); 

% Save data into database
dt = uluisik_koc_2011;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Baryshnikova~Myers, 2010
% DATA = baryshnikova_myers_2010
baryshnikova_myers_2010.pmid = 21076421;

phenotypes = {'growth'};
treatments = {'standard'};

% Load data
[data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2010_Baryshnikova~Myers/Supplementary_data_1_SMF_standard_100209.xlsx');

hits_orfs = data.raw(:,1);
hits_data = data.raw(:,2);

inds = find(cellfun(@isempty, hits_orfs) | cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

% Eliminate the data for non-deletion strains
t = regexp(hits_orfs,'_','split');
inds = find(cellfun(@length, t) > 1);
hits_orfs(inds) = [];
hits_data(inds,:) = [];

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

inds = find(~cellfun(@isnumeric, hits_data));
hits_data = cell2mat(hits_data);

hits_orfs = upper(strtrim(tested_orfs));
[t,t2] = grpstats(hits_data, hits_orfs,{'mean','gname'});

baryshnikova_myers_2010.orfs = t2;
baryshnikova_myers_2010.data = t;
baryshnikova_myers_2010.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2010_Baryshnikova~Myers/baryshnikova_myers_2010.mat','baryshnikova_myers_2010'); 

% Save data into database
dt = baryshnikova_myers_2010;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Gardocki~Lopes, 2005
% DATA = gardocki_lopes_2005
gardocki_lopes_2005.pmid = 15755922;

phenotypes = {'expression of PIS1'};
treatments = {'glucose 2%','glycerol 3%'};

% Load data
[data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Gardocki~Lopes/120_GENES_AFFECTING_PIS1.xlsx');

orfs = data.raw(8:end,2);
raw_data = data.raw(8:end,[11 14]);

inds = find(cellfun(@isempty, orfs) | cellfun(@isnumeric, orfs));
orfs(inds) = [];
raw_data(inds,:) = [];

inds = find(~strncmp('Y', orfs,1));
orfs(inds) = [];
raw_data(inds,:) = [];

% Transform symbols to numbers
% Glucose: + -> -1, - -> 0
% Glycerol: - -> 0, + -> +1

inds = find(strcmp('+', raw_data(:,1)));
raw_data(inds,1) = {-1};
inds = find(strcmp('-', raw_data(:,1)));
raw_data(inds,1) = {0};
inds = find(strcmp('+', raw_data(:,2)));
raw_data(inds,2) = {1};
inds = find(strcmp('-', raw_data(:,2)));
raw_data(inds,2) = {0};
inds = find(strcmp('NT', raw_data));
raw_data(inds) = {NaN};

raw_data = cell2mat(raw_data);

orfs = upper(strtrim(orfs));
[t,t2] = grpstats(raw_data, orfs,{'mean','gname'});

gardocki_lopes_2005.orfs = t2;
gardocki_lopes_2005.data = t;
gardocki_lopes_2005.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2005_Gardocki~Lopes/gardocki_lopes_2005.mat','gardocki_lopes_2005'); 

% Save data into database
dt = gardocki_lopes_2005;

% datasets = get_datasets_for_paper(dt);
% 
% [~,database_ix] = sortrows(datasets.names);
% [~,ph_ix] = sort(dt.ph);
ph_ix = 1:length(dt.ph);
% 
% % Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
% datasets.names(database_ix,:)
% dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, [575 576]);

%% Narayanaswamy~Marcotte, 2006
% DATA = narayanaswamy_marcotte_2006
narayanaswamy_marcotte_2006.pmid = 16507139;

phenotypes = {'expression of PIS1'};
treatments = {'standard'};

% Load data
[data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Narayanaswamy~Marcotte/gb-2006-7-1-r6-s2.xlsx');

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

save('Datasets/Phenotypes/2006_Narayanaswamy~Marcotte/narayanaswamy_marcotte_2006.mat','narayanaswamy_marcotte_2006'); 

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


%% Chamilos~Kontoyiannis, 2008
% DATA = chamilos_kontoyiannis_2008
chamilos_kontoyiannis_2008.pmid = 18212113;

phenotypes = {'growth'};
treatments = {'gliotoxin'};

% Load tested
tested_orfs = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Chamilos~Kontoyiannis/tested_orfs.txt','%s');

% Validate tested
expr = 'Y[A-P][RL][0-9]{3}[CW](-[ABC])*';
inds = find(cellfun(@isempty, regexpi(tested_orfs, expr)));

tested_orfs = unique(upper(strtrim(tested_orfs)));

% Load data
fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Chamilos~Kontoyiannis/data_genenames.txt','r');
C = textscan(fid, '%s %f','delimiter','\t');
fclose(fid);

genenames = C{1};
raw_data = C{2};

% Normalize to WT
inds = find(strcmp('WT', genenames));
raw_data = raw_data./raw_data(inds)-1;

genenames(inds) = [];
raw_data(inds) = [];

orfs = genename2orf_sgd(genenames);
orfs = upper(strtrim(orfs));

inds = find(cellfun(@isempty, orfs) | cellfun(@isnumeric, orfs));
orfs(inds) = [];
raw_data(inds,:) = [];

inds = find(~strncmp('Y', orfs,1));
orfs(inds) = [];
raw_data(inds,:) = [];

orfs = upper(strtrim(orfs));

missing_orfs = setdiff(orfs, tested_orfs);
tested_orfs = [tested_orfs; missing_orfs];  % Adding 2 missing strains to the list of tested strains.

chamilos_kostoyiannis_2008.orfs = tested_orfs;
chamilos_kostoyiannis_2008.data = zeros(length(tested_orfs),1);
[~,ind1,ind2] = intersect(tested_orfs, orfs);
chamilos_kostoyiannis_2008.data(ind1) = raw_data(ind2);
chamilos_kostoyiannis_2008.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2008_Chamilos~Kontoyiannis/chamilos_kontoyiannis_2008.mat','chamilos_kontoyiannis_2008'); 

% Save data into database
dt = chamilos_kontoyiannis_2008;

datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% % Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Ju~Xie, 2008
% DATA = ju_xie_2008
ju_xie_2008.pmid = 18070918;

phenotypes = {'abundance of Rpn4'};
treatments = {'standard'};

% Load tested
[data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Ju~Xie/mat_a_041902.xlsx','mat_a_041902');
tested_orfs = data.raw(3:end,2);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

% Validate tested
expr = 'Y[A-P][RL][0-9]{3}[CW](-[ABC])*';
inds = find(cellfun(@isempty, regexpi(tested_orfs, expr)));

tested_orfs(strcmp('YLR287-A', tested_orfs)) = {'YLR287C-A'};

tested_orfs = unique(upper(strtrim(tested_orfs)));

% Load data
fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Ju~Xie/hits.txt','r');
C = textscan(fid, '%s');
fclose(fid);

genenames = C{1};

orfs = genename2orf_sgd(genenames);

orfs(strcmp('cvt9',orfs)) = {'YPR049C'};

inds = find(cellfun(@isempty, orfs) | cellfun(@isnumeric, orfs));
orfs(inds) = [];

inds = find(~strncmp('Y', orfs,1));
orfs(inds) = [];

orfs = upper(strtrim(orfs));
data = ones(size(orfs));

missing_orfs = setdiff(orfs, tested_orfs);

ju_xie_2008.orfs = tested_orfs;
ju_xie_2008.data = zeros(length(tested_orfs),1);
[~,ind1,ind2] = intersect(tested_orfs, orfs);
ju_xie_2008.data(ind1) = data(ind2);
ju_xie_2008.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2008_Ju~Xie/ju_xie_2008.mat','ju_xie_2008'); 

% Save data into database
dt = ju_xie_2008;

datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% % Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

