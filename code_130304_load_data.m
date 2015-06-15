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


%% Cooper~Fields, 2010
% DATA = cooper_fields_2010

[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2010_Cooper~Fields/SupplementalTable4.xlsx');

cooper_fields_2010.orfs = upper(data.raw(2:end,1));
cooper_fields_2010.ph = data.raw(1,3:end)';

inds = find(strcmp('-', data.raw));
data.raw(inds) = {NaN};

cooper_fields_2010.data = cell2mat(data.raw(2:end,3:end));
cooper_fields_2010.desc = {'Supplemental Table 4. Log2 transformed ratios representing fold change compared to average for each identified amino acid in each of nearly 4500 samples.'};
cooper_fields_2010.pmid = 20610602;


% Identical ORFs with multiple values -> average
[t,t2] = grpstats(cooper_fields_2010.data, cooper_fields_2010.orfs, {'gname','mean'});
cooper_fields_2010.orfs = t;
cooper_fields_2010.data = t2;

save('Datasets/Phenotypes/2010_Cooper~Fields/cooper_fields_2010.mat','cooper_fields_2010');

% Save data into database
insert_data_into_db(cooper_fields_2010);

%% Begley~Samson, 2002
% DATA = begley_samson_2002
% NOTE = data sorting fixed (lower=slower)

begley_samson_2002.pmid = 12496357;

phenotype = {'Growth, dilution spot intensity'};

treatments = {'4NQO','MMS','t-BuOOH','UV'};
doses(1,:) = {'0.1 ug/ml','0.2 ug/ml','0.3 ug/ml','0.4 ug/ml'};
doses(2,:) = {'0.01%','0.02%','0.025%','0.03%'};
doses(3,:) = {'0.25 mM','0.5 mM','0.75 mM','1 mM'};
doses(4,:) = {'25 J/m2','50 J/m2','75 J/m2','100 J/m2'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2002_Begley~Samson/ORIG130404_Begley2001raw.xlsx','Sheet1');

% Find the column with the systematic ORF name
ind_orf = strmatch('ORF', data.txt(1,:));

% Find the column with the name of the condition
ind_condition = strmatch('Treatment', data.txt(1,:));

% Find the columns with the relevant data
% This is the ratio of the growth of this spot on the untreated plate to 
% the growth on a treated plate divided by the same ratio for the highest 
% growing WT colony. In the paper, a value of greater than 1.5 was
% considered sensitive; a value of less than .66 was considered resistant.

ind_data1 = strmatch('exp 1 highest control/dose 1', data.txt(1,:)); ind_data1 = [ind_data1:ind_data1+3];
ind_data2 = strmatch('exp 1 lowest control/dose 1', data.txt(1,:)); ind_data2 = [ind_data2:ind_data2+3];

ind_data3 = strmatch('exp 2 highest control/dose 1', data.txt(1,:)); ind_data3 = [ind_data3:ind_data3+3];
ind_data4 = strmatch('exp 2 lowest control/dose 1', data.txt(1,:)); ind_data4 = [ind_data4:ind_data4+3];

ind_data5 = strmatch('exp 3 highest control/dose 1', data.txt(1,:)); ind_data5 = [ind_data5:ind_data5+3];
ind_data6 = strmatch('exp 3 lowest control/dose 1', data.txt(1,:)); ind_data6 = [ind_data6:ind_data6+3];

% Extract the most relevant information, ignore the rest.
all_strings = cellfun(@num2str, data.raw(:,ind_orf),'UniformOutput',0);
ind_not_empty = setdiff(1:size(data.raw,1), strmatch('NaN', all_strings,'exact'))';

data2.orfs = data.raw(ind_not_empty,ind_orf);
data2.conditions = data.raw(ind_not_empty, ind_condition);
data2.data = data.raw(ind_not_empty, [ind_data1 ind_data2 ind_data3 ind_data4 ind_data5 ind_data6]);
data2.experiments = data.raw(1,[ind_data1 ind_data2 ind_data3 ind_data4 ind_data5 ind_data6]);

% Eliminate the first row (table headers)
data2.orfs(1) = [];
data2.conditions(1) = [];
data2.data(1,:) = [];

% Find weird values, before converting the data cell array into a matrix
t = find(~cellfun(@isnumeric, data2.data));
data2.data(t) = {NaN};

% Convert the data cell array into a matrix
data2.data = cell2mat(data2.data);

% Take the reciprocal of the values, such that higher numbers indicate
% higher growth, and viceversa
data2.data = 1./data2.data;


% Average the 3 replicates (exp 1-3), as well as the data calculated with
% respect to the highest and lowest controls
data2.data_avg(:,1) = nanmean(data2.data(:,1:4:end),2);
data2.data_avg(:,2) = nanmean(data2.data(:,2:4:end),2);
data2.data_avg(:,3) = nanmean(data2.data(:,3:4:end),2);
data2.data_avg(:,4) = nanmean(data2.data(:,4:4:end),2);


% Divide by condition

for ic = 1 : length(treatments)
    
    inds_cond = strmatch(treatments{ic}, data2.conditions);
    
    data3(ic).orfs = upper(data2.orfs(inds_cond));
    data3(ic).data = data2.data_avg(inds_cond,:);
    
    % Average data for identical ORFs that appear multiple times in each
    % dataset
    [t,t2] = grpstats(data3(ic).data, data3(ic).orfs, {'gname','mean'});
    data4(ic).orfs = t;
    data4(ic).data = t2;
    
    % Make sure the final list of all ORFs is comprehensive
    if ic == 1
        begley_samson_2002.orfs = data4(ic).orfs;
        begley_samson_2002.data = zeros(length(data4(ic).orfs),length(doses(:)))+NaN;
        begley_samson_2002.ph = cell(length(doses(:)),1);
    end
    inds = find(~ismember(data4(ic).orfs, begley_samson_2002.orfs));
    begley_samson_2002.orfs = [begley_samson_2002.orfs; data4(ic).orfs(inds)];
    
    % Now it is possible to intersect this dataset with the comprehensive
    % list, align and append all the values
    [C, ia,ib] = intersect(begley_samson_2002.orfs, data4(ic).orfs);
    begley_samson_2002.data(ia,4*(ic-1)+1:4*(ic-1)+4) = data4(ic).data(ib,:);
    begley_samson_2002.ph(4*(ic-1)+1:4*(ic-1)+4) = strcat(phenotype, {'; '}, treatments(ic),{', '}, doses(ic,:))';
    
    
end

% A few genenames to rename
genenames = {'RAD14','REV1','MAG1'};
orfs = genename2orf(genenames,'noannot');
[t,ind1,ind2] = intersect(begley_samson_2002.orfs, genenames);
begley_samson_2002.orfs(ind1) = orfs(ind2);

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', begley_samson_2002.orfs,1));
begley_samson_2002.orfs(inds) = [];
begley_samson_2002.data(inds,:) = [];


save('Datasets/Phenotypes/2002_Begley~Samson/begley_samson_2002.mat','begley_samson_2002');


% Save data into database
dt = begley_samson_2002;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Warringer~Blomberg, 2003
% DATA = warringer_blomberg_2003

warringer_blomberg_2003.pmid = 14676322;

phenotypes = {'Growth, lag phase';'Growth, exponential growth rate';'Growth, saturation level'};
treatments = {'NaCl, 0.85 M'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2003_Warringer~Blomberg/ORIG130305_LPI NaCl.xlsx','LPI');

data2.orfs = upper(data.raw(5:end,1));
data2.data = cell2mat(data.raw(5:end,2:4));

% Eliminate anything that doesn't look like an ORF
inds = setdiff(1:length(data2.orfs),strmatch('Y', data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
warringer_blomberg_2003.orfs = t;
warringer_blomberg_2003.data = t2;
warringer_blomberg_2003.ph = strcat(phenotypes, {'; '}, treatments);

save('Datasets/Phenotypes/2003_Warringer~Blomberg/warringer_blomberg_2003.mat','warringer_blomberg_2003');

% Save data into database
datasets = get_datasets_for_paper(warringer_blomberg_2003);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).shortname;
    datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(warringer_blomberg_2003.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.

insert_data_into_db(warringer_blomberg_2003, ph_ix, datasets_ids(database_ix));


%% Hartman~Tippery, 2004
% DATA = hartman_tippery_2004

hartman_tippery_2004.pmid = 15239834;
hartman_tippery_2004.desc = {'The loaded values are Z-scores with respect to the wild-type in a given condition (UNT or HU). The final values (normalized to UNT) are the difference between treated and untreated z-scores.'};

phenotypes = {'Growth, AUGC'};
treatments = {'UNT';'HU, 50 mM';'HU, 150 mM'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2004_Hartman~Tippery/gb-2004-5-7-r49-s7.xlsx','data');

ind_orf = strmatch('ORF', data.raw(1,:));
data2.orfs = data.raw(2:end, ind_orf);

ind_data1 = strmatch('No HU- Growth Index', data.raw(1,:));
ind_data2 = strmatch('50mM HU Growth Index', data.raw(1,:));
ind_data3 = strmatch('150 mM HU Growth Index', data.raw(1,:));

data2.data = data.raw(2:end, [ind_data1 ind_data2 ind_data3]);

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

inds = setdiff(1:length(data2.orfs),strmatch('Y', data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

data2.orfs = upper(data2.orfs);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Normalize by UNTREATED sample
data2.data(:,2) = data2.data(:,2) - data2.data(:,1);
data2.data(:,3) = data2.data(:,3) - data2.data(:,1);
data2.data(:,1) = [];
treatments(1) = [];


% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
hartman_tippery_2004.orfs = t;
hartman_tippery_2004.data = t2;
hartman_tippery_2004.ph = strcat(phenotypes, {'; '}, treatments);

save('Datasets/Phenotypes/2004_Hartman~Tippery/hartman_tippery_2004.mat','hartman_tippery_2004');

% Save data into database
dt = hartman_tippery_2004;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).shortname;
    datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Dunn~Jensen, 2006
% DATA = dunn_jensen_2006


dunn_jensen_2006.source = {'http://www.molbiolcell.org/content/suppl/2005/11/02/E05-06-0585.DC1/Supp_Table3.xls'};
dunn_jensen_2006.downloaddate = {'2013-03-06'};
dunn_jensen_2006.pmid = 16267274;

phenotypes = {'Growth, log2 hybridization ratio'};
treatments = {'EtBr, 25 ug/ml'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2006_Dunn~Jensen/Supp_Table3.xlsx','Sheet1');

ind_orf = strmatch('Systemic name', data.raw(12,:));
data2.orfs = upper(data.raw(14:end, ind_orf));

ind_data1 = strmatch('(-EtBr/+EtBr)', data.raw(12,:));

data2.data = data.raw(14:end, ind_data1);

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

inds = setdiff(1:length(data2.orfs),strmatch('Y', data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);


% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
dunn_jensen_2006.orfs = t;
dunn_jensen_2006.data = 1./t2;  % reverse the ratio so that the lower the value, the sicker the mutant.
dunn_jensen_2006.ph = strcat(phenotypes, {'; '}, treatments);

save('Datasets/Phenotypes/2006_Dunn~Jensen/dunn_jensen_2006.mat','dunn_jensen_2006');

% Save data into database
dt = dunn_jensen_2006;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).shortname;
    datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Yoshikawa~Shimizu, 2009
% DATA = yoshikawa_shimizu_2009

yoshikawa_shimizu_2009.pmid = 19054128;
yoshikawa_shimizu_2009.desc = {'The final values are treated divided by untreated.'};

phenotypes = {'Growth, exponential growth rate'};
treatments = {'UNT';'EtOH, 5%';'EtOH, 8%';'NaCl, 1 M'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2009_Yoshikawa~Shimizu/FYR_456_sm_tableS1.xlsx','data');

data2.orfs = upper(data.raw(3:end, 1));

ind_data = [3 4 6 7 9 10 12 13];
data2.data = data.raw(3:end, ind_data);

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

inds = setdiff(1:length(data2.orfs),strmatch('Y', data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average the 2 replicates
data3.orfs = data2.orfs;
data3.data(:,1) = nanmean(data2.data(:,1:2),2);
data3.data(:,2) = nanmean(data2.data(:,3:4),2);
data3.data(:,3) = nanmean(data2.data(:,5:6),2);
data3.data(:,4) = nanmean(data2.data(:,7:8),2);

% Normalize by UNTREATED sample
data3.data(:,2) = data3.data(:,2)./data3.data(:,1);
data3.data(:,3) = data3.data(:,3)./data3.data(:,1);
data3.data(:,4) = data3.data(:,4)./data3.data(:,1);
data3.data(:,1) = [];
treatments(1) = [];

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data3.data, data3.orfs, {'gname','mean'});
yoshikawa_shimizu_2009.orfs = t;
yoshikawa_shimizu_2009.data = t2;
yoshikawa_shimizu_2009.ph = strcat(phenotypes, {'; '}, treatments);

save('Datasets/Phenotypes/2009_Yoshikawa~Shimizu/yoshikawa_shimizu_2009.mat','yoshikawa_shimizu_2009');

% Save data into database
dt = yoshikawa_shimizu_2009;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).shortname;
    datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix([2 3 1]),:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix([2 3 1])));

%% Alamgir~Golshani, 2010
% DATA = alamgir_golshani_2010

alamgir_golshani_2010.source = {'http://www.biomedcentral.com/content/supplementary/1472-6769-10-6-s1.xls'};
alamgir_golshani_2010.downloaddate = {'2013-03-06'};
alamgir_golshani_2010.pmid = 20691087;

phenotypes = {'Growth, colony size'};
treatments = {'3-AT, 22 mg/ml';'cycloheximide, 45 ng/ml';'streptomycin, 40 mg/ml';'neomycin, 5.5 mg/ml'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2010_Alamgir~Golshani/1472-6769-10-6-s1.xlsx','Raw genome-wide data');

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,1)),strmatch('Y', data.raw(:,1)));
data.raw(inds,:) = [];

data2.orfs = upper(data.raw(:, 1));
data2.data = data.raw(:, 3:14);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Eliminate zeros
data2.data(data2.data == 0) = NaN;

% Average the 3 replicates
data3.orfs = data2.orfs;
data3.data(:,1) = nanmean(data2.data(:,1:3),2);
data3.data(:,2) = nanmean(data2.data(:,4:6),2);
data3.data(:,3) = nanmean(data2.data(:,7:9),2);
data3.data(:,4) = nanmean(data2.data(:,10:12),2);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data3.data, data3.orfs, {'gname','mean'});
alamgir_golshani_2010.orfs = t;
alamgir_golshani_2010.data = t2;
alamgir_golshani_2010.ph = strcat(phenotypes, {'; '}, treatments);

save('Datasets/Phenotypes/2010_Alamgir~Golshani/alamgir_golshani_2010.mat','alamgir_golshani_2010');

% Save data into database
dt = alamgir_golshani_2010;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).shortname;
    datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Chavel~Cullen, 2010
% DATA = chavel_cullen_2010

chavel_cullen_2010.source = {'http://www.plosgenetics.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pgen.1000883.s011'};
chavel_cullen_2010.downloaddate = {'2013-03-06'};
chavel_cullen_2010.pmid = 20333241;

phenotypes = {'Msb2p secretion'};
treatments = {''};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2010_Chavel~Cullen/journal.pgen.1000883.s011.xlsx','Complete Screen');

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,2)),strmatch('Y', data.raw(:,2)));
data.raw(inds,:) = [];

% Eliminate false positives
inds = find(cellfun(@isempty, data.raw(:,7))==0);
data.raw(inds,6) = NaN;


data2.orfs = upper(data.raw(:,2));
data2.data = data.raw(:,6);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Eliminate zeros
data2.data(data2.data == 0) = NaN;

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
chavel_cullen_2010.orfs = t;
chavel_cullen_2010.data = t2;
chavel_cullen_2010.ph = phenotypes;

save('Datasets/Phenotypes/2010_Chavel~Cullen/chavel_cullen_2010.mat','chavel_cullen_2010');

% Save data into database
dt = chavel_cullen_2010;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).shortname;
    datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Kaeberlein~Kennedy, 2005
% DATA = kaeberlein_kennedy_2005


kaeberlein_kennedy_2005.source = {'http://www.sciencemag.org/content/suppl/2005/11/15/310.5751.1193.DC1/Kaeberlein.SOM.pdf'};
kaeberlein_kennedy_2005.notes = {'SOM downloaded from Science website, then converted to Word and table S1 saved to Excel.'};
kaeberlein_kennedy_2005.downloaddate = {'2013-02-11'};
kaeberlein_kennedy_2005.pmid = 16293764;

phenotypes = {'Replicative life span'};
treatments = {''};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2005_Kaeberlein~Kennedy/Kaeberlein.SOM_tableS1.xlsx','Sheet1');

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

inds = strmatch('(YHR034W)', data.raw(:,1));
data.raw(inds,1) = {'YHR034W'};

inds = setdiff(1:length(data.raw(:,1)),strmatch('Y', data.raw(:,1)));
data.raw(inds,:) = [];


data2.orfs = upper(data.raw(:,1));
data2.data = data.raw(:,4:7);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Eliminate data with no info about replicates (N)
inds = find(isnan(data2.data(:,2)));
data2.data(inds,:) = NaN;


% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
kaeberlein_kennedy_2005.orfs = t;
kaeberlein_kennedy_2005.data = t2(:,1)./t2(:,3);
kaeberlein_kennedy_2005.ph = phenotypes;

save('Datasets/Phenotypes/2005_Kaeberlein~Kennedy/kaeberlein_kennedy_2005.mat','kaeberlein_kennedy_2005');

% Save data into database
dt = kaeberlein_kennedy_2005;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).shortname;
    datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Mozdy~Cech, 2008
% DATA = mozdy_cech_2008

mozdy_cech_2008.source = {'http://mcb.asm.org/content/suppl/2008/05/16/28.12.4152.DC1/MozdySupplementalTable1.pdf'};
mozdy_cech_2008.notes = {'SOM downloaded from MCB website, then converted to Word and saved to Excel.'};
mozdy_cech_2008.downloaddate = {'2013-03-07'};
mozdy_cech_2008.pmid = 18411302;

phenotypes = {'TLC1 abundance (relative to U2)';'TLC1 abundance (relative to ACT1)'};
treatments = {''};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2008_Mozdy~Cech/MozdySupplementalTable1.xlsx','Sheet1');

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,2)),strmatch('Y', data.raw(:,2)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,2) = cellfun(@strtrim, data.raw(:,2),'UniformOutput',0);

data2.orfs = upper(data.raw(:,2));
data2.data = data.raw(:,7:8);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
mozdy_cech_2008.orfs = t;
mozdy_cech_2008.data = t2;
mozdy_cech_2008.ph = phenotypes;

save('Datasets/Phenotypes/2008_Mozdy~Cech/mozdy_cech_2008.mat','mozdy_cech_2008');

% Save data into database
dt = mozdy_cech_2008;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).shortname;
    datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Fillingham~Andrews, 2009
% DATA = fillingham_andrews_2009


fillingham_andrews_2009.source = {'http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276509004614/1-s2.0-S1097276509004614-mmc3.xls/272198/FULL/S1097276509004614/c47432a0acc5d8c66fa4ae36e559c4f9/mmc3.xls'};
fillingham_andrews_2009.downloaddate = {'2013-03-07'};
fillingham_andrews_2009.pmid = 19683497;

phenotypes = {'HTA1 expression, z-score of log2 GFP:RFP'};
treatments = {''};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2009_Fillingham~Andrews/mmc3.xlsx','Sheet1');

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,1)),strmatch('Y', data.raw(:,1)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,1) = cellfun(@strtrim, data.raw(:,1),'UniformOutput',0);

data2.orfs = upper(data.raw(:,1));
data2.data = data.raw(:,2);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
fillingham_andrews_2009.orfs = t;
fillingham_andrews_2009.data = t2;
fillingham_andrews_2009.ph = phenotypes;

save('Datasets/Phenotypes/2009_Fillingham~Andrews/fillingham_andrews_2009.mat','fillingham_andrews_2009');

% Save data into database
dt = fillingham_andrews_2009;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).shortname;
    datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Hess~Caudy, 2009
% DATA = hess_caudy_2009

hess_caudy_2009.source = {'http://www.plosgenetics.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pgen.1000407.s005'};
hess_caudy_2009.downloaddate = {'2013-03-07'};
hess_caudy_2009.pmid = 19300474;

phenotypes = {'Petite frequency (fraction)';'Growth, doubling time (hours); glycerol';'Growth, doubling time (hours); glucose'; ...
    'Growth, saturation level (OD); glycerol';'Growth, saturation level (OD); glucose'};
treatments = {''};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2009_Hess~Caudy/journal.pgen.1000407.s005.xlsx','Sheet1');

% Get indices of the data columns
inds = find(cellfun(@isnumeric, data.raw(28,:)));
data.raw(28,inds) = {'NaN'};
ind_data1 = strmatch('Mean Petite Frequency with Respect to Wild-type', data.raw(28,:));
ind_data2 = strmatch('Doubling Time Mean (Glycerol) [hours]', data.raw(28,:));
ind_data3 = strmatch('Doubling Time Mean (Glucose) [hours]', data.raw(28,:));
ind_data4 = strmatch('Staturation Mean (Glycerol) [OD]', data.raw(28,:));
ind_data5 = strmatch('Staturation Mean (Glucose) [OD]', data.raw(28,:));
ind_data = [ind_data1; ind_data2; ind_data3; ind_data4; ind_data5];

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,2)),strmatch('Y', data.raw(:,2)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,2) = cellfun(@strtrim, data.raw(:,2),'UniformOutput',0);

data2.orfs = upper(data.raw(:,2));
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

data2.data(:,1) = data2.data(:,1)./100; % Convert percentages in fractions


% Normalize by GLUCOSE
data2.data(:,2) = data2.data(:,3)./data2.data(:,2);     % Taking the reciprocal of glycerol/glucose such that, at the end, lower numbers correspond to slower growth.
data2.data(:,4) = data2.data(:,4)./data2.data(:,5);     % For saturation OD, this is not necessary. Lower OD means slower/sicker mutant.
data2.data(:,[3,5]) = [];
phenotypes([3,5]) = [];

phenotypes(2) = {'Growth, doubling time (h); glucose/glycerol'};
phenotypes(3) = {'Growth, saturation level (OD); glycerol/glucose'};

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
hess_caudy_2009.orfs = t;
hess_caudy_2009.data = t2;
hess_caudy_2009.ph = phenotypes;

save('Datasets/Phenotypes/2009_Hess~Caudy/hess_caudy_2009.mat','hess_caudy_2009');

% Save data into database
dt = hess_caudy_2009;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).shortname;
    datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix([1 3 2]),:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix([1 3 2])));


%% Bleackley~MacGillivray, 2011
% DATA = bleackley_macgillivray_2011


bleackley_macgillivray_2011.source = {'Mark Bleackley, Ross MacGillivray'};
bleackley_macgillivray_2011.notes = {'File received from the authors because online SOM didn''t have ORFs, only genenames.'};
bleackley_macgillivray_2011.downloaddate = {'2013-03-11'};
bleackley_macgillivray_2011.pmid = 21212869;

phenotypes = {'Growth, colony size'};
treatments = {'Iron, Fe(NH4)2(SO4)2 (10 mM)';'Copper, CuCl2 (7 mM)';...
    'Manganese, MnCl2 (4 mM)';'Nickel, NiCl2 (3 mM)';'Zinc, ZnCl2 (7 mM)';...
    'Cobalt, CoCl2 (2.5 mM)'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2011_Bleackley~MacGillivray/metallomicsbleackley raw data.xls','rawdata');

% Get indices of the data columns
ind_data = 3:2:13;


% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,1)),strmatch('Y', data.raw(:,1)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,1) = cellfun(@strtrim, data.raw(:,1),'UniformOutput',0);

data2.orfs = upper(data.raw(:,1));
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
bleackley_macgillivray_2011.orfs = t;
bleackley_macgillivray_2011.data = t2;
bleackley_macgillivray_2011.ph = strcat(phenotypes, {'; '}, treatments);

save('Datasets/Phenotypes/2011_Bleackley~MacGillivray/bleackley_macgillivray_2011.mat','bleackley_macgillivray_2011');

% Save data into database
dt = bleackley_macgillivray_2011;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).short_name;
    datasets_names{i,3} = datasets(i).dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Yoshikawa~Shimizu, 2011
% DATA = yoshikawa_shimizu_2011


yoshikawa_shimizu_2011.source = {'http://onlinelibrary.wiley.com/store/10.1002/yea.1843/asset/supinfo/yea_1843_supportinginforTS1.xls?v=1&s=ebd728c627c24d863139e83ed47ca3e94b22c39e'};
yoshikawa_shimizu_2011.downloaddate = {'2013-03-08'};
yoshikawa_shimizu_2011.pmid = 21341307;

phenotypes = {'Growth, exponential growth rate (h^-1)'};
treatments = {''};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2011_Yoshikawa~Shimizu/yea_1843_supportinginforTS1.xlsx','Deletion');

% Get indices of the data columns
ind_data = 5; % "Average"

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,2)),strmatch('Y', data.raw(:,2)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,2) = cellfun(@strtrim, data.raw(:,2),'UniformOutput',0);

data2.orfs = upper(data.raw(:,2));
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
yoshikawa_shimizu_2011.orfs = t;
yoshikawa_shimizu_2011.data = t2;
yoshikawa_shimizu_2011.ph = phenotypes;

save('Datasets/Phenotypes/2011_Yoshikawa~Shimizu/yoshikawa_shimizu_2011.mat','yoshikawa_shimizu_2011');

% Save data into database
dt = yoshikawa_shimizu_2011;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).reporter;
    datasets_names{i,3} = datasets(i).short_name;
    datasets_names{i,4} = datasets(i).dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3 4]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Svensson~Samson, 2011
% DATA = svensson_samson_2011


svensson_samson_2011.source = {'http://www.biomedcentral.com/content/supplementary/1752-0509-5-157-s1.xls'};
svensson_samson_2011.downloaddate = {'2013-03-08'};
svensson_samson_2011.pmid = 21978764;

phenotypes = {'Growth, GI50'};
treatments = {'MMS'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2011_Svensson~Samson/1752-0509-5-157-s1.xlsx','2. Gi50 and R2 all strains');

% Get indices of the data columns
ind_data = 4;

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,2)),strmatch('Y', data.raw(:,2)));
data.raw(inds,:) = [];

% Separate deletions from DAMP strains. Personal communication from
% Peter Svensson: deletions are on plates 1-57, DAMPs are on plates
% 301-311

dels = 'ABCDEFGH';
plate = data.raw(:,1);
for i = 1 : length(dels)
    plate = strtok(plate, dels(i));
end
plate = cellfun(@str2num, plate);

inds = find(plate > 57);
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,2) = cellfun(@strtrim, data.raw(:,2),'UniformOutput',0);

data2.orfs = upper(data.raw(:,2));
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
svensson_samson_2011.orfs = t;
svensson_samson_2011.data = t2;
svensson_samson_2011.ph = strcat(phenotypes, {'; '} , treatments);

save('Datasets/Phenotypes/2011_Svensson~Samson/svensson_samson_2011.mat','svensson_samson_2011');

% Save data into database
dt = svensson_samson_2011;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).reporter;
    datasets_names{i,3} = datasets(i).short_name;
    datasets_names{i,4} = datasets(i).dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3 4]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));


%% Peyroche~Plateau, 2012
% DATA = peyroche_plateau_2012


peyroche_plateau_2012.source = {'http://www.plosone.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pone.0036343.s004'};
peyroche_plateau_2012.downloaddate = {'2013-03-08'};
peyroche_plateau_2012.pmid = 22586468;
peyroche_plateau_2012.desc = {'Relative fitness defect: rf = log2(wt(Se)/mut(Se) - wt/mut + 1), where wt and mut are the generation times of the WT and mutant strains with and without Na2Se.'};

phenotypes = {'Growth, log2 ratio'};
treatments = {'Sodium selenide, 1 uM 16 h'; 'Sodium selenide, 2 uM 16 h'; ...
    'Sodium selenide, 1 uM 27 h'; 'Sodium selenide, 2 uM 27 h'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2012_Peyroche~Plateau/journal.pone.0036343.s004.xlsx','data');

% Get indices of the data columns
ind_data = 5:8;

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,1)),strmatch('Y', data.raw(:,1)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,1) = cellfun(@strtrim, data.raw(:,1),'UniformOutput',0);

data2.orfs = upper(data.raw(:,1));
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
peyroche_plateau_2012.orfs = t;
peyroche_plateau_2012.data = t2;
peyroche_plateau_2012.ph = strcat(phenotypes, {'; '}, treatments);

save('Datasets/Phenotypes/2012_Peyroche~Plateau/peyroche_plateau_2012.mat','peyroche_plateau_2012');

% Save data into database
dt = peyroche_plateau_2012;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).reporter;
    datasets_names{i,3} = datasets(i).conditionset;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));


%% Kapitzky~Krogan, 2010

kapitzky_krogan_2010.source = {'http://www.nature.com/msb/journal/v6/n1/extref/msb2010107-s2.xls'};
kapitzky_krogan_2010.downloaddate = {'2013-03-09'};
kapitzky_krogan_2010.pmid = 21179023;

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2010_Kapitzky~Krogan/msb2010107-s2.xlsx','S.cerevisiae D-scores');

phenotypes = {'Growth, colony size'};
treatments = data.raw(2,2:end);
doses = strcat({' ('}, cellfun(@num2str, data.raw(3,2:end),'UniformOutput',0), {' uM)'});

% Get indices of the data columns
ind_data = 2:length(compounds)+1;

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

orfs = genename2orf(data.raw(:,1),'noannot');
orfs(strcmp('SOY1',orfs)) = {'YBR194W'};    % Manually adjust the ambiguous genes.

inds = find(~strncmp('Y', orfs,1));
orfs(inds) = [];
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
orfs = cellfun(@strtrim, orfs,'UniformOutput',0);

data2.orfs = upper(orfs);
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
kapitzky_krogan_2010.orfs = t;
kapitzky_krogan_2010.data = t2;
kapitzky_krogan_2010.ph = strcat(phenotypes, {'; '}, treatments, doses)';

save('Datasets/Phenotypes/2010_Kapitzky~Krogan/kapitzky_krogan_2010.mat','kapitzky_krogan_2010');

% Save data into database
dt = kapitzky_krogan_2010;

datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[1 2 3 4]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
[~,adj_ix] = sort([13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,2,3,1,8,4,5,6,7,9,10,11,12,33,34,35]);
datasets.names(database_ix(adj_ix),:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix(adj_ix)));

%% Kim~Cunningham, 2012
% DATA = kim_cunningham_2012

kim_cunningham_2012.pmid = 22511765;

phenotypes = {'Cell death, frequency of dead cells'};
treatments = {'tunicamycin, 2.5 ug/ml (Z-score ln)'; 'tunicamycin + FK506, 2.5 ug/ml + 1 ug/ml (Z-score ln)'; ...
    'dithiothreitol, 4 mM (Z-score ln)'; 'dithiothreitol, 4 mM, + FK506, 1 ug/ml (Z-score ln)'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2012_Kim~Cunningham/jbc.M112.363390-1.xlsx','sup table 1');

% Get indices of the data columns
ind_data = [14:15 17:18];

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,5)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,5)),strmatch('Y', data.raw(:,5)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,5) = cellfun(@strtrim, data.raw(:,5),'UniformOutput',0);

data2.orfs = upper(data.raw(:,5));
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
kim_cunningham_2012.orfs = t;
kim_cunningham_2012.data = t2;
kim_cunningham_2012.ph = strcat(phenotypes, {'; '}, treatments);

save('Datasets/Phenotypes/2012_Kim~Cunningham/kim_cunningham_2012.mat','kim_cunningham_2012');

% Save data into database
dt = kim_cunningham_2012;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    if isempty(datasets(i).reporter)
        datasets_names{i,2} = '';
    else
        datasets_names{i,2} = datasets(i).reporter;
    end
    datasets_names{i,3} = datasets(i).conditionset;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix([2 1 4 3]),:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix(([2 1 4 3]))));

%% Breslow~Weissman, 2008
% DATA = breslow_weissman_2008

breslow_weissman_2008.source = {'http://www.nature.com/nmeth/journal/v5/n8/extref/nmeth.1234-S4.xls'};
breslow_weissman_2008.downloaddate = {'2013-03-12'};
breslow_weissman_2008.pmid = 18622397;

phenotypes = {'Growth, flow cytometry'};
treatments = {''};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2008_Breslow~Weissman/nmeth.1234-S5.xlsx','Deletion growth rates');

% Get indices of the data columns
ind_data = find(strcmp('Median Growth Rate', data.raw(1,:)));

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,2)),strmatch('Y', data.raw(:,2)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,2) = cellfun(@strtrim, data.raw(:,2),'UniformOutput',0);

data2.orfs = upper(data.raw(:,2));
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
breslow_weissman_2008.orfs = t;
breslow_weissman_2008.data = t2;
breslow_weissman_2008.ph = phenotypes;

save('Datasets/Phenotypes/2008_Breslow~Weissman/breslow_weissman_2008.mat','breslow_weissman_2008');

% Save data into database
dt = breslow_weissman_2008;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    if isempty(datasets(i).reporter)
        datasets_names{i,2} = '';
    else
        datasets_names{i,2} = datasets(i).reporter;
    end
    datasets_names{i,3} = datasets(i).conditionset;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Ohya~Morishita, 2005
% DATA = ohya_morishita_2005


ohya_morishita_2005.source = {'http://scmd.gi.k.u-tokyo.ac.jp/datamine/pnas/mutant_analysis_2011_10_20.tab'};
ohya_morishita_2005.downloaddate = {'2013-03-12'};
ohya_morishita_2005.pmid = 16365294;

phenotypes = {'Growth, flow cytometry'};
treatments = {''};

[param.num,param.txt,param.raw] = xlsread('Datasets/Phenotypes/2005_Ohya~Morishita/parameter_names.xlsx','Sheet1');

data = importdata('Datasets/Phenotypes/2005_Ohya~Morishita/mutant_analysis_2011_10_20.tab');
data.orfs = data.textdata(2:end,1);
data.params = data.textdata(1,2:end)';

[t,ind1,ind2] = intersect(data.params, param.raw(:,1));
data.params_name(ind1,1) = param.raw(ind2,2);


% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.orfs));
data.orfs(inds) = [];
data.data(inds,:) = [];

inds = setdiff(1:length(data.orfs),strmatch('Y', data.orfs));
data.orfs(inds) = [];
data.data(inds,:) = [];

% Eliminate white spaces before/after ORF
data.orfs = cellfun(@strtrim, data.orfs,'UniformOutput',0);
data.orfs = upper(data.orfs);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data.data, data.orfs, {'gname','mean'});
ohya_morishita_2005.orfs = t;
ohya_morishita_2005.data = t2;
ohya_morishita_2005.ph = data.params_name;

save('Datasets/Phenotypes/2005_Ohya~Morishita/ohya_morishita_2005.mat','ohya_morishita_2005');


% NEED TO SEND EMAIL TO YOSHI OHYA WITH QUESTIONS
% % Save data into database
% dt = ohya_morishita_2005;
% 
% datasets = get_datasets_for_paper(dt);
% datasets_ids = zeros(length(datasets),1);
% datasets_names = cell(length(datasets),3);
% for i = 1 : length(datasets)
%     datasets_ids(i,1) = datasets(i).id;
%     datasets_names{i,1} = datasets(i).name;
%     if isempty(datasets(i).reporter)
%         datasets_names{i,2} = '';
%     else
%         datasets_names{i,2} = datasets(i).reporter;
%     end
%     datasets_names{i,3} = datasets(i).conditionset;
% end
% 
% [~,database_ix] = sortrows(datasets_names,[1 2 3]);
% [~,ph_ix] = sort(dt.ph);
% 
% % Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
% datasets_names(database_ix,:)
% dt.ph(ph_ix)
% 
% insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Xie~Huang, 2005
% DATA = xie_huang_2005


xie_huang_2005.source = {'Jing Huang'};
xie_huang_2005.downloaddate = {'2013-03-21'};
xie_huang_2005.pmid = 15883373;

phenotypes = {'Growth, spot intensity on cell microarray'};
treatments = {'Rapamycin, 10 nM'; 'Rapamycin, 30 nM'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2005_Xie~Huang/rapa_cellarray_all.xls','Rapa10nm');
[data2.num,data2.txt,data2.raw] = xlsread('Datasets/Phenotypes/2005_Xie~Huang/rapa_cellarray_all.xls','Rapa30nm');

% Get indices of the data columns
ind_data = find(strcmp('Rapa/DMSO', data.raw(1,:)));
ind_data2 = find(strcmp('Rapa/DMSO', data2.raw(1,:)));

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,3)));
data.raw(inds,:) = [];
inds = find(cellfun(@isnumeric, data2.raw(:,3)));
data2.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,3)),strmatch('Y', data.raw(:,3)));
data.raw(inds,:) = [];
inds = setdiff(1:length(data2.raw(:,3)), strmatch('Y', data2.raw(:,3)));
data2.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,3) = cellfun(@strtrim, data.raw(:,3),'UniformOutput',0);
data2.raw(:,3) = cellfun(@strtrim, data2.raw(:,3),'UniformOutput',0);

data3.orfs = upper(data.raw(:,3));
data3.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data3.data)==0);
data3.data(inds) = {NaN};
data3.data = cell2mat(data3.data);

data4.orfs = upper(data2.raw(:,3));
data4.data = data2.raw(:,ind_data2);
inds = find(cellfun(@isnumeric, data4.data)==0);
data4.data(inds) = {NaN};
data4.data = cell2mat(data4.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data3.data, data3.orfs, {'gname','mean'});
xie_huang_2005.orfs = t;
xie_huang_2005.data = t2;
xie_huang_2005.ph = strcat(phenotypes, '; ', treatments);

[t,t2] = grpstats(data4.data, data4.orfs, {'gname','mean'});

[t3,ind1,ind2] = intersect(xie_huang_2005.orfs, t);
xie_huang_2005.data(ind1,end+size(t2,2)) = t2(ind2,:);

% The list contains both non-essential deletions and essential heterozygous
% mutants. Eliminate the essential genes

load essential_genes_100908;
[t,ind1,ind2] = intersect(xie_huang_2005.orfs, essential_genes);
xie_huang_2005.orfs(ind1) = [];
xie_huang_2005.data(ind1,:) = [];

save('Datasets/Phenotypes/2005_Xie~Huang/xie_huang_2005.mat','xie_huang_2005');

% Save data into database
dt = xie_huang_2005;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    if isempty(datasets(i).reporter)
        datasets_names{i,2} = '';
    else
        datasets_names{i,2} = datasets(i).reporter;
    end
    datasets_names{i,3} = datasets(i).conditionset;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));


%% Botet~Santos, 2007
% DATA = botet_santos_2007
% NOTES = data unnormalized; potentially batch, row/col normalization needed?


botet_santos_2007.source = {'Javier Botet'};
botet_santos_2007.downloaddate = {'2013-04-18'};
botet_santos_2007.pmid = 17873082;

phenotypes = {'Growth, OD'};
treatments = {'SMM, 77 h'; 'SMM, 120 h'; 'Sulfanilamide, 0.1 mg/ml, 77 h';'Sulfanilamide, 0.1 mg/ml, 120 h'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2007_Botet~Santos/1ScreenSULFA&MS&MS+PABA.xlsx','DATA');

% Get indices of the data columns
ind_data = 15:18;

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,10)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,10)),strmatch('Y', data.raw(:,10)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,10) = cellfun(@strtrim, data.raw(:,10),'UniformOutput',0);

data2.orfs = upper(data.raw(:,10));
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Normalize by UNTREATED
data2.data(:,3) = data2.data(:,3)./data2.data(:,1);
data2.data(:,4) = data2.data(:,4)./data2.data(:,2);
data2.data(:,1:2) = [];
treatments(1:2) = [];

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
botet_santos_2007.orfs = t;
botet_santos_2007.data = t2;
botet_santos_2007.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2007_Botet~Santos/botet_santos_2007.mat','botet_santos_2007');

% Save data into database
dt = botet_santos_2007;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    if isempty(datasets(i).reporter)
        datasets_names{i,2} = '';
    else
        datasets_names{i,2} = datasets(i).reporter;
    end
    datasets_names{i,3} = datasets(i).conditionset;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Galvan Marquez~Smith, 2013
% DATA = galvan_marquez_smith_2013


galvan_marquez_smith_2013.source = {'Imelda Galvan Marquez'};
galvan_marquez_smith_2013.downloaddate = {'2013-05-19'};
galvan_marquez_smith_2013.pmid = 23624539;

phenotypes = {'Growth, colony size'};
treatments = {'Chitosan'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2013_GalvanMarquez~Smith/Chitosan Effect on GDA (raw data). Exp. 1 to 3. Imelda Galvan, 2013-1.xlsx','Sheet1');

% Get indices of the data columns
ind_data = find(strcmp('Ratio', data.raw(1,:)));

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,2)),strmatch('Y', data.raw(:,2)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,3) = cellfun(@strtrim, data.raw(:,2),'UniformOutput',0);

% Make sure all the data are numbers
t = data.raw(:,ind_data);
inds = find(cellfun(@isnumeric, t)==0);
t(inds) = {NaN};

data2.orfs = upper(data.raw(:,2));
data2.data = nanmean(cell2mat(t),2);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
galvan_marquez_smith_2013.orfs = t;
galvan_marquez_smith_2013.data = t2;
galvan_marquez_smith_2013.ph = strcat(phenotypes, '; ', treatments);


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(galvan_marquez_smith_2013.orfs, essential_genes);
galvan_marquez_smith_2013.orfs(ind1) = [];
galvan_marquez_smith_2013.data(ind1,:) = [];

save('Datasets/Phenotypes/2013_GalvanMarquez~Smith/galvan_marquez_smith_2013.mat','galvan_marquez_smith_2013');

% Save data into database
dt = galvan_marquez_smith_2013;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    if isempty(datasets(i).reporter)
        datasets_names{i,2} = '';
    else
        datasets_names{i,2} = datasets(i).reporter;
    end
    datasets_names{i,3} = datasets(i).conditionset;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Cai~Becker, 2006
% DATA = cai_becker_2006

cai_becker_2006.source = {'Jeffrey M. Becker'};
cai_becker_2006.downloaddate = {'2013-03-16'};
cai_becker_2006.pmid = 16361226;

phenotypes = {'Growth, OD'};

flds = dir('Datasets/Phenotypes/2006_Cai~Becker/Sceeningdeletionmutant_rawdata/');
flds(1:3) = [];

% PLATE MAPS: Get the data from the Excel files
map = [];
for i = 1 : length(flds)
    tmp = regexp(flds(i).name, '_', 'split');
    plate = str2num(tmp{1});
    
    files = dir(['Datasets/Phenotypes/2006_Cai~Becker/Sceeningdeletionmutant_rawdata/' flds(i).name '/']);
    f = find(cellfun(@isempty, strfind({files.name},'all.xls'))==0);
    filename = ['Datasets/Phenotypes/2006_Cai~Becker/Sceeningdeletionmutant_rawdata/' flds(i).name '/' files(f).name];
    [status,sheets] = xlsfinfo(filename);
    
    for sht = 1 : 3
        [mapdata.num,mapdata.txt,mapdata.raw] = xlsread(filename,sheets{sht});
        
        if length(find(cellfun(@isnumeric,mapdata.raw(:,4)))) < 10
            map = [map; mapdata.raw(:,1:5)];
            break;
        end
    end
end

inds = find(cellfun(@ischar, map(:,1))==0);
map(inds,1) = {'EMPTY'};

map(:,1) = upper(map(:,1));

% PLATE MAPS: A few rows to eliminate
original = {'STRAIN','FILTER'};
for i = 1 : length(original)
    inds = find(strncmp(original{i}, map(:,1), length(original{i})));
    map(inds,:) = [];
end
inds = find(cellfun(@isnumeric, map(:,4))); % Rows with "NaN" coordinates
map(inds,:) = [];
inds = find(strcmp(' ', map(:,4))); % Rows with " " (empty) coordinates
map(inds,:) = [];



% PLATE MAPS: A few rows to rename
original = {'NONE','PTR1','PTR2','CUP9'};
new = {'EMPTY','YGR184C','YKR093W','YPL177C'};
for i = 1 : length(original)
    inds = find(strncmp(original{i}, map(:,1), length(original{i})));
    map(inds,1) = new(i);
end

% PLATE MAPS: Change row letters into numbers
letters = 'ABCDEFGH';
for i = 1 : size(map, 1)
    map{i,4} = strfind(letters, map{i,4});
end


% PLATE MAPS: Keep ORFs and background as cell array, while plate, row &
% columns as numbers
map_txt = map(:,1:2);
map_num = cell2mat(map(:,3:5));


% TIMEPOINTS (1): TO DO ONLY ONCE. Since this information is not standardized, it is necessary to generate a SUMMARY FILE & process it manually
fid = fopen('Datasets/Phenotypes/2006_Cai~Becker/summary.txt','w');
for i = 1 : length(flds)
 
    tmp = regexp(flds(i).name, '_', 'split');
    plate = str2num(tmp{1});
    if length(replicates) < plate
        replicates(plate) = 1;
    else
        replicates(plate) = replicates(plates) + 1;
    end
    
    files = dir(['Datasets/Phenotypes/2006_Cai~Becker/Sceeningdeletionmutant_rawdata/' flds(i).name '/']);
    
    f = find(cellfun(@isempty, strfind({files.name},'.dat'))==0);
    
    for j = 1 : length(f)
        
        filename = ['Datasets/Phenotypes/2006_Cai~Becker/Sceeningdeletionmutant_rawdata/' flds(i).name '/' files(f(j)).name];
        
        C = textread(filename, '%s','delimiter', '\n');

        if length(C) >= 9
            header = C(1);
            fprintf(fid,'%s\t%s\t%d\t%d\t%s\n', flds(i).name, files(f(j)).name, plate, replicates(plate_red), header{h});
        end
        
        if length(C) >= 19
            header = C(11);
            fprintf(fid,'%s\t%s\t%d\t%d\t%s\n', flds(i).name, files(f(j)).name, plate, replicates(plate_red), header{h});
        end
        
        if length(C) >= 29
            header = C(21);
            fprintf(fid,'%s\t%s\t%d\t%d\t%s\n', flds(i).name, files(f(j)).name, plate, replicates(plate_red), header{h});
        end
        
        clear C; 
    end
    
end
fclose(fid);

% TIMEPOINTS (2): The summary.txt file was processed manually to extract information about 1) Medium, 2) Date, 3) Timepoint.

% TIMEPOINTS (3): Load in the new summary2.txt file
C = importdata('Datasets/Phenotypes/2006_Cai~Becker/summary2.txt','\t');
fields = regexp(C{1}, '\t', 'split');
fields(end) = {'TMP'};

for i = 2 : length(C)
    tmp = regexp(C{i}, '\t', 'split');
    for j = 1 : length(fields)
        summary.(fields{j})(i,1) = tmp(j);
    end
end

% DATA (1): Go over all the data and save it properly

flds = dir('Datasets/Phenotypes/2006_Cai~Becker/Sceeningdeletionmutant_rawdata/');
flds(1:3) = [];

data = struct('plateid',[0], ...
    'medium',{'x'},...
    'plate',[0],...
    'date',{'x'},...
    'hour',{'x'},...
    'timepoint',[0],...
    'replicates',[0],...
    'orfs',{'x'},...
    'rows',[0],...
    'cols',[0],...
    'data',[0],...
    'folder',{'x'},...
    'file',{'x'},...
    'header',{'x'});

plateid = 0;    %unique plate identifier

for i = 1 : length(flds)
    
    % To identify a set of data you need folder, file and header.

    inds1 = find(strcmp(flds(i).name, summary.Folder)); % FOLDER
   
    tmp = regexp(flds(i).name, '_', 'split');
    plate = str2num(tmp{1});
    
    % Find the ORFs
    inds = find(map_num(:,1) == plate);
    orfs = map_txt(inds,1);
    rows = map_num(inds,2);
    cols = map_num(inds,3);
    
    inds_data = sub2ind([8 12], rows, cols);
 
    files = dir(['Datasets/Phenotypes/2006_Cai~Becker/Sceeningdeletionmutant_rawdata/' flds(i).name '/']);
    
    f = find(cellfun(@isempty, strfind({files.name},'.dat'))==0);
    
    for j = 1 : length(f)
        
        inds2 = find(strcmp(files(f(j)).name, summary.File)); % FILE
        
        filename = ['Datasets/Phenotypes/2006_Cai~Becker/Sceeningdeletionmutant_rawdata/' flds(i).name '/' files(f(j)).name];
        
        C = textread(filename, '%s','delimiter', '\n');
                
        data_row_num = [1 11 21];
        
        for d = 1 : 3
            
            if length(C) > data_row_num(d)
                
                header = C(data_row_num(d));
                inds3 = find(strcmp(header{1}, summary.Header)); % HEADER
                inds = intersect(intersect(inds1,inds2),inds3); % IDENTIFY DATA
            
                A = dlmread(filename,'\t',[data_row_num(d) 0 data_row_num(d)+7 11]);
            
                if ~isempty(summary.Date{inds}) & ~strcmp('#VALUE!', summary.Hour{inds})
                    
                    plateid = plateid+1;

                    data.plateid = [data.plateid; zeros(length(orfs),1)+plateid];
                    data.medium = [data.medium; repmat(summary.Medium(inds),length(orfs),1)];
                    data.plate = [data.plate; repmat(plate, length(orfs),1)];
                    data.date = [data.date; repmat(summary.Date(inds), length(orfs),1)];
                    data.hour = [data.hour; repmat(summary.Hour(inds), length(orfs),1)];
                    data.timepoint = [data.timepoint; repmat(summary.Timepoint(inds), length(orfs),1)];
                    data.replicates = [data.replicates; repmat(summary.Replicates(inds), length(orfs),1)];
                    data.orfs = [data.orfs; orfs];
                    data.rows = [data.rows; rows];
                    data.cols = [data.cols; cols];
                    data.data = [data.data; A(inds_data)];
                    data.folder = [data.folder; repmat(summary.Folder(inds), length(orfs),1)];
                    data.file = [data.file; repmat(summary.File(inds), length(orfs),1)];
                    data.header = [data.header; repmat(header, length(orfs),1)];
                    
                end
            end
        end
    end
end


fields = fieldnames(data);
for i = 1 : length(fields)
    data.(fields{i})(1) = [];   % delete the first initiating value.
end

% DATA (2): Transform a few values into numbers
data.timepoint = cellfun(@str2num, data.timepoint); % timepoints
data.replicates = cellfun(@str2num, data.replicates); % replicates
[tmp,ia,data.medium_num] = unique(data.medium); % medium

[tmp_u,ia,ib] = unique(strcat(data.date, {' '}, data.hour)); % date & time
tmp_u_vc = datevec(tmp_u);
tmp_u_vc(tmp_u_vc(:,2)>1,1) = 2003; tmp_u_vc(tmp_u_vc(:,2)==1,1) == 2004;
tmp_u_nm = datenum(tmp_u_vc);
data.datenum = tmp_u_nm(ib);

% DATA (3): Calculate proper time intervals between timepoints
plate_rep_medium = unique([data.plate data.replicates data.medium_num],'rows');
for i = 1 : size(plate_rep_medium,1)
    
    inds = find(data.plate == plate_rep_medium(i,1) & data.replicates == plate_rep_medium(i,2) & data.medium_num == plate_rep_medium(i,3));
    [timepoints,ia,ib] = unique(data.timepoint(inds));
    datetime = data.datenum(inds(ia));
    timeintervals = [0; diff(datetime)] *24;    % number of hours between consecutive timepoints
    
    data.timeintervals(inds) = timeintervals(ib);
        
end

         
% % QC: Visualize the data to see if there're any biases.
% inds = find(data.plateid == randsample(unique(data.plateid),1)); % grab a random plate
% otherplates = find(data.plate == data.plate(inds(1)) & data.replicates == data.replicates(inds(1)) & strcmp(data.medium{inds(1)}, data.medium)); % find the other timepoints for this plate, this medium and this replicate
% 
% plateids = unique(data.plateid(otherplates));
% figure();
% for i = 1 : length(plateids)
%     subplot(3,3,i);
%     plate_map = zeros(8,12);
%     inds = find(data.plateid == plateids(i));
%     data_inds = sub2ind([8 12], data.rows(inds), data.cols(inds));
%     plate_map(data_inds) = data.data(inds);
% 
%     imagesc(plate_map,[0 0.3]);
%     title([data.medium{inds(1)} ', Plate ' num2str(data.plate(inds(1))) ', Replicate ' num2str(data.replicates(inds(1))) ', Date/time ' data.date{inds(1)} '/' data.hour{inds(1)}]);
%     colorbar;
%     axis image;
% end
                
% DATA (4): Combine all timepoints into a single measure. Arbitrary = Area
% under the growth curve

% Prepare summary matrix ORFS x Medium x Replicates
all_orfs = unique(data.orfs);
all_media = unique(data.medium);

plate_rep_medium = unique([data.plate data.replicates data.medium_num],'rows');
data_matrix = zeros(length(all_orfs), length(all_media),size(plate_rep_medium,1))+NaN;
curr_rep = 1;
for i = 1 : size(plate_rep_medium,1)
    inds = find(data.plate == plate_rep_medium(i,1) & data.replicates == plate_rep_medium(i,2) & data.medium_num == plate_rep_medium(i,3));
    [timeintervals,ia,ib] = unique(data.timeintervals(inds));
    [orfs,i1,i2] = unique(data.orfs(inds));
    auc = zeros(length(orfs),1);
    for o = 1 : length(orfs)
        inds2 = find(strcmp(orfs{o}, data.orfs(inds)));
        [t,ix] = sort(data.timepoint(inds(inds2)));
        dt = data.data(inds(inds2(ix)))-data.data(inds(inds2(ix(1))));  % difference in growth between now and the first timepoint
        auc(o) = nansum(data.timeintervals(inds(inds2(ix))).*dt);
        if isnan(auc(o))
            break;
        end
    end
    % Normalize by the average of the entire plate
    auc_norm = auc ./ nanmean(auc);
    
    [t,ind1,ind2] = intersect(all_orfs, orfs);
    data_matrix(ind1,plate_rep_medium(i,3),curr_rep) = auc_norm(ind2);
    curr_rep = curr_rep + 1;
end
    
% Check how many replicates ORFs have in general
n = sum(~isnan(data_matrix),3);

% Average data for ORFs across all replicates
data_matrix2 = nanmean(data_matrix,3);

cai_becker_2006.orfs = all_orfs;
cai_becker_2006.ph = strcat(phenotypes, '; ', all_media);
cai_becker_2006.data = data_matrix2;

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', cai_becker_2006.orfs,1));
cai_becker_2006.orfs(inds) = [];
cai_becker_2006.data(inds,:) = [];



% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(cai_becker_2006.orfs, essential_genes);
cai_becker_2006.orfs(ind1) = [];
cai_becker_2006.data(ind1,:) = [];

save('Datasets/Phenotypes/2006_Cai~Becker/cai_becker_2006.mat','cai_becker_2006');

% Save data into database
dt = cai_becker_2006;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    if isempty(datasets(i).reporter)
        datasets_names{i,2} = '';
    else
        datasets_names{i,2} = datasets(i).reporter;
    end
    datasets_names{i,3} = datasets(i).conditionset;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

%% Chesi~Gitler, 2012
% DATA = chesi_gitler_2012

chesi_gitler_2012.pmid = 22457822;

phenotypes = {'Growth, colony size'};
treatments = {'Manganese'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2012_Chesi~Gitler/yeast deletions Mn.xlsx','single deletion');

crr = zeros(size(data.raw,2)-1,3);  % Concentration, within-round replicate, round
for i = 2 : length(data.raw(1,:))
    t = textscan(data.raw{1,i},'%d %s size-%d %d');
    crr(i-1,:) = cell2mat(t([1 3 4]));   
end

% Get indices of the data columns
ind_data = 2:19;


% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

inds = find(~strncmp('Y', data.raw(:,1),1));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,1) = cellfun(@strtrim, data.raw(:,1),'UniformOutput',0);

% Make sure all the data are numbers
t = data.raw(:,ind_data);
inds = find(~cellfun(@isnumeric, t));
t(inds) = {NaN};

ind_data_conc1 = [1 2 7 8 13 14];
ind_data_conc2 = [3 4 9 10 15 16];
ind_data_conc3 = [5 6 11 12 17 18];

data2.orfs = upper(data.raw(:,1));
data2.data(:,1) = nanmean(cell2mat(t(:,ind_data_conc1)),2);
data2.data(:,2) = nanmean(cell2mat(t(:,ind_data_conc2)),2);
data2.data(:,3) = nanmean(cell2mat(t(:,ind_data_conc3)),2);

% Normalize to no drug
data2.data(:,2) = data2.data(:,2) ./ data2.data(:,1);
data2.data(:,3) = data2.data(:,3) ./ data2.data(:,1);
data2.data(:,1) = [];

concentrations = unique(crr(:,1));
concentrations(1) = [];

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
chesi_gitler_2012.orfs = t;
chesi_gitler_2012.data = t2;
chesi_gitler_2012.ph = strcat(phenotypes, '; ', treatments, '; ', num2str(concentrations), ' mM');

save('Datasets/Phenotypes/2012_Chesi~Gitler/chesi_gitler_2012.mat','chesi_gitler_2012');

% Save data into database
dt = chesi_gitler_2012;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Askree~McEachern, 2004
% DATA = askree_mceachern_2004


askree_mceachern_2004.source = {'http://www.pnas.org/content/suppl/2004/05/20/0401263101.DC1/01263Table3.xls; Michael McEachern'};
askree_mceachern_2004.downloaddate = {'2013-07-17'};
askree_mceachern_2004.pmid = 15161972;

phenotypes = {'Telomere length'};
treatments = {''};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2004_Askree~McEachern/01263Table3.xlsx','Sheet1');

% Get indices of the data columns
ind_data = 3:4;


% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

% Set data for inconsistent mutants to NaN
inds = find(strncmp('*Y', data.raw(:,1),2));
data.raw(inds,ind_data) = {NaN};

% Eliminate the asterisk before the ORF name
temp_orfs = char(data.raw(inds,1));
data.raw(inds,1) = cellstr(temp_orfs(:,2:end));

inds = find(~strncmp('Y', data.raw(:,1),1));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,1) = cellfun(@strtrim, data.raw(:,1),'UniformOutput',0);

% Make sure all the data are numbers
t = data.raw(:,ind_data);
inds = find(~cellfun(@isnumeric, t));
t(inds) = {NaN};

data2.orfs = upper(data.raw(:,1));
data2.data(:,1) = nanmean(cell2mat(t),2);

% Load the tested genes
[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2004_Askree~McEachern/S.cerftpmata.xlsx','Sheet1');
tested_orfs = data.raw(:,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds,:) = [];

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
data3.orfs = t;
data3.data = t2;

% Put everything together
askree_mceachern_2004.orfs = tested_orfs;
askree_mceachern_2004.ph = phenotypes;
askree_mceachern_2004.data = zeros(length(tested_orfs),length(phenotypes));
[t,ind1,ind2] = intersect(data3.orfs, askree_mceachern_2004.orfs);
askree_mceachern_2004.data(ind2,:) = data3.data(ind1,:);



% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(askree_mceachern_2004.orfs, essential_genes);
askree_mceachern_2004.orfs(ind1) = [];
askree_mceachern_2004.data(ind1,:) = [];

save('Datasets/Phenotypes/2004_Askree~McEachern/askree_mceachern_2004.mat','askree_mceachern_2004');

% Save data into database
dt = askree_mceachern_2004;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Brett~Rao, 2011
% DATA = brett_rao_2011

brett_rao_2011.source = {'http://www.plosone.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pone.0017619.s003'};
brett_rao_2011.downloaddate = {'2013-07-18'};
brett_rao_2011.pmid = 21423800;

phenotypes = {'Growth, OD600';'Vacuolar pH (pHv)'};
treatments = {'external pH 2.7';'external pH 4.0';'external pH 7.0'};

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2011_Brett~Rao/journal.pone.0017619.s003.xlsx','Unsorted Data');

% Get indices of the data columns
ind_data = [4:6 7:9];


% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

inds = find(~strncmp('Y', data.raw(:,2),1));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,2) = cellfun(@strtrim, data.raw(:,2),'UniformOutput',0);

% Make sure all the data are numbers
t = data.raw(:,ind_data);
inds = find(~cellfun(@isnumeric, t));
t(inds) = {NaN};

data2.orfs = upper(data.raw(:,2));
data2.data = cell2mat(t);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
brett_rao_2011.orfs = t;
brett_rao_2011.data = t2;
brett_rao_2011.ph = [strcat(phenotypes{1}, '; ', treatments); strcat(phenotypes{2},'; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(brett_rao_2011.orfs, essential_genes);
brett_rao_2011.orfs(ind1) = [];
brett_rao_2011.data(ind1,:) = [];

save('Datasets/Phenotypes/2011_Brett~Rao/brett_rao_2011.mat','brett_rao_2011');

% Save data into database
dt = brett_rao_2011;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Rieger~Maundrell, 1999
% DATA = rieger_maundrell_1999

rieger_maundrell_1999.source = {'manuscript PDF'};
rieger_maundrell_1999.downloaddate = {'2014-01-31'};
rieger_maundrell_1999.pmid = 10407277;

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/1999_Rieger~Maundrell/rieger_maundrell_1999_data.xlsx','Data');

phenotypes = {'growth'};
treatments = data.txt(2:end,1);
phenotype_severity = {'HS','S','R'};
phenotype_severity_num = [-2,-1,1];

% Eliminate white spaces before/after ORF
data.raw(1,2:end) = cellfun(@strtrim, data.raw(1,2:end),'UniformOutput',0);

% Replace string scores with numbers
t = data.raw(2:end,2:end);
for i = 1 : length(phenotype_severity)
    inds = find(strcmp(phenotype_severity{i}, t));
    t(inds) = {phenotype_severity_num(i)};
end

rieger_maundrell_1999.orfs = upper(data.raw(1,2:end))';
rieger_maundrell_1999.data = cell2mat(t');
rieger_maundrell_1999.data(isnan(rieger_maundrell_1999.data)) = 0;
rieger_maundrell_1999.ph = [strcat(phenotypes{1}, '; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(rieger_maundrell_1999.orfs, essential_genes);
rieger_maundrell_1999.orfs(ind1) = [];
rieger_maundrell_1999.data(ind1,:) = [];

save('Datasets/Phenotypes/1999_Rieger~Maundrell/rieger_maundrell_1999.mat','rieger_maundrell_1999');

% Save data into database
dt = rieger_maundrell_1999;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
adj_ix = [1 3 4 7 8 9 15 21 19 30 34 35 37 39 40 10 45 43 46 55 59 67 11 12 13 14 16 17 18 20 22 23 24 25 26 27 28 29 31 32 33 36 38 41 42 44 47 48 49 50 51 52 53 54 56 57 58 2 5 6 60 61 62 63 64 65 66];
[~,adj_ix] = sort(adj_ix);

datasets.names(database_ix(adj_ix),:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix(adj_ix)));

%% Chan~Zheng, 2000
% DATA = chan_zheng_2000
% NEED = tested genes

chan_zheng_2000.source = {'Supplementary Table 3'};
chan_zheng_2000.downloaddate = {'2014-01-31'};
chan_zheng_2000.pmid = 11078525;

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2000_Chan~Zheng/chan_zheng_2000_HAP.xlsx','Sheet1');

phenotypes = {'growth'};
treatments = {'rapamycin, 25 nM'};

% Eliminate white spaces before/after ORF
data.raw(:,1) = cellfun(@strtrim, data.raw(:,1),'UniformOutput',0);

chan_zheng_2000.orfs = upper(data.raw(:,1));
chan_zheng_2000.data = cell2mat(data.raw(:,2));
chan_zheng_2000.ph = [strcat(phenotypes, '; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(chan_zheng_2000.orfs, essential_genes);
chan_zheng_2000.orfs(ind1) = [];
chan_zheng_2000.data(ind1,:) = [];

save('Datasets/Phenotypes/2000_Chan~Zheng/chan_zheng_2000.mat','chan_zheng_2000');

% Save data into database
dt = chan_zheng_2000;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Ooi~Boeke, 2001
% DATA = ooi_boeke_2001
% NEED = tested genes

ooi_boeke_2001.source = {'Main PDF'};
ooi_boeke_2001.downloaddate = {'2014-02-03'};
ooi_boeke_2001.pmid = 11701889;

hits = textread('Datasets/Phenotypes/2001_Ooi~Boeke/ooi_boeke_2001.txt','%s');
hits = lower(hits);
hits(strcmp('lig4', hits)) = {'dnl4'};
hits(strcmp('gpe2', hits)) = {'yal056w'};

phenotypes = {'NHEJ defect'};
treatments = {''};

% Transform gene names into ORFs
hits_orf = genename2orf(hits);
tmp = split_by_delimiter('_',hits_orf);
hits_orf = tmp(:,1);

ooi_boeke_2001.orfs = upper(hits_orf);
ooi_boeke_2001.data = -ones(size(hits_orf));
ooi_boeke_2001.ph = [strcat(phenotypes, '; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(ooi_boeke_2001.orfs, essential_genes);
ooi_boeke_2001.orfs(ind1) = [];
ooi_boeke_2001.data(ind1,:) = [];

save('Datasets/Phenotypes/2001_Ooi~Boeke/ooi_boeke_2001.mat','ooi_boeke_2001');

% % Save data into database
% dt = ooi_boeke_2001;
% datasets = get_datasets_for_paper(dt);
% 
% [~,database_ix] = sortrows(datasets.names,[1 2 3]);
% [~,ph_ix] = sort(dt.ph);
% 
% % Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
% datasets.names(database_ix,:)
% dt.ph(ph_ix)
% 
% insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Bianchi~Frontali, 2001
% DATA = bianchi_frontali_2001

bianchi_frontali_2001.pmid = 11746602;

phenotypes = {'Growth'};

[treatments_name, treatments_dose, treatments_abbr] = ...
    textread('Datasets/Phenotypes/2001_Bianchi~Frontali/bianchi_frontali_2001_conditions.txt','%s %s %s','delimiter','\t');
tested_genes = textread('Datasets/Phenotypes/2001_Bianchi~Frontali/bianchi_frontali_2001_tested.txt','%s');
[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2001_Bianchi~Frontali/bianchi_frontali_2001_data.xlsx','Sheet1');

phenotype_severity = {'HS','S','R'};
phenotype_severity_num = [-2, -1, 1];

bianchi_frontali_2001.orfs = upper(tested_genes);
bianchi_frontali_2001.orfs = cellfun(@strtrim, bianchi_frontali_2001.orfs,'UniformOutput',0);

bianchi_frontali_2001.ph = strcat(phenotypes, '; ', treatments_name, ', ', treatments_dose);

bianchi_frontali_2001.data = zeros(length(bianchi_frontali_2001.orfs), length(bianchi_frontali_2001.ph));

inds = find(cellfun(@isnumeric, data.raw));
data.raw(inds) = {''};

for i = 1 : size(data.raw,1)
    
    inds_orf = strcmp(upper(data.raw{i,1}), bianchi_frontali_2001.orfs);
       
    for j = 2 : 4
        tmp = regexp(data.raw{i,j},', ', 'split');
        for k = 1 : length(tmp)
            if ~isempty(tmp{k})
                inds_treat = find(strcmp(tmp{k}, treatments_abbr));
                bianchi_frontali_2001.data(inds_orf,inds_treat) = phenotype_severity_num(j-1);
            end
        end
    end
    
end

% Check
n = sum(abs(bianchi_frontali_2001.data)>0,2);
length(find(n>0))   % should be 163

n = sum(abs(bianchi_frontali_2001.data)>0,1);
length(find(n==0))  % should be 0

save('Datasets/Phenotypes/2001_Bianchi~Frontali/bianchi_frontali_2001.mat','bianchi_frontali_2001');

% Save data into database
dt = bianchi_frontali_2001;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);
ph_ix = ph_ix([7 5 6 9 20 18 22 21 10 11 28 27 12 32 1 33 2 13 14 15 16 17 19 23 24 25 26 29 30 31 34 3 8 4 35 36 37]);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
adj_ix = [7 5 6 9 20 18 22 21 10 11 12 32 1 33 2 13 14 15 16 17 19 23 24 25 26 28 27 29 30 31 34 3 8 4 35 36 37];
[~,adj_ix] = sort(adj_ix);

datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix(adj_ix)));

%% Akache~Turcotte, 2002
% DATA = akache_turcotte_2002

akache_turcotte_2002.source = {'manuscript PDF'};
akache_turcotte_2002.downloaddate = {'2014-02-10'};
akache_turcotte_2002.pmid = 11943786;

[data.num,data.txt,data.raw] = xlsread('Datasets/Phenotypes/2002_Akache~Turcotte/akache_turcotte_2002_data.xlsx','Sheet1');

phenotypes = {'growth'};
treatments = data.txt(1,2:end)';
phenotype_severity = {'S';'SS';'R';'RR'};
phenotype_severity_num = [-1,-2,1,2];

% Eliminate white spaces before/after ORF
data.raw(2:end,1) = cellfun(@strtrim, data.raw(2:end,1),'UniformOutput',0);

% Replace string scores with numbers
t = data.raw(2:end,2:end);
for i = 1 : length(phenotype_severity)
    inds = find(strcmp(phenotype_severity{i}, t));
    t(inds) = {phenotype_severity_num(i)};
end
t(cellfun(@isnan,t))={0};

akache_turcotte_2002.orfs = upper(data.raw(2:end,1));
akache_turcotte_2002.data = cell2mat(t);
akache_turcotte_2002.ph = [strcat(phenotypes{1}, '; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(akache_turcotte_2002.orfs, essential_genes);
akache_turcotte_2002.orfs(ind1) = [];
akache_turcotte_2002.data(ind1,:) = [];

save('Datasets/Phenotypes/2002_Akache~Turcotte/akache_turcotte_2002.mat','akache_turcotte_2002');

% Save data into database
dt = akache_turcotte_2002;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));
 
%% Desmoucelles~Daignan-Fornier, 2002
% DATA = desmoucelles_daignan_fornier_2002
% NOTE = 1) typo in data: YKR065W should be YKR065C 2) YCR002C is in the
% result, but not in tested.

desmoucelles_daignan_fornier_2002.source = {'manuscript PDF'};
desmoucelles_daignan_fornier_2002.downloaddate = {'2014-02-10'};
desmoucelles_daignan_fornier_2002.pmid = 12016207;

[data.num,data.txt,data.raw] = ...
    xlsread('Datasets/Phenotypes/2002_Desmoucelles~Daignan-Fornier/desmoucelles_daignan_fornier_2002_data.xlsx','data.txt');
[tested.num,tested.txt,tested.raw] = ...
    xlsread('Datasets/Phenotypes/2002_Desmoucelles~Daignan-Fornier/euroscarf list.xlsx','1_1.xlwb');

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, tested.raw(:,2)));
tested.raw(inds,:) = [];

inds = find(~strncmp('Y', tested.raw(:,2),1));
tested.raw(inds,:) = [];

tested.raw(:,2) = cellfun(@strtrim, tested.raw(:,2),'UniformOutput',0);

phenotypes = {'growth'};
treatments = {'MPA'};
phenotype_severity = {'S','SS','R','RR'};
phenotype_severity_num = [-1,-2,1,2];

% Eliminate white spaces before/after ORF
data.raw(1:end,1) = cellfun(@strtrim, data.raw(1:end,1),'UniformOutput',0);
inds = find(strcmp('YKR065W', data.raw(:,1)));
data.raw(inds,1) = {'YKR065C'};

% Replace string scores with numbers
t = data.raw(1:end,2);
for i = 1 : length(phenotype_severity)
    inds = find(strcmp(phenotype_severity{i}, t));
    t(inds) = {phenotype_severity_num(i)};
end
t(cellfun(@isnan,t))={0};

desmoucelles_daignan_fornier_2002.orfs = upper(tested.raw(:,2));
desmoucelles_daignan_fornier_2002.data = zeros(length(desmoucelles_daignan_fornier_2002.orfs),1);

[tmp,ind1,ind2] = intersect(desmoucelles_daignan_fornier_2002.orfs, upper(data.raw(1:end,1)));
desmoucelles_daignan_fornier_2002.data(ind1) = cell2mat(t(ind2));

desmoucelles_daignan_fornier_2002.ph = [strcat(phenotypes{1}, '; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(desmoucelles_daignan_fornier_2002.orfs, essential_genes);
desmoucelles_daignan_fornier_2002.orfs(ind1) = [];
desmoucelles_daignan_fornier_2002.data(ind1,:) = [];

save('Datasets/Phenotypes/2002_Desmoucelles~Daignan-Fornier/desmoucelles_daignan_fornier_2002.mat','desmoucelles_daignan_fornier_2002'); 

% Save data into database
dt = desmoucelles_daignan_fornier_2002;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Blackburn~Avery, 2003
% DATA = blackburn_avery_2003
% TESTED = not available

blackburn_avery_2003.source = {'manuscript PDF'};
blackburn_avery_2003.downloaddate = {'2014-02-12'};
blackburn_avery_2003.pmid = 12543677;

[data.num,data.txt,data.raw] = ...
    xlsread('Datasets/Phenotypes/2003_Blackburn~Avery/blackburn_avery_2003_data.xlsx','data.txt');

phenotypes = {'growth (MIC)'};
treatments = data.raw(1,2:8)';

% Eliminate white spaces before/after ORF
data.raw(2:end,1) = cellfun(@strtrim, data.raw(2:end,1),'UniformOutput',0);

% Replace 'Inf' with Inf
t = data.raw(2:end,2:end);
t(~cellfun(@isnumeric, t)) = {Inf};

blackburn_avery_2003.orfs = upper(data.raw(2:end,1));
blackburn_avery_2003.data = cell2mat(t);

blackburn_avery_2003.ph = [strcat(phenotypes{1}, '; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(blackburn_avery_2003.orfs, essential_genes);
blackburn_avery_2003.orfs(ind1) = [];
blackburn_avery_2003.data(ind1,:) = [];

save('Datasets/Phenotypes/2003_Blackburn~Avery/blackburn_avery_2003.mat','blackburn_avery_2003'); 

% Save data into database
dt = blackburn_avery_2003;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Parsons~Boone, 2004
% DATA = parsons_boone_2004
% TESTED = not available

parsons_boone_2004.source = {'http://www.nature.com/nbt/journal/v22/n1/extref/nbt919-S2.xls'};
parsons_boone_2004.downloaddate = {'2014-02-20'};
parsons_boone_2004.pmid = 14661025;

[data.num,data.txt,data.raw] = ...
    xlsread('Datasets/Phenotypes/2004_Parsons~Boone/nbt919-S2.xlsx','Sheet1');

inds = find(strcmp('ORF', data.raw(:,1)));

phenotypes = {'growth (colony size)'};
treatments = data.raw(inds,3:end)';

% Eliminate white spaces before/after ORF
data.raw(inds+1:end,1) = cellfun(@strtrim, data.raw(inds+1:end,1),'UniformOutput',0);

% Replace 'Inf' with Inf
t = data.raw(inds+1:end,3:end);
t(cellfun(@isnan, t)) = {0};

parsons_boone_2004.orfs = upper(data.raw(inds+1:end,1));
parsons_boone_2004.data = cell2mat(t);

% Flip the sign of the values, such that negative = sensitive, positive =
% resistant
parsons_boone_2004.data = -parsons_boone_2004.data;

parsons_boone_2004.ph = [strcat(phenotypes{1}, '; ', treatments)];


% Eliminate the essential genes
% load essential_genes_100908;
% [t,ind1,ind2] = intersect(parsons_boone_2004.orfs, essential_genes);
% parsons_boone_2004.orfs(ind1) = [];
% parsons_boone_2004.data(ind1,:) = [];

save('Datasets/Phenotypes/2004_Parsons~Boone/parsons_boone_2004.mat','parsons_boone_2004'); 

% Save data into database
dt = parsons_boone_2004;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
adj_ix = [4 6 8 1 2 3 5 7 9 10 11 12];
[~,adj_ix] = sort(adj_ix);

datasets.names(database_ix(adj_ix),:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix(adj_ix)));

%% Serrano~Arino, 2004
% DATA = serrano_arino_2004

serrano_arino_2004.source = {'main PDF'};
serrano_arino_2004.downloaddate = {'2014-03-03'};
serrano_arino_2004.pmid = 14993228;

[data.num,data.txt,data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2004_Serrano~Arino/serrano_arino_2004.xlsx','Sheet1');

phenotypes = {'growth (colony size)'};
treatments = {'pH 6-2.-7.5'};

% Eliminate white spaces before/after ORF
data.raw(:,1) = cellfun(@strtrim, data.raw(:,1),'UniformOutput',0);

% Translate genenames to ORF
data.raw(:,4) = genename2orf(data.raw(:,1),'noannot');

% Manually adjust the genenames that couldn't be matched
data.raw(strcmpi('cwh36', data.raw(:,4)),4) = {'YCL007C'};
data.raw(strcmpi('lys7', data.raw(:,4)),4) = {'YMR038C'};
data.raw(strcmpi('rcs1', data.raw(:,4)),4) = {'YGL071W'};

hits_orfs = data.raw(:,4);
scores = cell2mat(data.raw(:,2));

% Flip the scores such that -5 is the most sensitive and -1 is the least
% sensitive
scores = scores - 6;

% Load tested genes
[tested.num,tested.txt,tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2004_Serrano~Arino/BY4741.xlsx','Tabelle1');
tested_orfs = unique(upper(tested.raw(2:end,2)));

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
serrano_arino_2004.orfs = tested_orfs;
serrano_arino_2004.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(serrano_arino_2004.orfs, hits_orfs);
serrano_arino_2004.data(ind1,:) = scores(ind2,:);

serrano_arino_2004.ph = [strcat(phenotypes{1}, '; ', treatments)];


save('Datasets/Phenotypes/2004_Serrano~Arino/serrano_arino_2004.mat','serrano_arino_2004'); 

% Save data into database
dt = serrano_arino_2004;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Giorgini~Muchowski, 2005
% DATA = giorgini_muchowski_2005

giorgini_muchowski_2005.pmid = 15806102;

hits = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Giorgini~Muchowski/giorgini_muchowski_2005_hits.txt','%s');

phenotypes = {'growth (pooled CFU)'};
treatments = {'Htt103Q'};

% Eliminate white spaces before/after ORF
hits(:,1) = cellfun(@strtrim, hits,'UniformOutput',0);

% Translate genenames to ORF
hits_orfs = genename2orf(hits,'noannot');

% Manually adjust the genenames that couldn't be matched
hits_orfs(strcmpi('arg7', hits_orfs)) = {'YMR062C'};
hits_orfs(strcmpi('YBR089C-A', hits_orfs)) = {'YBR090C-A'}; % the current version YBR089C-A is not in the tested_orfs list

% Hits = LOF suppressors of Htt103Q sensitivity, so positive score.
scores = ones(length(hits_orfs),1);

% Load tested genes
[tested.num,tested.txt,tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Giorgini~Muchowski/Mat_a_obs_v2(1).0.xlsx','DATA');
inds = find(cellfun(@isnumeric, tested.raw(:,2)));
tested.raw(inds,:) = [];
tested_orfs = unique(upper(tested.raw(2:end,2)));

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
giorgini_muchowski_2005.orfs = tested_orfs;
giorgini_muchowski_2005.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(giorgini_muchowski_2005.orfs, hits_orfs);
giorgini_muchowski_2005.data(ind1,:) = scores(ind2,:);

giorgini_muchowski_2005.ph = [strcat(phenotypes{1}, '; ', treatments)];


save('Datasets/Phenotypes/2005_Giorgini~Muchowski/giorgini_muchowski_2005.mat','giorgini_muchowski_2005'); 

% Save data into database
dt = giorgini_muchowski_2005;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Hellauer~Turcotte, 2005
% DATA = hellauer_turcotte_2005

hellauer_turcotte_2005.source = {'main PDF'};
hellauer_turcotte_2005.downloaddate = {'2014-03-04'};
hellauer_turcotte_2005.pmid = 16061773;

[hits_orfs, scores] = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Hellauer~Turcotte/hits_all.txt','%s %d');

phenotypes = {'growth (colony size)'};
treatments = {'tirapazamine'};

% Eliminate white spaces before/after ORF
hits_orfs(:,1) = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

% Load tested genes
tested_orfs = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Hellauer~Turcotte/ORF.txt','%s');

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
hellauer_turcotte_2005.orfs = tested_orfs;
hellauer_turcotte_2005.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(hellauer_turcotte_2005.orfs, hits_orfs);
hellauer_turcotte_2005.data(ind1,:) = scores(ind2,:);

hellauer_turcotte_2005.ph = [strcat(phenotypes{1}, '; ', treatments)];


save('Datasets/Phenotypes/2005_Hellauer~Turcotte/hellauer_turcotte_2005.mat','hellauer_turcotte_2005'); 

% Save data into database
dt = hellauer_turcotte_2005;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

% %% Proszynski~Walch-Solimena, 2005
% % DATA = proszynski_walch_solimena_2005
% % TESTED = waiting
% 
% proszynski_walch_solimena_2005.source = {'http://www.pnas.org/content/suppl/2005/11/29/0509107102.DC1/09107Table2.xls'};
% proszynski_walch_solimena_2005.downloaddate = {'2014-03-07'};
% proszynski_walch_solimena_2005.pmid = 16330752;
% 
% [data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Proszynski~Walch-Solimena/09107Table2.xlsx','Sheet1');
% 
% % There're 2 classes of mutants, with 2 different patterns of altered
% % mislocalization. As a first pass, I'll treat them the same just as a
% % "mislocalized phenotype".
% phenotypes = {'localization to cell surface'};
% treatments = {''};
% 
% % Find the ORF column
% [r,c] = find(strcmp('ORF', data.raw));
% 
% hits_orfs = data.raw(r+1:end,c);
% inds = find(cellfun(@isnumeric, hits_orfs));
% hits_orfs(inds) = [];
% 
% % Eliminate white spaces before/after ORF
% hits_orfs(:,1) = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

% % Load tested genes
% tested_orfs = ...
%     textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Hellauer~Turcotte/ORF.txt','%s');
% 
% % Eliminate anything that doesn't look like an ORF
% inds = find(~strncmp('Y', tested_orfs,1));
% tested_orfs(inds) = [];
% 
% % Check if all the hits are in the tested space
% [missing,inds] = setdiff(hits_orfs, tested_orfs);
% 
% % Create dataset
% hellauer_turcotte_2005.orfs = tested_orfs;
% hellauer_turcotte_2005.data = zeros(length(tested_orfs), length(treatments));
% [t,ind1,ind2] = intersect(hellauer_turcotte_2005.orfs, hits_orfs);
% hellauer_turcotte_2005.data(ind1,:) = scores(ind2,:);
% 
% hellauer_turcotte_2005.ph = [strcat(phenotypes{1}, '; ', treatments)];
% 
% 
% save('Datasets/Phenotypes/2005_Hellauer~Turcotte/hellauer_turcotte_2005.mat','hellauer_turcotte_2005');


% %% Serviene~Nagy, 2005
% % DATA = serviene_nagy_2005
% % TESTED = waiting
% 
% serviene_nagy_2005.source = {'main PDF'};
% serviene_nagy_2005.downloaddate = {'2014-03-07'};
% serviene_nagy_2005.pmid = 16027361;
% 
% hits1_genenames = {'ctl1', 'met22', 'xrn1', 'ubp3','hur1'};
% hits1_scores = [2,2,2,2,2];
% 
% hits2_genenames = {'pep7','ipk1','cho2','dci1'};
% hits2_scores = [-2,-2,-2,-2];
% 
% hits3_genenames = {'vam7','pth1','vps29','vps35','ngg1'};
% hits3_scores = [-1,-1,-1,-1,-1];
% 
% hits = [hits1_genenames'; hits2_genenames'; hits3_genenames'];
% hits_scores = [hits1_scores'; hits2_scores'; hits3_scores'];
% 
% hits_orfs = genename2orf(hits,'noannot');
% hits_orfs(strcmpi('xrn1', hits_orfs)) = {'YGL173C'};
% 
% phenotypes = {'viral RNA recombination'};
% treatments = {''};
% 
% 
% % Eliminate white spaces before/after ORF
% hits_orfs(:,1) = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

% % Load tested genes
% tested_orfs = ...
%     textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Hellauer~Turcotte/ORF.txt','%s');
% 
% % Eliminate anything that doesn't look like an ORF
% inds = find(~strncmp('Y', tested_orfs,1));
% tested_orfs(inds) = [];
% 
% % Check if all the hits are in the tested space
% [missing,inds] = setdiff(hits_orfs, tested_orfs);
% 
% % Create dataset
% hellauer_turcotte_2005.orfs = tested_orfs;
% hellauer_turcotte_2005.data = zeros(length(tested_orfs), length(treatments));
% [t,ind1,ind2] = intersect(hellauer_turcotte_2005.orfs, hits_orfs);
% hellauer_turcotte_2005.data(ind1,:) = scores(ind2,:);
% 
% hellauer_turcotte_2005.ph = [strcat(phenotypes{1}, '; ', treatments)];
% 
% 
% save('Datasets/Phenotypes/2005_Hellauer~Turcotte/hellauer_turcotte_2005.mat','hellauer_turcotte_2005');

%% Gatbonton~Bedalov, 2006
% DATA = gatbonton_bedalov_2006

gatbonton_bedalov_2006.source = {'main PDF'};
gatbonton_bedalov_2006.downloaddate = {'2014-03-10'};
gatbonton_bedalov_2006.pmid = 16552446;

[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Gatbonton~Bedalov/gatbonton_bedalov_2006_hits.xlsx','Sheet1');

hits_orfs = genename2orf(data.raw(:,1),'noannot');

% Adjust the genenames that couldn't be mapped automatically
hits_orfs(strcmpi('apg17', hits_orfs)) = {'YLR423C'};
hits_orfs(strcmpi('fyv13', hits_orfs)) = {'YGR160W'};
hits_orfs(strcmpi('sig1', hits_orfs)) = {'YER068W'};
hits_orfs(strcmpi('srb10', hits_orfs)) = {'YPL042C'};
hits_orfs(strcmpi('srb9', hits_orfs)) = {'YDR443C'};
hits_orfs(strcmpi('ssn6', hits_orfs)) = {'YBR112C'};
hits_orfs(strcmpi('vps22', hits_orfs)) = {'YPL002C'};
hits_orfs(strcmpi('vps23', hits_orfs)) = {'YCL008C'};

% Eliminate white spaces before/after ORF
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

hits_scores = cell2mat(data.raw(:,3));

% Converting the scores such that: 
% a) 1 = weakest phenotype, 3 = strongest phenotype
% b) long telomeres = positive scores, short telomeres = negative scores
hits_scores = abs(hits_scores-4);
inds = find(strcmp('S', data.raw(:,2)));
hits_scores(inds) = -hits_scores(inds);

phenotypes = {'telomere length'};
treatments = {''};


% Load tested genes
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Gatbonton~Bedalov/genelist_altered_020806.xlsx','mat alpha copy.txt');

tested_orfs = data.raw(2:end,1);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

tested_orfs = unique(tested_orfs);

% Create dataset
gatbonton_bedalov_2006.orfs = tested_orfs;
gatbonton_bedalov_2006.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(gatbonton_bedalov_2006.orfs, hits_orfs);
gatbonton_bedalov_2006.data(ind1,:) = hits_scores(ind2,:);

gatbonton_bedalov_2006.ph = [strcat(phenotypes{1}, '; ', treatments)];


save('Datasets/Phenotypes/2006_Gatbonton~Bedalov/gatbonton_bedalov_2006.mat','gatbonton_bedalov_2006');

% Save data into database
dt = gatbonton_bedalov_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Stepchenkova~Pavlov, 2005
% DATA = stepchenkova_pavlov_2005

stepchenkova_pavlov_2005.source = {'main PDF'};
stepchenkova_pavlov_2005.downloaddate = {'2014-03-10'};
stepchenkova_pavlov_2005.pmid = 15932646;

[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Stepchenkova~Pavlov/stepchenkova_pavlov_2005_hits.xlsx','Sheet1');

hits_genenames = data.raw(:,1);
inds = find(cellfun(@isnumeric, hits_genenames));
hits_genenames(inds) = [];
data.raw(inds,:) = [];

hits_orfs = genename2orf(hits_genenames,'noannot');

% Eliminate white spaces before/after ORF
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

hits_scores = cell2mat(data.raw(:,2:3));

% Normalize to WT
inds = find(strcmp('WT', hits_orfs));
hits_scores = hits_scores ./ repmat(hits_scores(inds,:), length(hits_orfs),1);
hits_scores(inds,:) = [];
hits_orfs(inds) = [];

phenotypes = {'growth';'mutagenicity'};
treatments = {'HAP','HAP'};


% Load tested genes
[t1,tested_orfs, t2,t3,t4,t5,t6,t7] = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Stepchenkova~Pavlov/mat_alpha_041902.txt', ...
    '%s %s %s %s %s %s %s %s','delimiter','\t');

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

tested_orfs = unique(upper(tested_orfs));

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
stepchenkova_pavlov_2005.orfs = tested_orfs;
stepchenkova_pavlov_2005.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(stepchenkova_pavlov_2005.orfs, hits_orfs);
stepchenkova_pavlov_2005.data(ind1,:) = hits_scores(ind2,:);

stepchenkova_pavlov_2005.ph = [strcat(phenotypes, '; ', treatments)];


save('Datasets/Phenotypes/2005_Stepchenkova~Pavlov/stepchenkova_pavlov_2005.mat','stepchenkova_pavlov_2005');

% Save data into database
dt = stepchenkova_pavlov_2005;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Outten~Culotta, 2005
% DATA = outten_culotta_2005

outten_culotta_2005.source = {'http://www.biochemj.org/bj/388/bj3880093add.pdf'};
outten_culotta_2005.downloaddate = {'2014-03-10'};
outten_culotta_2005.pmid = 15641941;

[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Outten~Culotta/outten_culotta_2005.xlsx','Sheet1');

hits_orfs = data.raw(:,1);

% Eliminate white spaces before/after ORF
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

hits_scores = cell2mat(data.raw(:,2));

% Adjust scores such that -1 = weakest phenotype, -4 = strongest phenotype.
hits_scores = hits_scores - 5;

% The 6 hits with no score (not sure why), set to 0
hits_scores(isnan(hits_scores)) = 0;

phenotypes = {'growth'};
treatments = {'hyperoxia'};


% Load tested genes
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Outten~Culotta/Yeast Knockout -BY4741.xlsx','mat_a_060701.txt');

tested_orfs = data.raw(2:end,2);

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

tested_orfs = unique(upper(tested_orfs));

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% From supplement: MED2 (YDL005C) not available in Mat-a background, was used Mat-alpha
% instead. So it should be added it to the list of tested
tested_orfs = [tested_orfs; {'YDL005C'}];

% YHR039C-A is not in the list of tested set, but it has the alias
% YHR039C-B, which is in the tested set. So, these are likely to be the same gene, so YHR039C-A should be renamed.
hits_orfs(strcmp('YHR039C-A', hits_orfs)) = {'YHR039C-B'};

% Create dataset
outten_culotta_2005.orfs = tested_orfs;
outten_culotta_2005.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(outten_culotta_2005.orfs, hits_orfs);
outten_culotta_2005.data(ind1,:) = hits_scores(ind2,:);

outten_culotta_2005.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2005_Outten~Culotta/outten_culotta_2005.mat','outten_culotta_2005');

% Save data into database
dt = outten_culotta_2005;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Reiner~Schneiter, 2006
% DATA = reiner_schneiter_2006

reiner_schneiter_2006.pmid = 16251356;

phenotypes = {'growth'};
treatments = {'hypoxia'};

% Load tested genes
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Reiner~Schneiter/BY4741-MATa COLLECTION.xls','chr11_1yes');

tested_orfs = data.raw(2:end,2);

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

tested_orfs = unique(upper(tested_orfs));

% Load hits
hits = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Reiner~Schneiter/reiner_scheiter_2006_hits.txt','%s');

% This list of ORFs is lacking the last character (published that way), so
% we have to match it to the list of tested strains.
for i = 1 : length(hits)
    inds = find(strncmp(hits{i}, tested_orfs, length(hits{i})));
    if length(inds) == 1
        hits_orfs(i) = tested_orfs(inds);
    else
        fprintf('%s\t%d\n', hits{i}, length(inds));
    end
end

% Two ORFs (YBR039W and YNL243W) could not be found in the list of tested
% strains. So we have to remove them.
hits_orfs(strncmp('YBR039', hits, length('YBR039')) | strncmp('YNL243', hits, length('YNL243'))) = [];

% Eliminate white spaces before/after ORF
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

% Assign a -1 to all the strains.
hits_scores = zeros(length(hits_orfs),1)-1;

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
reiner_schneiter_2006.orfs = tested_orfs;
reiner_schneiter_2006.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(reiner_schneiter_2006.orfs, hits_orfs);
reiner_schneiter_2006.data(ind1,:) = hits_scores(ind2,:);

reiner_schneiter_2006.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2006_Reiner~Schneiter/reiner_schneiter_2006.mat','reiner_schneiter_2006');

% Save data into database
dt = reiner_schneiter_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Galvan~Smith, 2008
% DATA = galvan_smith_2008

galvan_smith_2008.pmid = 17950387;

phenotypes = {'growth'};
treatments = {'chimaphilin'};

% Load plate maps
[map.txt, map.num, map.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Galvan~Smith/yGDA-Master_Plate_list_Combined(New).xlsx','Sheet1');
map.raw(1,:) = [];
map.platerowcol = cell2mat(map.raw(:,4:6));

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Galvan~Smith/IMELDA 08Feb2006GDAraw data.xls','Sheet1');
data.raw(1,:) = [];
data.platerowcol = cell2mat(data.raw(:,1:3));
data.orfs = cell(size(data.raw,1),1);

[C,ia,ib] = intersect(data.platerowcol,map.platerowcol,'rows');
data.orfs(ia) = map.raw(ib,2);

data.scores = cell2mat(data.raw(:,10:11));
data.scores_norm = data.scores(:,2)./data.scores(:,1);  % Test normalized divided by control normalized

inds = find(cellfun(@isnumeric, data.orfs));
data.orfs(inds) = [];
data.scores_norm(inds) = [];

data.orfs = cellfun(@strtrim, data.orfs,'UniformOutput',0);

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', data.orfs,1));
data.orfs(inds) = [];
data.scores_norm(inds) = [];

data.orfs = upper(data.orfs);


% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data.scores_norm, data.orfs, {'gname','mean'});
galvan_smith_2008.orfs = t;
galvan_smith_2008.data = t2;
galvan_smith_2008.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2008_Galvan~Smith/galvan_smith_2008.mat','galvan_smith_2008');

% Save data into database
dt = galvan_smith_2008;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Gonzalez-Ramos~Daran, 2013
% DATA = gonzalez_ramos_daran_2003

gonzalez_ramos_daran_2003.pmid = 23552365;

phenotypes = {'growth'};
treatments = {'butanol'};

% Load plate maps
[map.txt, map.num, map.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2013_Gonzalez-Ramos~Daran/Knockout collection map.xls','DATA');
map.raw(1,:) = [];
inds = find(cellfun(@isnumeric, map.raw(:,6)));
map.raw(inds,:) = [];
[t,ia,map.rowNum] = unique(map.raw(:,6));
map.platerowcol = [cell2mat(map.raw(:,5)) map.rowNum cell2mat(map.raw(:,7))];

% Load data
datafile = '/Users/Anastasia/Laboratory/Datasets/Phenotypes/2013_Gonzalez-Ramos~Daran/Butanol tolerance. Screening knockout collection.xlsx';
[status,sheets] = xlsfinfo(datafile);
all_data.plate = [];
all_data.rows = [];
all_data.cols = [];
all_data.norm_data = [];

for i = 1 : length(sheets)
    control_plate = zeros(8,12,2)+NaN;
    test_plate = zeros(8,12,2)+NaN;
    plateNum = str2num(sheets{i}(3:min(length(sheets{i}),4)));
    
    [data.txt, data.num, data.raw] = ...
        xlsread(datafile,sheets{i},'C7:N14');
    control_plate(:,:,1) = cell2mat(data.raw);
 
    [data.txt, data.num, data.raw] = ...
        xlsread(datafile,sheets{i},'Q7:AB14');
    control_plate(:,:,2) = cell2mat(data.raw);
    
    [data.txt, data.num, data.raw] = ...
        xlsread(datafile,sheets{i},'C19:N26');
    test_plate(:,:,2) = cell2mat(data.raw);
    
    [data.txt, data.num, data.raw] = ...
        xlsread(datafile,sheets{i},'Q19:AB26');
    test_plate(:,:,2) = cell2mat(data.raw);
    
    norm_data = nanmean(test_plate,3)./nanmean(control_plate,3);
    [cols,rows] = meshgrid(1:12,1:8);
    
    all_data.plate = [all_data.plate; zeros(length(norm_data(:)),1)+plateNum];
    all_data.rows = [all_data.rows; rows(:)];
    all_data.cols = [all_data.cols; cols(:)];
    all_data.norm_data = [all_data.norm_data; norm_data(:)];
    i
end

all_data.platerowcol = [all_data.plate all_data.rows all_data.cols];

% Average replicates
all_data.platerowcol2 = cell(length(all_data.plate),1);
for i = 1 : length(all_data.plate)
    all_data.platerowcol2{i} = [num2str(all_data.plate(i)),'-', num2str(all_data.rows(i)),'-', num2str(all_data.cols(i))];
end
[t,t2] = grpstats(all_data.norm_data, all_data.platerowcol2, {'gname','mean'});

all_data2.norm_data = t2;

tmp = regexp(t,'-','split');
tmp2 = vertcat(tmp{:});
all_data2.platerowcol = cellfun(@str2num, tmp2);

% Not all positions that are in the DATA have corresponding entries in the
% MAP. But that is somewhat expected, because not all positions on the
% plates are filled with strains.

[C,ia,ib] = intersect(all_data2.platerowcol,map.platerowcol,'rows');
all_data2.orfs = cell(size(all_data2.platerowcol,1),1);
all_data2.orfs(ia) = map.raw(ib,2);

% Eliminate empty positions
inds = find(cellfun(@isempty, all_data2.orfs));
all_data2.orfs(inds) = [];
all_data2.platerowcol(inds,:) = [];
all_data2.norm_data(inds) = [];

inds = find(cellfun(@isnumeric, all_data2.orfs));
all_data2.orfs(inds) = [];
all_data2.platerowcol(inds,:) = [];
all_data2.norm_data(inds) = [];

all_data2.orfs = cellfun(@strtrim, all_data2.orfs,'UniformOutput',0);

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', all_data2.orfs,1));
all_data2.orfs(inds) = [];
all_data2.platerowcol(inds,:) = [];
all_data2.norm_data(inds) = [];

all_data2.orfs = upper(all_data2.orfs);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(all_data2.norm_data, all_data2.orfs, {'gname','mean'});
gonzalez_ramos_daran_2003.orfs = t;
gonzalez_ramos_daran_2003.data = t2;
gonzalez_ramos_daran_2003.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2013_Gonzalez-Ramos~Daran/gonzalez_ramos_daran_2003.mat','gonzalez_ramos_daran_2003');

% Save data into database
dt = gonzalez_ramos_daran_2003;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Huang~Paulovich, 2013
% DATA = huang_paulovich_2013

huang_paulovich_2013.pmid = 23382077;

phenotypes = {'growth'};
treatments = {'MMS'};

% Load plate maps
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2013_Huang~Paulovich/Mat_a_obs_v4 0.xls','DATA');
tested_orfs = tested.raw(2:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

tested_orfs = unique(upper(tested_orfs));


% Load data
hits_genenames = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2013_Huang~Paulovich/huang_paulovich_2013_hits.txt','%s');
hits_orfs = genename2orf(hits_genenames,'noannot');
hits_scores = -ones(length(hits_orfs),1);

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
huang_paulovich_2013.orfs = tested_orfs;
huang_paulovich_2013.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(huang_paulovich_2013.orfs, hits_orfs);
huang_paulovich_2013.data(ind1,:) = hits_scores(ind2,:);

huang_paulovich_2013.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2013_Huang~Paulovich/huang_paulovich_2013.mat','huang_paulovich_2013');

% Save data into database
dt = huang_paulovich_2013;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Hirasawa~Shimizu, 2013
% DATA = hirasawa_shimizu_2013

hirasawa_shimizu_2013.pmid = 23665193;

phenotypes = {'lactate production'};
treatments = {''};

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2013_Hirasawa~Shimizu/data_from_S1_PDF.xlsx','Mutants');
data.raw(1,:) = [];
data2.genenames = [data.raw(:,1); data.raw(:,4); data.raw(:,7); data.raw(:,10)];
data2.data = [data.raw(:,2:3); data.raw(:,5:6); data.raw(:,8:9); data.raw(:,11:12)];

inds = find(cellfun(@isnumeric, data2.genenames));
data2.genenames(inds) = [];
data2.data(inds,:) = [];

data2.genenames = cellfun(@strtrim, data2.genenames,'UniformOutput',0);

data2.orfs = genename2orf(data2.genenames,'noannot');

% Adjustments
data2.orfs(strcmp('ARR4', data2.orfs)) = {'YDL100C'};
data2.orfs(strcmp('FMP12', data2.orfs)) = {'YHL021C'};
data2.orfs(strcmp('FMP13', data2.orfs)) = {'YKR016W'};
data2.orfs(strcmp('FMP14', data2.orfs)) = {'YPL099C'};
data2.orfs(strcmp('FMP22', data2.orfs)) = {'YHR198C'};
data2.orfs(strcmp('FMP24', data2.orfs)) = {'YMR115W'};
data2.orfs(strcmp('FMP26', data2.orfs)) = {'YJR080C'};
data2.orfs(strcmp('FMP29', data2.orfs)) = {'YER080W'};
data2.orfs(strcmp('FMP34', data2.orfs)) = {'YHR199C'};
data2.orfs(strcmp('FMP36', data2.orfs)) = {'YDR493W'};
data2.orfs(strcmp('FMP38', data2.orfs)) = {'YOR205C'};
data2.orfs(strcmp('FMP39', data2.orfs)) = {'YMR157C'};
data2.orfs(strcmp('FMP50', data2.orfs)) = {'YKR027W'};
data2.orfs(strcmp('FMP51', data2.orfs)) = {'YBR262C'};
data2.orfs(strcmp('FUN34', data2.orfs)) = {'YNR002C'};
data2.orfs(strcmp('GRD19', data2.orfs)) = {'YOR357C'};
data2.orfs(strcmp('MDM39', data2.orfs)) = {'YGL020C'};
data2.orfs(strcmp('MSU1', data2.orfs)) = {'YMR287C'};
data2.orfs(strcmp('PRP12', data2.orfs)) = {'YMR302C'};
data2.orfs(strcmp('RCS1', data2.orfs)) = {'YGL071W'};
data2.orfs(strcmp('RMD7', data2.orfs)) = {'YER083C'};
data2.orfs(strcmp('SOY1', data2.orfs)) = {'YBR194W'};
data2.orfs(strcmp('ZMS1', data2.orfs)) = {'YJR127C'};

data2.orfs = upper(data2.orfs);

inds = find(~strncmp('Y', data2.orfs,1));
data2.orfs(inds) = [];
data2.genenames(inds) = [];
data2.data(inds,:) = [];

data2.data = cell2mat(data2.data);
data2.data = nanmean(data2.data,2);

% Load controls
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2013_Hirasawa~Shimizu/data_from_S1_PDF.xlsx','CTRL');
data.raw(1,:) = [];
ctrl_data = [data.raw(:,2:3); data.raw(:,5:6); data.raw(:,8:9); data.raw(:,11:12)];
ctrl_data = reshape(ctrl_data,[],1);
inds = find(~cellfun(@isnumeric, ctrl_data));
ctrl_data(inds) = [];
ctrl_data = cell2mat(ctrl_data);
ctrl_data = nanmean(ctrl_data);

data2.data_norm = data2.data ./ ctrl_data;

% Average replicates
[t,t2] = grpstats(data2.data_norm, data2.orfs, {'gname','mean'});
hirasawa_shimizu_2013.orfs = t;
hirasawa_shimizu_2013.data = t2;
hirasawa_shimizu_2013.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2013_Hirasawa~Shimizu/hirasawa_shimizu_2013.mat','hirasawa_shimizu_2013');

% Save data into database
dt = hirasawa_shimizu_2013;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));



%% Schmidt~Boyer, 2012
% DATA = schmidt_boyer_2012

schmidt_boyer_2012.pmid = 22902726;

phenotypes = {'growth (OD)';'growth (colony size)';'growth (MIC)'};
treatments = {'boric acid'};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2012_Schmidt~Boyer/BA sensitivity comprehensive data file.xlsx','Strain list');
tested_orfs = tested.raw(2:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2012_Schmidt~Boyer/BA sensitivity comprehensive data file.xlsx','Sheet1');

% Dataset1: Tested = all; hits = liquid
hits_orfs = data.raw(3:end,1);
inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_orfs = unique(upper(hits_orfs));
hits_scores = -ones(length(hits_orfs),1);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

schmidt_boyer_2012.orfs = tested_orfs;
schmidt_boyer_2012.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
schmidt_boyer_2012.data(ind2,1) = hits_scores(ind1);

% Dataset2: Tested = dataset1, hits = solid
hits_orfs = data.raw(3:end,3);
inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_orfs = unique(upper(hits_orfs));
hits_scores = -ones(length(hits_orfs),1);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

inds = find(schmidt_boyer_2012.data(:,1)==0);
schmidt_boyer_2012.data(inds,2) = NaN;
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
schmidt_boyer_2012.data(ind2,2) = hits_scores(ind1);


% Dataset3: Tested = dataset2, data = MIC
hits_orfs = data.raw(3:end,5);
hits_scores = data.raw(3:end,7);

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = unique(upper(hits_orfs));
hits_scores = -cell2mat(hits_scores); %MIC values transformed to -MIC such that lower values correspond to the most sensitive strains (high original MIC)

[missing, ix] = setdiff(hits_orfs, tested_orfs);

inds = find(schmidt_boyer_2012.data(:,2)==0 | isnan(schmidt_boyer_2012.data(:,2)));
schmidt_boyer_2012.data(inds,3) = NaN;
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
schmidt_boyer_2012.data(ind2,3) = hits_scores(ind1);

schmidt_boyer_2012.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2012_Schmidt~Boyer/schmidt_boyer_2012.mat','schmidt_boyer_2012');


% Save data into database
dt = schmidt_boyer_2012;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
adj_ix = [1 3 2];
[~, adj_ix] = sort(adj_ix);

datasets.names(database_ix(adj_ix),:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix(adj_ix)));

%% Cheng~Bakalinsky, 2007
% DATA = cheng_bakalinsky_2007

cheng_bakalinsky_2007.pmid = 17644632;

phenotypes = {'growth (MIC)'};
treatments = {'oxalic acid'};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2007_Cheng~Bakalinsky/YSC1054Y.copy.xlsx','mat_alpha_obs');
tested_orfs = tested.raw(2:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2007_Cheng~Bakalinsky/Supplementary_Table_1.xlsx','Table 1');

% Dataset1: Tested = all; hits = liquid
hits_orfs = data.raw(4:end,2);
hits_scores = data.raw(4:end,3);
inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = upper(hits_orfs);
hits_scores = cellfun(@str2num, hits_scores);

% Many orfs are listed multiple times, but their scores are always the
% same. So, I can easily take the first of their values

[hits_orfs, ia,ic] = unique(hits_orfs);
hits_scores = hits_scores(ia);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

cheng_bakalinsky_2007.orfs = tested_orfs;
cheng_bakalinsky_2007.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
cheng_bakalinsky_2007.data(ind2,1) = hits_scores(ind1);

cheng_bakalinsky_2007.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2007_Cheng~Bakalinsky/cheng_bakalinsky_2007.mat','cheng_bakalinsky_2007');

% Save data into database
dt = cheng_bakalinsky_2007;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Gresham~Botstein, 2011
% DATA = gresham_botstein_2011

gresham_botstein_2011.pmid = 20944018;

phenotypes = {'Death rate (%/hr)'};
treatments = {'phosphate starvation';'leucine starvation'};

% Dataset #1
[data.num,data.txt,data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Gresham~Botstein/TABLES4.xlsx','phoAbs.txt');

data.genenames = data.raw(2:end, 1);
data.data = cell2mat(data.raw(2:end, 2));
data.genenames_noannot = cell(size(data.genenames));
% Eliminate the "_p" suffix from the genenames
for i = 1 : length(data.genenames)
    C = regexp(data.genenames{i},'_','split');
    data.genenames_noannot{i} = C{1};
end

data.orfs = genename2orf(data.genenames_noannot,'noannot');
% Eliminate the occasional suffix added by the renaming script
for i = 1 : length(data.orfs)
    C = regexp(data.orfs{i},'_', 'split');
    data.orfs{i} = C{1};
end
inds = find(~strncmp('Y', data.orfs,1));

unique(data.orfs(inds))


% Manual genename-to-ORF fixes
fid =fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Gresham~Botstein/genename_to_orf_fixes.txt');
C = textscan(fid, '%s\t%s\n');
fclose(fid);
for i = 1 : length(C{1})
    inds = strcmp(C{1}(i), data.orfs);
    data.orfs(inds) = C{2}(i);
end

inds = find(cellfun(@isnumeric, data.orfs));

data.orfs = cellfun(@strtrim, data.orfs,'UniformOutput',0);




% Dataset #2
[data2.num,data2.txt,data2.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Gresham~Botstein/TABLES5.xlsx','leuAbs.txt');

data2.genenames = data2.raw(2:end, 1);
data2.data = cell2mat(data2.raw(2:end, 2));
data2.genenames_noannot = cell(size(data2.genenames));
% Eliminate the "_p" suffix from the genenames
for i = 1 : length(data2.genenames)
    C = regexp(data2.genenames{i},'_','split');
    data2.genenames_noannot{i} = C{1};
end

data2.orfs = genename2orf(data2.genenames_noannot,'noannot');
% Eliminate the occasional suffix added by the renaming script
for i = 1 : length(data2.orfs)
    C = regexp(data2.orfs{i},'_', 'split');
    data2.orfs{i} = C{1};
end
inds = find(~strncmp('Y', data2.orfs,1));

unique(data2.orfs(inds))

% Manual genename-to-ORF fixes
fid =fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Gresham~Botstein/genename_to_orf_fixes.txt');
C = textscan(fid, '%s\t%s\n');
fclose(fid);
for i = 1 : length(C{1})
    inds = strcmp(C{1}(i), data2.orfs);
    data2.orfs(inds) = C{2}(i);
end

inds = find(cellfun(@isnumeric, data2.orfs));

data2.orfs = cellfun(@strtrim, data2.orfs,'UniformOutput',0);

% Average data for identical ORFs that appear multiple times
[data.orfs_u,data.data_avg] = grpstats(data.data, data.orfs, {'gname','mean'});
[data2.orfs_u,data2.data_avg] = grpstats(data2.data, data2.orfs, {'gname','mean'});

gresham_botstein_2011.orfs = unique([data.orfs_u; data2.orfs_u]);
gresham_botstein_2011.data = nan(length(gresham_botstein_2011.orfs),2);

[~,ind1,ind2] = intersect(gresham_botstein_2011.orfs, data.orfs_u);
gresham_botstein_2011.data(ind1,1) = data.data_avg(ind2);
[~,ind1,ind2] = intersect(gresham_botstein_2011.orfs, data2.orfs_u);
gresham_botstein_2011.data(ind1,2) = data2.data_avg(ind2);

gresham_botstein_2011.ph = strcat(phenotypes, {'; '}, treatments);

save('Datasets/Phenotypes/2011_Gresham~Botstein/gresham_botstein_2011.mat','gresham_botstein_2011');

% Save data into database
dt = gresham_botstein_2011;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Sambade~Kane, 2005
% DATA = sambade_kane_2005

sambade_kane_2005.pmid = 15937126;

phenotypes = {'growth (colony size)'};
treatments = {'pH [7.5] CaCl2 [60 mM]'};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Sambade~Kane/ResGen 384 well set 14 plates.xlsx','ResGen MATa -384');
tested_orfs = tested.raw(4:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
hits_genenames = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Sambade~Kane/hits_genenames.txt','%s');

hits_genenames = upper(hits_genenames);
hits_orfs = genename2orf(hits_genenames,'noannot');

hits_orfs(strcmpi('ada3', hits_orfs)) = {'YDR176W'};
hits_orfs(strcmpi('cak1', hits_orfs)) = {'YFL029C'};
hits_orfs(strcmpi('cwh36', hits_orfs)) = {'YCL007C'};
hits_orfs(strcmpi('rcs1', hits_orfs)) = {'YGL071W'};
hits_orfs(strcmpi('rmd7', hits_orfs)) = {'YER083C'};
hits_orfs(strcmpi('vma1', hits_orfs)) = {'YDL185W'};
hits_orfs(strcmpi('vma12', hits_orfs)) = {'YKL119C'};
hits_orfs(strcmpi('vma16', hits_orfs)) = {'YHR026W'};
hits_orfs(strcmpi('vma3', hits_orfs)) = {'YEL027W'};

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

hits_scores = -ones(length(hits_orfs),1);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = upper(hits_orfs);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
hits_orfs(strcmpi('YHR039C-A', hits_orfs)) = {'YHR039C-B'};

[missing, ix] = setdiff(hits_orfs, tested_orfs);
tested_orfs = [tested_orfs; missing];   % 1 orfs to be added

sambade_kane_2005.orfs = tested_orfs;
sambade_kane_2005.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
sambade_kane_2005.data(ind2,1) = hits_scores(ind1);

sambade_kane_2005.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2005_Sambade~Kane/sambade_kane_2005.mat','sambade_kane_2005');

% Save data into database
dt = sambade_kane_2005;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Bishop~Avery, 2007
% DATA = bishop_avery_2007

bishop_avery_2007.pmid = 17176259;

phenotypes = {'growth (colony size)'};
treatments = {'Ni(NO3)2 [2.5-4 mM]'};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2007_Bishop~Avery/Mata_DeletionArray+slow_growers.xlsx','96');
tested_orfs = tested.raw(3:end,2);
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2007_Bishop~Avery/Mata_DeletionArray+slow_growers.xlsx','slow growers');
tested_orfs = [tested_orfs; tested.raw(3:end,2)];
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
hits_genenames_resistant = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2007_Bishop~Avery/hits_genenames_resistant.txt','%s');

hits_genenames_resistant = upper(hits_genenames_resistant);
hits_orfs_resistant = genename2orf(hits_genenames_resistant,'noannot');
hits_orfs_resistant = cellfun(@strtrim, hits_orfs_resistant,'UniformOutput',0);
hits_orfs_resistant = unique(upper(hits_orfs_resistant));

hits_scores_resistant = ones(length(hits_orfs_resistant),1);

inds = find(~strncmp('Y', hits_orfs_resistant,1));
hits_orfs_resistant(inds) = [];
hits_scores_resistant(inds) = [];


hits_genenames_sensitive = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2007_Bishop~Avery/hits_genenames_sensitive.txt','%s');

hits_genenames_sensitive = upper(hits_genenames_sensitive);
hits_orfs_sensitive = genename2orf(hits_genenames_sensitive,'noannot');

% Adjustments
hits_orfs_sensitive(strcmpi('FMP18', hits_orfs_sensitive)) = {'YKR065C'};

hits_orfs_sensitive = cellfun(@strtrim, hits_orfs_sensitive,'UniformOutput',0);
hits_orfs_sensitive = unique(upper(hits_orfs_sensitive));

hits_scores_sensitive = -ones(length(hits_orfs_sensitive),1);

inds = find(~strncmp('Y', hits_orfs_sensitive,1));
hits_orfs_sensitive(inds) = [];
hits_scores_sensitive(inds) = [];

overlapping = intersect(hits_orfs_resistant, hits_orfs_sensitive);

hits_orfs = [hits_orfs_resistant; hits_orfs_sensitive];
hits_scores = [hits_scores_resistant; hits_scores_sensitive];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
tested_orfs = [tested_orfs; missing];   % 5 orfs to be added

bishop_avery_2007.orfs = tested_orfs;
bishop_avery_2007.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
bishop_avery_2007.data(ind2,1) = hits_scores(ind1);

bishop_avery_2007.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2007_Bishop~Avery/bishop_avery_2007.mat','bishop_avery_2007');

% Save data into database
dt = bishop_avery_2007;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Serero~Boiteux, 2008
% DATA = serero_boiteux_2008

serero_boiteux_2008.pmid = 18514590;

phenotypes = {'growth (colony size)'};
treatments = {'CdCl2 [100 uM]'};

% Load tested

% ATTEMPT #1: go through the EXCEL files. PROBLEM: too many genes obtained
% (5617).
% home_dir = '/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Serero~Boiteux/Deletion_mutants_list/';
% folders = dir(home_dir);
% folders_names = {folders.name};
% folders_names(strncmp('.', folders_names,1)) = [];
% 
% tested_orfs = [];
% for i = 1 : length(folders_names)
%     % Find all Excel files in the folder
%     xls_files = dir([home_dir folders_names{i} '/*.XLS']);
%     xls_files_names = {xls_files.name};
%     
%     for j = 1 : length(xls_files_names)
%         filename = [home_dir folders_names{i} '/' xls_files_names{j}];
%         [status,sheets] = xlsfinfo(filename);
%         
%         for k = 1 : length(sheets)
%             [tested.txt, tested.num, tested.raw] = ...
%                 xlsread(filename,sheets{k});
%             if ~isempty(tested.raw)
%                 inds = find(strcmp('BY4741', tested.raw(:,3)));
%                 tested_orfs = [tested_orfs; tested.raw(inds,2)];
%             end
%         end
%     end
%     i
% end

% ATTEMPT #2: GO through the DOC files.
% cd ~/Laboratory/Datasets/Phenotypes/2008_Serero~Boiteux/Deletion_mutants_list/
% Covert all DOC files into TXT files by running: sudo textutil -convert txt */*.DOC
% Read the TXT files that end with "1"

home_dir = '/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Serero~Boiteux/Deletion_mutants_list/';
folders = dir(home_dir);
folders_names = {folders.name};
folders_names(strncmp('.', folders_names,1)) = [];

tested_orfs = [];
for i = 1 : length(folders_names)
    % Find all TXT files in the folder
    txt_files = dir([home_dir folders_names{i} '/*.txt']);
    txt_files_names = {txt_files.name};
    rtf_files = dir([home_dir folders_names{i} '/*.RTF']);
    txt_files_names = [txt_files_names; {rtf_files.name}];
    
    for j = 1 : length(txt_files_names)
        % If filenames ends in "~1" or "a", load the list
        t = regexp(txt_files_names{j},'\.','split');
        if strcmp(t{1}(end),'1') | strcmp(t{1}(end),'a')
            tst = textread([home_dir folders_names{i} '/' txt_files_names{j}],'%s');
            inds = find(strncmp('Y', tst,1));
            tested_orfs = [tested_orfs; tst(inds)];
        end
    end
    i
end
        

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
tested_orfs = upper(tested_orfs);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(tested_orfs);

% Load data
[hits_genenames, hits_scores_txt] = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Serero~Boiteux/hits_genenames.txt','%s %s','delimiter','\t');

hits_scores = cellfun(@length, hits_scores_txt)+1;
hits_scores = -hits_scores;

hits_genenames = upper(hits_genenames);
hits_orfs = genename2orf(hits_genenames,'noannot');

% Adjustments
hits_orfs(strcmpi('lys7', hits_orfs)) = {'YMR038C'};
hits_orfs(strcmpi('trf4', hits_orfs)) = {'YOL115W'};

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);


serero_boiteux_2008.orfs = tested_orfs;
serero_boiteux_2008.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
serero_boiteux_2008.data(ind2,1) = hits_scores(ind1);

serero_boiteux_2008.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2008_Serero~Boiteux/serero_boiteux_2008.mat','serero_boiteux_2008');

% Save data into database
dt = serero_boiteux_2008;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Matsufuji~Nakagawa, 2008
% DATA = matsufuji_nakagawa_2008

matsufuji_nakagawa_2008.pmid = 19061187;

phenotypes = {'growth (colony size)'};
treatments = {'acetaldehyde [35 mM]'};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Matsufuji~Nakagawa/S.c 5000.xlsx','remake');
tested_orfs = tested.raw(3:end,3);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
hits_orfs = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Matsufuji~Nakagawa/hits_orfs.txt','%s');

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
hits_orfs = unique(upper(hits_orfs));

% Correction: ORF YHR041C was in the hit list twice. One instance should
% have been YDR378C.

hits_orfs = [hits_orfs; {'YDR378C'}];

hits_scores = -ones(length(hits_orfs),1);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
tested_orfs = [tested_orfs; missing];   % 1 orfs to be added

matsufuji_nakagawa_2008.orfs = tested_orfs;
matsufuji_nakagawa_2008.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
matsufuji_nakagawa_2008.data(ind2,1) = hits_scores(ind1);

matsufuji_nakagawa_2008.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2008_Matsufuji~Nakagawa/matsufuji_nakagawa_2008.mat','matsufuji_nakagawa_2008');

% Save data into database
dt = matsufuji_nakagawa_2008;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Kemmer~Roberge, 2009
% DATA = kemmer_roberge_2009

kemmer_roberge_2009.pmid = 19144191;

phenotypes = {'growth (colony size)'};
treatments = {'dhMotC [60 uM]'};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Kemmer~Roberge/haploid set.xlsx','haploid set');
tested_orfs = tested.raw(6:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
hits_genenames = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Kemmer~Roberge/hits_genenames.txt','%s');

hits_genenames = cellfun(@strtrim, hits_genenames,'UniformOutput',0);
hits_orfs = genename2orf(hits_genenames,'noannot');
hits_orfs = unique(upper(hits_orfs));

hits_scores = -ones(length(hits_orfs),1);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
hits_orfs(strcmp('YHR039C-A', hits_orfs)) = {'YHR039C-B'};

kemmer_roberge_2009.orfs = tested_orfs;
kemmer_roberge_2009.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
kemmer_roberge_2009.data(ind2,1) = hits_scores(ind1);

kemmer_roberge_2009.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2009_Kemmer~Roberge/kemmer_roberge_2009.mat','kemmer_roberge_2009');

% Save data into database
dt = kemmer_roberge_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Xia~Flores-Rozas, 2007
% DATA = xia_flores_rozas_2007

xia_flores_rozas_2007.pmid = 18056469;

phenotypes = {'growth (colony size)'};
treatments = {'doxorubicin [20 umol/L]'};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2007_Xia~Flores-Rozas/Mat_a.xlsx','mat_a_041902');
tested_orfs = tested.raw(3:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[hits_genenames, hits_scores_txt] = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2007_Xia~Flores-Rozas/hits_genenames.txt','%s %s','delimiter','\t');

hits_genenames = cellfun(@strtrim, hits_genenames,'UniformOutput',0);
hits_orfs = genename2orf(hits_genenames,'noannot');
hits_scores = -cellfun(@length, hits_scores_txt);

length(unique(upper(hits_orfs)))

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);


xia_flores_rozas_2007.orfs = tested_orfs;
xia_flores_rozas_2007.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
xia_flores_rozas_2007.data(ind2,1) = hits_scores(ind1);

xia_flores_rozas_2007.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2007_Xia~Flores-Rozas/xia_flores_rozas_2007.mat','xia_flores_rozas_2007');

% Save data into database
dt = xia_flores_rozas_2007;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Dos Santos~Sa-Correia, 2011
% DATA = dos_santos_sa_correia_2011

dos_santos_sa_correia_2011.pmid = 21960436;

phenotypes = {'growth [spot assay]'};
treatments = {'quinine [1.5-1.7 g/L]'};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Dos Santos~Sa-Correia/List of strains tested.xlsx','Tabelle2');
tested_orfs = tested.raw(2:end,1);
% slow_growers = tested.raw(2:end,2);
% inds = find(~cellfun(@isnumeric, slow_growers));
% tested_orfs(inds) = [];

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
hits_genenames_HS = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Dos Santos~Sa-Correia/hits_genenames_hs.txt','%s');

hits_genenames_HS = cellfun(@strtrim, hits_genenames_HS,'UniformOutput',0);
hits_orfs_HS = genename2orf(hits_genenames_HS,'noannot');
hits_scores_HS = zeros(length(hits_orfs_HS),1)-2;

length(unique(upper(hits_orfs_HS)))

hits_genenames_S = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Dos Santos~Sa-Correia/hits_genenames_s.txt','%s');

hits_genenames_S = cellfun(@strtrim, hits_genenames_S,'UniformOutput',0);
hits_orfs_S = genename2orf(hits_genenames_S,'noannot');
hits_scores_S = zeros(length(hits_orfs_S),1)-1;

length(unique(upper(hits_orfs_S)))

hits_genenames_R = ...
    textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Dos Santos~Sa-Correia/hits_genenames_r.txt','%s');

hits_genenames_R = cellfun(@strtrim, hits_genenames_R,'UniformOutput',0);
hits_orfs_R = genename2orf(hits_genenames_R,'noannot');
hits_scores_R = zeros(length(hits_orfs_R),1)+1;

length(unique(upper(hits_orfs_R)))

% Adjust the overlapping strains as follows: eliminate HS from S; eliminate
% HS ^ R and S ^ R from both
[~,~,ind2] = intersect(hits_orfs_HS, hits_orfs_S);
hits_orfs_S(ind2) = []; hits_scores_S(ind2) = [];

[~,ind1,ind2] = intersect(hits_orfs_HS, hits_orfs_R);
hits_orfs_HS(ind1) = []; hits_scores_HS(ind1) = [];
hits_orfs_R(ind2) = []; hits_scores_R(ind2) = [];

[~,ind1,ind2] = intersect(hits_orfs_S, hits_orfs_R);
hits_orfs_S(ind1) = []; hits_scores_S(ind1) = [];
hits_orfs_R(ind2) = []; hits_scores_R(ind2) = [];

hits_orfs = [hits_orfs_HS; hits_orfs_S; hits_orfs_R];
hits_scores = [hits_scores_HS; hits_scores_S; hits_scores_R];

length(unique(hits_orfs))

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
hits_orfs(strcmp('YHR039C-A', hits_orfs)) = {'YHR039C-B'};
hits_orfs(strcmp('YBR089C-A', hits_orfs)) = {'YBR090C-A'};

dos_santos_sa_correia_2011.orfs = tested_orfs;
dos_santos_sa_correia_2011.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
dos_santos_sa_correia_2011.data(ind2,1) = hits_scores(ind1);

dos_santos_sa_correia_2011.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2011_Dos Santos~Sa-Correia/dos_santos_sa_correia_2011.mat','dos_santos_sa_correia_2011');

% Save data into database
dt = dos_santos_sa_correia_2011;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Ding~Bakalinsky, 2013
% DATA = ding_bakalinsky_2013

ding_bakalinsky_2013.pmid = 23828602;

phenotypes = {'growth [pooled CFU]'};
treatments = {'acetic acid [122.5 mM]'};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2013_Ding~Bakalinsky/YSC1054Y.copy.xlsx','mat_alpha_obs');
tested_orfs = tested.raw(2:end,2);

inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2013_Ding~Bakalinsky/hits_genenames_R.txt');
hits_genenames_R = textscan(fid,'%s');
fclose(fid);

hits_genenames_R = hits_genenames_R{1};

hits_genenames_R = cellfun(@strtrim, hits_genenames_R,'UniformOutput',0);
hits_genenames_R = unique(upper(hits_genenames_R));
hits_orfs_R = unique(genename2orf(hits_genenames_R,'noannot'));
hits_scores_R = zeros(length(hits_orfs_R),1)+1;

inds = find(~strncmp('Y', hits_orfs_R,1));
hits_orfs_R(inds) = [];
hits_scores_R(inds) = [];

[missing, ix] = setdiff(hits_orfs_R, tested_orfs);

ding_bakalinsky_2013.orfs = tested_orfs;
ding_bakalinsky_2013.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs_R, tested_orfs);
ding_bakalinsky_2013.data(ind2,1) = hits_scores_R(ind1);

ding_bakalinsky_2013.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2013_Ding~Bakalinsky/ding_bakalinsky_2013.mat','ding_bakalinsky_2013');

% Save data into database
dt = ding_bakalinsky_2013;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Postma~Ralser, 2009
% DATA = postma_ralser_2009

postma_ralser_2009.pmid = 20157578;

phenotypes = {'growth [spot assay]'};
treatments = {'hibernation'};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Postma~Ralser/Mat_a_obs_v2.0.xlsx','Mat_a_obs_v2.0.txt');
tested_orfs = tested.raw(2:end,1);


inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Postma~Ralser/hits_orfs.txt');
hits_orfs = textscan(fid,'%s');
fclose(fid);

hits_orfs = unique(hits_orfs{1});
hits_scores = zeros(length(hits_orfs),1)+1;

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustment
tested_orfs = [tested_orfs; missing];   % 1 ORF added

postma_ralser_2009.orfs = tested_orfs;
postma_ralser_2009.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
postma_ralser_2009.data(ind2,1) = hits_scores(ind1);

postma_ralser_2009.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2009_Postma~Ralser/postma_ralser_2009.mat','postma_ralser_2009');

% Save data into database
dt = postma_ralser_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Thevissen~Francois, 2007
% DATA = thevissen_francois_2007

thevissen_francois_2007.pmid = 17553796;

phenotypes = {'growth [MIC]'};
treatments = {'miconazole [0.025-12.5 ug/ml]'};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2007_Thevissen~Francois/Euroscarf library.xlsx','Tabelle1');
tested_orfs = tested.raw(2:end,2);


inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2007_Thevissen~Francois/hits_orfs_scores.txt');
hits = textscan(fid,'%s %d');
hits_orfs = upper(hits{1});
hits_scores = -hits{2};
fclose(fid);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustment
tested_orfs = [tested_orfs; missing];   % 2 ORFs added

thevissen_francois_2007.orfs = tested_orfs;
thevissen_francois_2007.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
thevissen_francois_2007.data(ind2,1) = hits_scores(ind1);

thevissen_francois_2007.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2007_Thevissen~Francois/thevissen_francois_2007.mat','thevissen_francois_2007');

% Save data into database
dt = thevissen_francois_2007;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Mir~Cashikar, 2009
% DATA = mir_cashikar_2009

mir_cashikar_2009.pmid = 18936161;

phenotypes = {'growth [MIC]'};
treatments = {'heat stress (temperature [50ºC], duration [30 min])'};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Mir~Cashikar/Mat_a.xlsx','mat_a_041902');
tested_orfs = tested.raw(3:end,2);


inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Mir~Cashikar/hits_orfs_scores.txt');
hits = textscan(fid,'%s %d');
hits_orfs = upper(hits{1});
hits_scores = -hits{2};
fclose(fid);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

mir_cashikar_2009.orfs = tested_orfs;
mir_cashikar_2009.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
mir_cashikar_2009.data(ind2,1) = hits_scores(ind1);

mir_cashikar_2009.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2009_Mir~Cashikar/mir_cashikar_2009.mat','mir_cashikar_2009');

% Save data into database
dt = mir_cashikar_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% van Voorst~Bradt, 2006
% DATA = van_voorst_bradt_2006
% TESTED = not available

van_voorst_bradt_2006.pmid = 16598687;

phenotypes = {'growth [streaks on agar]'};
treatments = {'ethanol [6%]'};

% Load data
fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_van_Voorst~Brandt/hits_genenames.txt');
hits = textscan(fid,'%s');
hits_genenames = hits{1};
fclose(fid);

hits_orfs = genename2orf(hits_genenames,'noannot');

% Adjustments
hits_orfs(strcmpi('ada3', hits_orfs)) = {'YDR176W'};
hits_orfs(strcmpi('vps39', hits_orfs)) = {'YDL077C'};

hits_scores = -ones(length(hits_orfs),1);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

van_voorst_bradt_2006.orfs = hits_orfs;
van_voorst_bradt_2006.data = hits_scores;

van_voorst_bradt_2006.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2006_van_Voorst~Brandt/van_voorst_bradt_2006.mat','van_voorst_bradt_2006');

% Save data into database
dt = van_voorst_bradt_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Kawahata~Iefuji, 2006
% DATA = kawahata_iefuji_2006
kawahata_iefuji_2006.pmid = 16911514;

phenotypes = {'growth [spot assay]'};
treatments = {'lactic acid [5.1% w/v], pH [2.7]'; ...
    'lactic acid [3.1% w/v], pH [2.9]'; ...
    'HCl [0.28% w/v], pH [2.4]';...
    'HCl [0.24% w/v], pH [2.6]'; ...
    'acetic acid [0.5% w/v], pH [4.2]'; ...
    'acetic acid [0.4% w/v], pH [4.3]'};

% Load resistant
[hits_resistant.txt, hits_resistant.num, hits_resistant.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Kawahata~Iefuji/hits.xlsx','Resistant');

hits_resistant_orfs = unique([hits_resistant.raw(:,1); hits_resistant.raw(:,6)]);

inds = cellfun(@isnumeric, hits_resistant_orfs);
hits_resistant_orfs(inds) = [];
hits_resistant_orfs = cellfun(@strtrim, hits_resistant_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_resistant_orfs,1);
hits_resistant_orfs(inds) = [];
hits_resistant_orfs = unique(upper(hits_resistant_orfs));

hits_resistant_scores = zeros(length(hits_resistant_orfs),3);

% Scores
tmp_orfs = [hits_resistant.raw(strcmp('+', hits_resistant.raw(:,3)),1); hits_resistant.raw(strcmp('+', hits_resistant.raw(:,8)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_resistant_orfs);
hits_resistant_scores(ind2,1) = 1;

tmp_orfs = [hits_resistant.raw(strcmp('+', hits_resistant.raw(:,4)),1); hits_resistant.raw(strcmp('+', hits_resistant.raw(:,9)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_resistant_orfs);
hits_resistant_scores(ind2,2) = 1;

tmp_orfs = [hits_resistant.raw(strcmp('+', hits_resistant.raw(:,5)),1); hits_resistant.raw(strcmp('+', hits_resistant.raw(:,10)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_resistant_orfs);
hits_resistant_scores(ind2,3) = 1;

% Load sensitive
[hits_sensitive.txt, hits_sensitive.num, hits_sensitive.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Kawahata~Iefuji/hits.xlsx','Sensitive');

hits_sensitive_orfs = [hits_sensitive.raw(:,1); hits_sensitive.raw(:,6)];

inds = cellfun(@isnumeric, hits_sensitive_orfs);
hits_sensitive_orfs(inds) = [];
hits_sensitive_orfs = cellfun(@strtrim, hits_sensitive_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_sensitive_orfs,1);
hits_sensitive_orfs(inds) = [];
hits_sensitive_orfs = unique(upper(hits_sensitive_orfs));

hits_sensitive_scores = zeros(length(hits_sensitive_orfs),3);

% Scores
tmp_orfs = [hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,3)),1); hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,8)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_sensitive_orfs);
hits_sensitive_scores(ind2,1) = -1;

tmp_orfs = [hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,4)),1); hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,9)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_sensitive_orfs);
hits_sensitive_scores(ind2,2) = -1;

tmp_orfs = [hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,5)),1); hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,10)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_sensitive_orfs);
hits_sensitive_scores(ind2,3) = -1;


% Check overlap between resistant and sensitive
length(intersect(hits_resistant_orfs(hits_resistant_scores(:,3)>0), hits_sensitive_orfs(hits_sensitive_scores(:,3)<0)))

kawahata_iefuji_2006.orfs = unique([hits_resistant_orfs; hits_sensitive_orfs]);
kawahata_iefuji_2006.data = zeros(length(kawahata_iefuji_2006.orfs), length(phenotypes));
[~,ind1,ind2] = intersect(hits_resistant_orfs, kawahata_iefuji_2006.orfs);
kawahata_iefuji_2006.data(ind2,[1 3 5]) = hits_resistant_scores(ind1,:);
[~,ind1,ind2] = intersect(hits_sensitive_orfs, kawahata_iefuji_2006.orfs);
kawahata_iefuji_2006.data(ind2,[2 4 6]) = hits_sensitive_scores(ind1,:);


kawahata_iefuji_2006.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2006_Kawahata~Iefuji/kawahata_iefuji_2006.mat','kawahata_iefuji_2006');

% Save data into database
dt = kawahata_iefuji_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Arita~Costa, 2009
% DATA = arita_costa_2009
arita_costa_2009.pmid = 19917080;

phenotypes = {'growth [spot assay]'};
treatments = {'NiSO4 [0.75-1.25 mM]'};

% Load resistant
[hits.txt, hits.num, hits.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Arita~Costa/1471-2164-10-524-S1.XLS','Sheet1');
hits_sensitive_orfs = hits.raw(7:end,1);
hits_resistant_orfs = hits.raw(7:end,2);

inds = cellfun(@isnumeric, hits_resistant_orfs);
hits_resistant_orfs(inds) = [];
hits_resistant_orfs = cellfun(@strtrim, hits_resistant_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_resistant_orfs,1);
hits_resistant_orfs(inds) = [];
hits_resistant_orfs = unique(upper(hits_resistant_orfs));

hits_resistant_scores = ones(length(hits_resistant_orfs),1);

inds = cellfun(@isnumeric, hits_sensitive_orfs);
hits_sensitive_orfs(inds) = [];
hits_sensitive_orfs = cellfun(@strtrim, hits_sensitive_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_sensitive_orfs,1);
hits_sensitive_orfs(inds) = [];
hits_sensitive_orfs = unique(upper(hits_sensitive_orfs));

hits_sensitive_scores = -ones(length(hits_sensitive_orfs),1);



% Check overlap between resistant and sensitive
length(intersect(hits_resistant_orfs, hits_sensitive_orfs))

arita_costa_2009.orfs = unique([hits_resistant_orfs; hits_sensitive_orfs]);
arita_costa_2009.data = zeros(length(arita_costa_2009.orfs), length(phenotypes));
[~,ind1,ind2] = intersect(hits_resistant_orfs, arita_costa_2009.orfs);
arita_costa_2009.data(ind2,1) = hits_resistant_scores(ind1);
[~,ind1,ind2] = intersect(hits_sensitive_orfs, arita_costa_2009.orfs);
arita_costa_2009.data(ind2,1) = hits_sensitive_scores(ind1);


arita_costa_2009.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2009_Arita~Costa/arita_costa_2009.mat','arita_costa_2009');

% Save data into database
dt = arita_costa_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Smith~Bakalinsky, 2013
% DATA = smith_bakalinsky_2013
smith_bakalinsky_2013.pmid = 23144132;

phenotypes = {'growth [CFU]'};
treatments = {'AuNP [10-100 ug/ml]'};

% Load tested (Same as Ding~Bakalinksy, 2013)
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2013_Ding~Bakalinsky/YSC1054Y.copy.xlsx','mat_alpha_obs');
tested_orfs = tested.raw(2:end,2);

inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load resistant
fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2013_Smith~Bakalinsky/hits_orfs.txt');
hits = textscan(fid,'%s');
hits_orfs = hits{1};
fclose(fid);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_orfs = unique(upper(hits_orfs));

hits_scores = ones(length(hits_orfs),1);


smith_bakalinsky_2013.orfs = tested_orfs;
smith_bakalinsky_2013.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
smith_bakalinsky_2013.data(ind2,1) = hits_scores(ind1);

smith_bakalinsky_2013.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2013_Smith~Bakalinsky/smith_bakalinsky_2013.mat','smith_bakalinsky_2013');

% Save data into database
dt = smith_bakalinsky_2013;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Zhou~Costa, 2009
% DATA = zhou_costa_2009
% TESTED = not available
zhou_costa_2009.pmid = 19631266;

phenotypes = {'growth [spot assay]'};
treatments = {'NaAsO2 [0.075-1 mM]'};

% Load hits
[hits.txt, hits.num, hits.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Zhou~Costa/table.xlsx','table.csv');
hits_orfs = hits.raw(4:end,1);
hits_scores = hits.raw(4:end,4);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = upper(hits_orfs);
hits_scores(strcmp('Above 5', hits_scores)) = {5};
hits_scores = cell2mat(hits_scores);

[hits_orfs, hits_scores] = grpstats(hits_scores,hits_orfs,{'gname','mean'});

wt = 4.47;  %IC50 from paper
hits_scores = hits_scores - wt; % negative scores = more sensitive than WT, and viceversa.


zhou_costa_2009.orfs = hits_orfs;
zhou_costa_2009.data = hits_scores;

zhou_costa_2009.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2009_Zhou~Costa/zhou_costa_2009.mat','zhou_costa_2009');

% Save data into database
dt = zhou_costa_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Orij~Smits, 2012
% DATA = orij_smits_2012

orij_smits_2012.pmid = 23021432;

phenotypes = {'cytosolic pH'};
treatments = {'standard'};

% Load hits
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2012_Orij~Smits/pH Screen raw.xlsx','initial screens');
hits_orfs = data.raw(2:end,1);
hits_scores = data.raw(2:end,3:4);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds,:) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_scores(inds,:) = [];

hits_orfs = upper(hits_orfs);
hits_scores = cell2mat(hits_scores);
hits_scores = nanmean(hits_scores, 2);

[hits_orfs, hits_scores] = grpstats(hits_scores,hits_orfs,{'gname','mean'});

orij_smits_2012.orfs = hits_orfs;
orij_smits_2012.data = hits_scores;

orij_smits_2012.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2012_Orij~Smits/orij_smits_2012.mat','orij_smits_2012');

% Save data into database
dt = orij_smits_2012;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Garay~DeLuna, 2014
% DATA = garay_deluna_2014

garay_deluna_2014.pmid = 24586198;

phenotypes = {'chronological life span'};
treatments = {'standard'};

% Load hits
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2014_Garay~DeLuna/journal.pgen.1004168.s013.xlsx','Genome-wide CLS screen');
hits_orfs = data.raw(2:end,1);
hits_scores = data.raw(2:end,3);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds,:) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_scores(inds,:) = [];

hits_orfs = upper(hits_orfs);
hits_scores = cell2mat(hits_scores);

[hits_orfs, hits_scores] = grpstats(hits_scores,hits_orfs,{'gname','mean'});

garay_deluna_2014.orfs = hits_orfs;
garay_deluna_2014.data = hits_scores;

garay_deluna_2014.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2014_Garay~DeLuna/garay_deluna_2014.mat','garay_deluna_2014');

% Save data into database
dt = garay_deluna_2014;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

% %% Wolinski~Kohlwein, 2009
% % DATA = wolinski_kohlwein_2009
% 
% wolinski_kohlwein_2009.pmid = 19118449;
% 
% phenotypes = {'PST1 localization to peroxisome'};
% treatments = {'standard'};
% 
% % Load hits
% [data.txt, data.num, data.raw] = ...
%     xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Wolinski~Kohlwein/523strainquant_v2.xlsx','Tabelle1');
% hits_orfs = data.raw(2:end,1);
% hits_scores = data.raw(2:end,2);
% 
% tmp = regexp(hits_orfs,'-','split');
% for i = 1 : length(tmp)
%     hits_orfs(i) = tmp{i}(1);
% end
% 
% inds = cellfun(@isnumeric, hits_orfs);
% hits_orfs(inds) = [];
% hits_scores(inds,:) = [];
% hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
% inds = ~strncmp('Y', hits_orfs,1);
% hits_orfs(inds) = [];
% hits_scores(inds,:) = [];
% 
% hits_orfs = upper(hits_orfs);
% hits_scores = cell2mat(hits_scores);
% 
% [hits_orfs, hits_scores] = grpstats(hits_scores,hits_orfs,{'gname','mean'});
% 
% wolinski_kohlwein_2009.orfs = hits_orfs;
% wolinski_kohlwein_2009.data = hits_scores;
% wolinski_kohlwein_2009.ph = [strcat(phenotypes, '; ', treatments)];
% 
% save('Datasets/Phenotypes/2009_Wolinski~Kohlwein/wolinski_kohlwein_2009.mat','wolinski_kohlwein_2009');

%% Huang~O'Shea, 2005
% DATA = huang_oshea_2005

huang_oshea_2005.pmid = 15695358;

phenotypes = {'Pho5 phosphatase activity'};
treatments = {'phosphate starvation'};

plates = [1:51,70:71];
rows = 'ABCDEFGH';
cols = 1:12;
timepoints = [0:120:360];

[map.txt, map.num, map.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Huang~OShea/MATa Collection.xlsx','MATa Collection.xls');

map.orf = map.raw(5:end,2);
map.plate = cell2mat(map.raw(5:end, 5));
map.row = map.raw(5:end,6);
map.col = cell2mat(map.raw(5:end,7));

inds = find(cellfun(@isnan, map.row));
map.orf(inds) = [];
map.plate(inds) = [];
map.row(inds) = [];
map.col(inds) = [];

for i = 1 : length(map.row)
    map.row{i} = findstr(map.row{i}, rows);
end

map.row = cell2mat(map.row);
map.idx = sub2ind([8 12],map.row, map.col);
map.data = zeros(length(map.orf),4);

inds = cellfun(@isnumeric, map.orf);

map.orf = cellfun(@strtrim, map.orf,'UniformOutput',0);
inds = ~strncmp('Y', map.orf,1);


[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Huang~OShea/Initial Screen Data.xls','Sheet1');

[all_cols, all_rows] = meshgrid(1:12,1:8);
all_idx = sub2ind([8 12], all_rows, all_cols)';
[~,ix] = sort(all_idx(:));

for i = 1 : length(plates)
    inds = find(map.plate == plates(i));
    for t = 1 : length(timepoints)
        this_data = cell2mat(data.raw(6+(i-1)*7+(t-1),2:end));
        map.data(inds,t) = this_data(ix(map.idx(inds)));
    end
end

map.data_avg = nanmean(map.data,2);

% Average multiple occurences of the same ORF
[map.orf_u, map.data_avg_avg] = grpstats(map.data_avg,map.orf,{'gname','mean'});

huang_oshea_2005.orfs = map.orf_u;
huang_oshea_2005.data = map.data_avg_avg;
huang_oshea_2005.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2005_Huang~Oshea/huang_oshea_2005.mat','huang_oshea_2005');

% Save data into database
dt = huang_oshea_2005;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Freimoser~Amrhein, 2006
% DATA = freimoser_amrhein_2006

freimoser_amrhein_2006.pmid = 17107617;

phenotypes = {'inorganic polyphosphate abundance'};
treatments = {''};

folder = '/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Freimoser~Amrhein/';

[map.txt, map.num, map.raw] = xlsread([folder 'Strains_plates_annotation.xlsx']);

map.orf = map.raw(4:end, 4);
map.plate = map.raw(4:end,1);
map.row = map.raw(4:end,2);
map.col = map.raw(4:end,3);

inds = cellfun(@isnumeric, map.orf);
map.orf(inds) = [];
map.plate(inds) = [];
map.row(inds) = [];
map.col(inds) = [];

inds = ~strncmp('Y', map.orf,1);
map.orf(inds) = [];
map.plate(inds) = [];
map.row(inds) = [];
map.col(inds) = [];

map.plate = cell2mat(map.plate);
map.col = cell2mat(map.col);

files = dir([folder 'Archive/']);
file_names = {files.name};
file_names = file_names(strncmp('Plate', file_names,5));

fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Freimoser~Amrhein/columns_to_extract.txt');
C = textscan(fid,'%s %s %s %s %s %s','delimiter','\t');
fclose(fid);

all_data.plate = {};
all_data.row = {};
all_data.col = {};
all_data.data = {};
for i = 1 : length(file_names)
    [~, sheets] = xlsfinfo([folder 'Archive/' file_names{i}]);
    ind_file = find(strcmp(file_names{i}, C{1}));
    
    for j = 1 : length(sheets)
        ind_sheet = find(strcmp(sheets{j}, C{2}(ind_file)));
        
        [data.txt, data.num, data.raw] = xlsread([folder 'Archive/' file_names{i}],sheets{j});
        
        
                
        ind_plate = strfind('ABCDEFGHJIKLMNOP', C{3}{ind_file(ind_sheet)});
        ind_row = strfind('ABCDEFGHJIKLMNOP', C{4}{ind_file(ind_sheet)});
        ind_col = strfind('ABCDEFGHJIKLMNOP', C{5}{ind_file(ind_sheet)});
        ind_data = strfind('ABCDEFGHJIKLMNOP', C{6}{ind_file(ind_sheet)});
        
        if isempty(ind_plate) || isempty(ind_row) || isempty(ind_col)
            fprintf('%s\t%s\n', file_names{i}, sheets{j});
        else
            all_data.plate = [all_data.plate; data.raw(:,ind_plate)];
            all_data.row = [all_data.row; data.raw(:,ind_row)];
            all_data.col = [all_data.col; data.raw(:,ind_col)];
            all_data.data = [all_data.data; data.raw(:,end-2)];
        end
    end
end

all_data.plate = cell2mat(all_data.plate);
all_data.col = cell2mat(all_data.col);

inds = isnan(all_data.plate) | isnan(all_data.col) | cellfun(@isnumeric, all_data.row) | ~cellfun(@isnumeric, all_data.data);

all_data.plate(inds) = [];
all_data.row(inds) = [];
all_data.col(inds) = [];
all_data.data(inds) = [];

all_data.data = cell2mat(all_data.data);


% Map ORFs
all_data.orf = cell(size(all_data.data));
for i = 1 : length(all_data.plate)
    inds = find(map.plate == all_data.plate(i) & map.col == all_data.col(i) & strcmp(all_data.row{i}, map.row));
    if ~isempty(inds)
        all_data.orf(i,1) = map.orf(inds);
    end
end

inds = find(cellfun(@isempty, all_data.orf));
all_data.orf(inds) = [];
all_data.plate(inds) = [];
all_data.row(inds) = [];
all_data.col(inds) = [];
all_data.data(inds) = [];


% Average multiple occurences of the same ORF
[all_data.orf_u, all_data.data_avg] = grpstats(all_data.data,all_data.orf,{'gname','mean'});

freimoser_amrhein_2006.orfs = all_data.orf_u;
freimoser_amrhein_2006.data = all_data.data_avg;
freimoser_amrhein_2006.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2006_Freimoser~Amrhein/freimoser_amrhein_2006.mat','freimoser_amrhein_2006');

% Save data into database
dt = freimoser_amrhein_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Copic~Miller, 2009
% DATA = copic_miller_2009
copic_miller_2009.pmid = 19433630;

phenotypes = {'Kar2 secretion'};
treatments = {''};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Copic~Miller/mat_alpha_obs_v1.0.xlsx','mat_alpha_obs');
tested_orfs = tested.raw(2:end,1);

inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Copic~Miller/TableS1.xlsx','TABLE S1');
hits_orfs = data.raw(5:end,1);
hits_notes = data.raw(5:end,3);
hits_scores = data.raw(5:end,6);

inds = strcmp('mat-a only', hits_notes);
hits_orfs(inds) = [];
hits_scores(inds) = [];

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_scores(inds) = [];
hits_orfs = upper(hits_orfs);

hits_scores = cell2mat(hits_scores);

copic_miller_2009.orfs = tested_orfs;
copic_miller_2009.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
copic_miller_2009.data(ind2,1) = hits_scores(ind1);

copic_miller_2009.ph = [strcat(phenotypes, '; ', treatments)];

save('Datasets/Phenotypes/2009_Copic~Miller/copic_miller_2009.mat','copic_miller_2009');

% Save data into database
dt = copic_miller_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Kanki~Klionsky, 2009
% DATA = kanki_klionsky_2009
kanki_klionsky_2009.pmid = 19793921;

phenotypes = {'mitophagy'};
treatments = {'YPL'};

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Kanki~Klionsky/TableS2-2.xlsx','Sheet1');
hits_orfs = data.raw(4:end,1);
hits_scores = data.raw(4:end,2);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_scores(inds) = [];
hits_orfs = upper(hits_orfs);

inds = strcmp('x', hits_scores);
hits_scores(inds) = {NaN};
inds = strcmp('No', hits_scores);
hits_scores(inds) = {-1};
inds = strcmp('Yes', hits_scores);
hits_scores(inds) = {1};

hits_scores = cell2mat(hits_scores);

hits_orfs_u = unique(hits_orfs);
hits_scores_u = zeros(length(hits_orfs_u),1)+NaN;
for i = 1 : length(hits_orfs_u)
    inds = find(strcmp(hits_orfs_u{i}, hits_orfs));
    if length(inds) == 1
        hits_scores_u(i) = hits_scores(inds);
    elseif length(inds) == 2
        tmp = prod(hits_scores(inds));
        if isnan(tmp) || tmp < 0    % If at least one of the values is NaN or the values are of opposite sign, set them to NaN
            hits_scores_u(i) = NaN;
        else
            hits_scores_u(i) = hits_scores(inds(1));    % If the values are of the same sign, pick the first number (they are the same)
        end
    else
        fprintf('%s\t', hits_orfs_u{i});
        for j = 1 : length(inds)
            fprintf('%d,', hits_scores(inds(j)));
        end
        fprintf('\n');
    end
end

kanki_klionsky_2009.orfs = hits_orfs_u;
kanki_klionsky_2009.data = hits_scores_u;

kanki_klionsky_2009.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2009_Kanki~Klionsky/kanki_klionsky_2009.mat','kanki_klionsky_2009');

% Save data into database
dt = kanki_klionsky_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Hancock~Lopes, 2006
% DATA = hancock_lopes_2006
hancock_lopes_2006.pmid = 16582425;

phenotypes = {'Opi-'};
treatments = {''};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Hancock~Lopes/mat_alpha_061101.xlsx','mat_alpha_061101');
tested_orfs = tested.raw(4:end,2);

inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));


% Load data
fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Hancock~Lopes/hits_genenames.txt');
hits_genenames = textscan(fid, '%s');
hits_genenames = hits_genenames{1};
fclose(fid);

hits_orfs = genename2orf(hits_genenames,'noannot');

hits_scores = ones(length(hits_orfs),1);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_scores(inds) = [];
hits_orfs = upper(hits_orfs);

hancock_lopes_2006.orfs = tested_orfs;
hancock_lopes_2006.data = zeros(length(tested_orfs), length(phenotypes))+NaN;
[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
hancock_lopes_2006.data(ind2,1) = hits_scores(ind1,1);

hancock_lopes_2006.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2006_Hancock~Lopes/hancock_lopes_2006.mat','hancock_lopes_2006'); 

% Save data into database
dt = hancock_lopes_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Alamgir~Golshani, 2008
% DATA = alamgir_golshani_2008
alamgir_golshani_2008.pmid = 19055778;

phenotypes = {'growth [colony size]'};
treatments = {'paromomycin'};

% Load tested
[map.txt, map.num, map.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Alamgir~Golshani/Master plate list -- part 1.xlsx','Sheet2-Corrected');
map.orf = map.raw(2:end,2);
map.plate = cell2mat(map.raw(2:end,4));
map.row = cell2mat(map.raw(2:end,5));
map.col = cell2mat(map.raw(2:end,6));


[map2.txt, map2.num, map2.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Alamgir~Golshani/Master plate list -- part2.xlsx','Sheet2-Corrected');
[rinds,cinds] = find(~cellfun(@isnumeric, map2.raw(:,4:6)));
map2.raw(unique(rinds),:) = [];
map2.orf = map2.raw(:,2);
map2.plate = cell2mat(map2.raw(:,4));
map2.row = cell2mat(map2.raw(:,5));
map2.col = cell2mat(map2.raw(:,6));

map.orf = [map.orf; map2.orf];
map.plate = [map.plate; map2.plate];
map.row = [map.row; map2.row];
map.col = [map.col; map2.col];

map.inds = sub2ind([length(unique(map.plate)), length(unique(map.row)), length(unique(map.col))], map.plate, map.row, map.col);

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Alamgir~Golshani/Plate_Analyzer_Alamgir_13mgmLParomomycin_08Dec05(1).xlsx');
data.plate = cell2mat(data.raw(2:end,1));
data.row = cell2mat(data.raw(2:end,2));
data.col = cell2mat(data.raw(2:end,3));
data.data = cell2mat(data.raw(2:end,4));

data.inds = sub2ind([length(unique(data.plate)), length(unique(data.row)), length(unique(data.col))], data.plate, data.row, data.col);
data.orf = cell(length(data.inds),1);

[~,ind1,ind2] = intersect(data.inds, map.inds);     
data.orf(ind1) = map.orf(ind2);

inds = find(cellfun(@isempty, data.orf));       % 26 positions present in data don't have an ORF
data.orf(inds) = [];
data.data(inds) = [];

inds = cellfun(@isnumeric, data.orf);
data.orf(inds) = [];
data.data(inds) = [];

data.orf = cellfun(@strtrim, data.orf,'UniformOutput',0);

inds = ~strncmp('Y', data.orf,1);
data.orf(inds) = [];
data.data(inds) = [];
data.orf = upper(data.orf);

[data_avg.orf, data_avg.data] = grpstats(data.data,data.orf,{'gname','mean'});

alamgir_golshani_2008.orfs = data_avg.orf;
alamgir_golshani_2008.data = data_avg.data;

alamgir_golshani_2008.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2008_Alamgir~Golshani/alamgir_golshani_2008.mat','alamgir_golshani_2008'); 

% Save data into database
dt = alamgir_golshani_2008;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Luban~Schmidt, 2005
% DATA = luban_schmidt_2005
luban_schmidt_2005.pmid = 15908144;

phenotypes = {'petite','mtDNA intron slicing'};
treatments = {'Gly'};

% Load tested
fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Luban~Schmidt/list_of_used_knockouts_PhD_Thesis_Luban.txt');
C = textscan(fid, '%s\n');
fclose(fid);
tested_orfs = unique(upper(C{1}));

% Load data
fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Luban~Schmidt/list_of_pet_mutants.txt');
C = textscan(fid,'%s\n');
fclose(fid);
pet_mutants = unique(upper(C{1}));

fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Luban~Schmidt/list_of_intron_def_mutants.txt');
C = textscan(fid,'%s\n');
fclose(fid);
intron_mutants = unique(upper(C{1}));

inds = find(cellfun(@isempty, tested_orfs));

inds = find(cellfun(@isnumeric, tested_orfs));


tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', tested_orfs,1));

luban_schmidt_2005.orfs = tested_orfs;
luban_schmidt_2005.data = nan(length(tested_orfs),2);

missing = setdiff(pet_mutants, tested_orfs);

[~,ind1,ind2] = intersect(tested_orfs, pet_mutants);
luban_schmidt_2005.data(ind1,1) = 1;
luban_schmidt_2005.data(isnan(luban_schmidt_2005.data(:,1)),1) = 0;

luban_schmidt_2005.data(ind1,2) = 0;
[~,ind1,ind2] = intersect(tested_orfs, intron_mutants);
luban_schmidt_2005.data(ind1,2) = -1;

luban_schmidt_2005.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2005_Luban~Schmidt/luban_schmidt_2005.mat','luban_schmidt_2005'); 

% Save data into database
dt = luban_schmidt_2005;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Nyswaner~Garfinkel,2008
% DATA = nyswaner_garfinkel_2008
nyswaner_garfinkel_2008.pmid = 18202368;

phenotypes = {'increased Ty1 transposon mobility'};
treatments = {''};

% Load tested
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Nyswaner~Garfinkel/Matalphakos counted.xlsx','Sheet2');
tested_orfs = tested.raw(6:end,2);

inds = find(cellfun(@isempty, tested_orfs));
tested_orfs(inds) = [];

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Nyswaner~Garfinkel/genetics.107.082602-9.xlsx','Sheet1');
hits_orfs = data.raw(:,1);

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];

nyswaner_garfinkel_2008.orfs = tested_orfs;
nyswaner_garfinkel_2008.data = zeros(length(tested_orfs),1);

missing = setdiff(hits_orfs, tested_orfs);

% Adding 3 ORFs to the list of tested
tested_orfs = [tested_orfs; missing];

[~,ind1,ind2] = intersect(tested_orfs, hits_orfs);
nyswaner_garfinkel_2008.data(ind1,1) = 1;

nyswaner_garfinkel_2008.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2008_Nyswaner~Garfinkel/nyswaner_garfinkel_2008.mat','nyswaner_garfinkel_2008'); 

% Save data into database
dt = nyswaner_garfinkel_2008;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Dakshinamurthy~Garfinkel, 2010
% DATA = dakshinamurthy_garfinkel_2010
dakshinamurthy_garfinkel_2010.pmid = 20498295;

phenotypes = {'decreased Ty1 transposon mobility'};
treatments = {''};

load nyswaner_garfinkel_2008;
dakshinamurthy_garfinkel_2010.orfs = nyswaner_garfinkel_2008.orfs;      % Same screen, positive & negative phenotypes analyzed separately
dakshinamurthy_garfinkel_2010.data = zeros(length(dakshinamurthy_garfinkel_2010.orfs),1);

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2010_Dakshinamurthy~Garfinkel/TableS1-2.xlsx','Sheet2');
hits_orfs = data.raw(:,2);
hits_data = data.raw(:,3);

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds) = [];

hits_data = cellfun(@length, hits_data);

missing = setdiff(hits_orfs, dakshinamurthy_garfinkel_2010.orfs);

% Adding 4 ORFs to the list of tested
dakshinamurthy_garfinkel_2010.orfs = [dakshinamurthy_garfinkel_2010.orfs; missing];

[~,ind1,ind2] = intersect(dakshinamurthy_garfinkel_2010.orfs, hits_orfs);
dakshinamurthy_garfinkel_2010.data(ind1,1) = hits_data(ind2);

dakshinamurthy_garfinkel_2010.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2010_Dakshinamurthy~Garfinkel/dakshinamurthy_garfinkel_2010.mat','dakshinamurthy_garfinkel_2010'); 

% Save data into database
dt = dakshinamurthy_garfinkel_2010;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Lam~Conibear, 2006
% DATA = lam_conibear_2006
% TESTED = not available
lam_conibear_2006.pmid = 16818716;

phenotypes = {'polytopic membrane protein trafficking'};
treatments = {''};

% Load data
fid = fopen('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Lam~Conibear/hits_genes_data.txt');
C = textscan(fid, '%s\t%.3f\n');
fclose(fid);

hits_genes = C{1};
hits_data = C{2};

hits_orfs = genename2orf(hits_genes,'noannot');
hits_orfs(strcmp('CHS4', hits_orfs)) = {'YBL061C'};

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds) = [];

lam_conibear_2006.orfs = hits_orfs;
lam_conibear_2006.data = hits_data;

lam_conibear_2006.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2006_Lam~Conibear/lam_conibear_2006.mat','lam_conibear_2006'); 

% Save data into database
dt = lam_conibear_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Burston~Conibear, 2009
% DATA = burston_conibear_2009
% TESTED = not available
burston_conibear_2009.pmid = 19506040;

phenotypes = {'endocytosis (MatA)';'endocytosis (MatAlpha)'};
treatments = {''};

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Burston~Conibear/JCB_200811116_TS1.xlsx','TableS1');

hits_orfs = data.raw(6:end,2);
hits_data_a = data.raw(6:end,4);
hits_data_alpha = data.raw(6:end,5);

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data_a(inds) = [];
hits_data_alpha(inds) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data_a(inds) = [];
hits_data_alpha(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data_a(inds) = [];
hits_data_alpha(inds) = [];

burston_conibear_2009.orfs = hits_orfs;
burston_conibear_2009.data = cell2mat([hits_data_a hits_data_alpha]);

burston_conibear_2009.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2009_Burston~Conibear/burston_conibear_2009.mat','burston_conibear_2009'); 

% Save data into database
dt = burston_conibear_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Kwak~Weller, 2011
% DATA = kwak_weller_2011
kwak_weller_2011.pmid = 21193664;

phenotypes = {'growth (spot assay)';'growth (culture turbidity)'};
treatments = {'2,4-DAPG [200 ug/ml]';'2,4-DAPG [40 ug/ml]'};

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Kwak~Weller/Yeast liquid screen OD.xlsx','Set 1');

hits_orfs1 = data.raw(2:end,5);
hits_data1 = data.raw(2:end,8);

[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Kwak~Weller/Yeast liquid screen OD.xlsx','Set 2');

hits_orfs2 = data.raw(2:end,5);
hits_data2 = data.raw(2:end,8);


[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Kwak~Weller/Yeast liquid screen OD.xlsx','Sheet3');

hits_orfs3 = data.raw(2:end,5);
hits_data3 = data.raw(2:end,8);

hits_orfs = [hits_orfs1; hits_orfs2; hits_orfs3];
hits_data = [hits_data1; hits_data2; hits_data3];

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

inds = find(cellfun(@isnumeric, hits_orfs) | ~cellfun(@isnumeric, hits_data));
hits_orfs(inds) = [];
hits_data(inds) = [];

hits_data = cell2mat(hits_data);

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

hits_orfs = upper(hits_orfs);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds) = [];

[hits_orfs_a, hits_data_a] = grpstats(hits_data, hits_orfs, {'gname','mean'});

kwak_weller_2011.orfs = hits_orfs_a;
kwak_weller_2011.data = hits_data_a;

% Second dataset
[data.num, data.txt, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Kwak~Weller/CANDIDATE VARIFY.xlsx','Sheet1');
hits_orfs = data.raw(3:end,5);
inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
hits_orfs = unique(upper(hits_orfs));
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];

[missing, ix] = setdiff(hits_orfs, kwak_weller_2011.orfs);

[~,ind1,ind2] = intersect(hits_orfs, kwak_weller_2011.orfs);
kwak_weller_2011.data(ind2,2) = -1;

kwak_weller_2011.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2011_Kwak~Weller/kwak_weller_2011.mat','kwak_weller_2011'); 

% Save data into database
dt = kwak_weller_2011;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Abe~Minegishi, 2008
% DATA = abe_minegishi_2008
abe_minegishi_2008.pmid = 18245339;

phenotypes = {'growth (culture turbidity)'};
treatments = {'high pressure';'low temperature'};

% Load tested strains
[tested.txt, tested.num, tested.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Abe~Minegishi/mat_alpha_041902.xlsx','mat_alpha_041902.txt');

tested_orfs = tested.raw(4:end,3);

inds = find(cellfun(@isempty, tested_orfs));
tested_orfs(inds) = [];

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(upper(cellfun(@strtrim, tested_orfs,'UniformOutput',0)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

[data.num, data.txt, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Abe~Minegishi/Table1_Abe_Genetics.xlsx','Table 1 (2)');

hits_orfs = data.raw(7:end,3);
hits_data = data.raw(7:end,[26 31]);
inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
hits_orfs = unique(upper(hits_orfs));
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_data = cell2mat(hits_data);
hits_data = hits_data/100;  % transform percent into fractions

[missing, ix] = setdiff(hits_orfs, tested_orfs);

abe_minegishi_2008.orfs = tested_orfs;
abe_minegishi_2008.data = nan(length(tested_orfs),2);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
abe_minegishi_2008.data(ind2,:) = hits_data(ind1,:);

abe_minegishi_2008.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2008_Abe~Minegishi/abe_minegishi_2008.mat','abe_minegishi_2008'); 

% Save data into database
dt = abe_minegishi_2008;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Begley~Samson, 2004
% DATA = begley_samson_2004
begley_samson_2004.pmid = 15469827;

phenotypes = {'growth (spot assay)'};
treatments = {'MMS';'t-BuOOH';'4NQO';'UV score'};

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2004_Begley~Samson/Begley2003.xlsx','Sheet1');

hits_orfs = data.raw(2:end,1);
hits_data = data.raw(2:end,7:10);

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_orfs = upper(cellfun(@strtrim, hits_orfs,'UniformOutput',0));

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_data = -cell2mat(hits_data);   % The sensitivity scores are such that the higher the number, the more sensitive the mutant.

[t,t2] = grpstats(hits_data,hits_orfs,{'gname','mean'});
begley_samson_2004.orfs = t;
begley_samson_2004.data = t2;

begley_samson_2004.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2004_Begley~Samson/begley_samson_2004.mat','begley_samson_2004'); 

% Save data into database
dt = begley_samson_2004;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Lesuisse~Dancis, 2005
% DATA = lesuisse_dancis_2005
lesuisse_dancis_2005.pmid = 15489514;

% Part 1
phenotypes = {'iron updake','ferrireductase activity','Fe(II) uptake'};
treatments = {'standard'};

% Load data
[data.txt, data.num, data.raw] = ...
    xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Lesuisse~Dancis/EUROSCARF_SCREENING.xlsx','Sheet1');

hits_orfs = data.raw(2:end,1);
hits_data = data.raw(2:end,[16 12 14]);

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_orfs = upper(strtrim(hits_orfs));

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

inds = find(~cellfun(@isnumeric, hits_data));
hits_data(inds) = {NaN};

hits_data = cell2mat(hits_data);

% % Part 2
% [data.txt, data.num, data.raw] = ...
%     xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2005_Lesuisse~Dancis/Final_1screen-1.xlsx','Tabelle1');
% hits_orfs = data.raw(2:end,2);
% hits_data = data.raw(2:end,11:21);
% 
% inds = find(cellfun(@isempty, hits_orfs));
% hits_orfs(inds) = [];
% hits_data(inds,:) = [];
% 
% inds = find(cellfun(@isnumeric, hits_orfs));
% hits_orfs(inds) = [];
% hits_data(inds,:) = [];
% 
% hits_orfs = upper(strtrim(hits_orfs));
% 
% inds = find(~strncmp('Y', hits_orfs,1));
% hits_orfs(inds) = [];
% hits_data(inds,:) = [];
% 
% inds = find(cellfun(@isnumeric, hits_data));
% hits_data(inds) = {''};


% Combine part1 and part2
[t,t2] = grpstats(hits_data,hits_orfs,{'gname','mean'});
lesuisse_dancis_2005.orfs = t;
lesuisse_dancis_2005.data = t2;

lesuisse_dancis_2005.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2005_Lesuisse~Dancis/lesuisse_dancis_2005.mat','lesuisse_dancis_2005'); 

% % TODO: Save data into database
% dt = lesuisse_dancis_2005;
% datasets = get_datasets_for_paper(dt);
% 
% [~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
% [~,ph_ix] = sort(dt.ph);
% 
% % Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
% datasets.names(database_ix,:)
% dt.ph(ph_ix)
% 
% insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Alto~Dixon, 2006
% DATA = alto_dixon_2006
alto_dixon_2006.pmid = 16413487;

% Part 1
phenotypes = {'growth'};
treatments = {'IpgB2 effector protein'};

% Load data
hits_gn = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Alto~Dixon/hits.txt','%s');
hits_data = ones(size(hits_gn));

hits_orfs = genename2orf_sgd(hits_gn);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

% Load tested
tested_orfs = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2006_Alto~Dixon/FG_array_genes.txt','%s');

[missing, ix] = setdiff(hits_orfs, tested_orfs);

alto_dixon_2006.orfs = tested_orfs;
alto_dixon_2006.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
alto_dixon_2006.data(ind2,:) = hits_data(ind1,:);

alto_dixon_2006.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2006_Alto~Dixon/alto_dixon_2006.mat','alto_dixon_2006'); 

% TODO: Save data into database
dt = alto_dixon_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Merz~Westermann, 2009
% DATA = merz_westermann_2009
merz_westermann_2009.pmid = 19751518;

% Part 1
phenotypes = {'growth'};
treatments = {'Gly, 3%'};

% Load data
[data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Merz~Westermann/gb-2009-10-9-r95-s1.xlsx');
hits_orfs = data.raw(3:end,1);
hits_orfs = unique(strtrim(upper(hits_orfs)));

hits_data = -ones(size(hits_orfs));

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

% Load tested
[tested.txt, tested.num, tested.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2009_Merz~Westermann/pet-Screen.xlsx','mat_alpha_obs');
tested_orfs = tested.raw(2:end,2);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

merz_westermann_2009.orfs = tested_orfs;
merz_westermann_2009.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
merz_westermann_2009.data(ind2,:) = hits_data(ind1,:);

merz_westermann_2009.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2009_Merz~Westermann/merz_westermann_2009.mat','merz_westermann_2009'); 

% Save data into database
dt = merz_westermann_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Ralser~Lehrach, 2008
% DATA = ralser_lehrach_2008
ralser_lehrach_2008.pmid = 19004802;

phenotypes = {'growth'};
treatments = {'2-DG'};

% Load data
hits_orfs = textread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Ralser~Lehrach/hits_orfs.txt','%s');
hits_orfs = unique(strtrim(upper(hits_orfs)));

hits_data = ones(size(hits_orfs));

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

% Load tested
[tested.txt, tested.num, tested.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2008_Ralser~Lehrach/Mat_a_obs_v2.0.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

ralser_lehrach_2008.orfs = tested_orfs;
ralser_lehrach_2008.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
ralser_lehrach_2008.data(ind2,:) = hits_data(ind1,:);

ralser_lehrach_2008.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2008_Ralser~Lehrach/ralser_lehrach_2008.mat','ralser_lehrach_2008'); 

% Save data into database
dt = ralser_lehrach_2008;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

%% Gaupel~Tenniswood, 2014
% DATA = gaupel_tenniswood_2004
gaupel_tenniswood_2004.pmid = 24968945;

phenotypes = {'growth'};
treatments = {'CG-1521'};

% Load data (1)
[data_sens.txt, data_sens.num, data_sens.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2014_Gaupel~Tenniswood/1471-2164-15-528-s1.xlsx','Top sensitive strains');

hits_orfs1 = data_sens.raw(2:end,1);
hits_data1 = -cell2mat(data_sens.raw(2:end,3));

hits_orfs1 = strtrim(upper(hits_orfs1));

inds = find(~strncmp('Y', hits_orfs1,1));
hits_orfs1(inds) = [];
hits_data1(inds,:) = [];

[t,t2] = grpstats(hits_data1,hits_orfs1,{'gname','mean'});
hits_orfs1 = t;
hits_data1 = t2;

% Load data (2)
[data_res.txt, data_res.num, data_res.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2014_Gaupel~Tenniswood/1471-2164-15-528-s1.xlsx','Top resistant strains');

hits_orfs2 = data_res.raw(2:end,1);
hits_data2 = -cell2mat(data_res.raw(2:end,3));

hits_orfs2 = strtrim(upper(hits_orfs2));

inds = find(~strncmp('Y', hits_orfs2,1));
hits_orfs2(inds) = [];
hits_data2(inds,:) = [];

[t,t2] = grpstats(hits_data2,hits_orfs2,{'gname','mean'});
hits_orfs2 = t;
hits_data2 = t2;

% Eliminate overlap between sensitive and resistant strains?
[~,ind1,ind2] = intersect(hits_orfs1, hits_orfs2);
hits_orfs1(ind1) = []; hits_data1(ind1,:) = [];
hits_orfs2(ind2) = []; hits_data2(ind2,:) = [];

% Combine part1 and part2
hits_orfs = [hits_orfs1; hits_orfs2];
hits_data = [hits_data1; hits_data2];

% Load tested
[tested.txt, tested.num, tested.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2014_Gaupel~Tenniswood/CompleteDeletionLibrary.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

gaupel_tenniswood_2004.orfs = tested_orfs;
gaupel_tenniswood_2004.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
gaupel_tenniswood_2004.data(ind2,:) = hits_data(ind1,:);

gaupel_tenniswood_2004.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2014_Gaupel~Tenniswood/gaupel_tenniswood_2004.mat','gaupel_tenniswood_2004'); 

% Save data into database
dt = gaupel_tenniswood_2004;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


%% Suzuki~Yoshida, 2011
% DATA = suzuki_yoshida_2011
suzuki_yoshida_2011.pmid = 21601516;

phenotypes = {'glutathione abundance'};
treatments = {'standard'};

% Load data
[data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Suzuki~Yoshida/hits_data.xlsx');

hits_genes = data.raw(1:end,1);
hits_data = cell2mat(data.raw(1:end,2));

hits_genes = strtrim(lower(hits_genes));
hits_orfs = genename2orf_sgd(hits_genes);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

[t,t2] = grpstats(hits_data,hits_orfs,{'gname','mean'});
hits_orfs = t;
hits_data = t2;

% Load tested
[tested.txt, tested.num, tested.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2011_Suzuki~Yoshida/YKOmatalpha_GSH_list070508.xlsx');
tested_orfs = tested.raw(3:end,3);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

suzuki_yoshida_2011.orfs = tested_orfs;
suzuki_yoshida_2011.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
suzuki_yoshida_2011.data(ind2,:) = hits_data(ind1,:);

suzuki_yoshida_2011.ph = strcat(phenotypes, '; ', treatments);

save('Datasets/Phenotypes/2011_Suzuki~Yoshida/suzuki_yoshida_2011.mat','suzuki_yoshida_2011'); 

% Save data into database
dt = suzuki_yoshida_2011;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));