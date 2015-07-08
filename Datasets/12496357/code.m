%% Begley~Samson, 2002
% DATA = begley_samson_2002
function FILENAMES = code()
FILENAMES = {};
% NOTE = data sorting fixed (lower=slower)

begley_samson_2002.pmid = 12496357;

phenotype = {'Growth, dilution spot intensity'};

treatments = {'4NQO','MMS','t-BuOOH','UV'};
doses(1,:) = {'0.1 ug/ml','0.2 ug/ml','0.3 ug/ml','0.4 ug/ml'};
doses(2,:) = {'0.01%','0.02%','0.025%','0.03%'};
doses(3,:) = {'0.25 mM','0.5 mM','0.75 mM','1 mM'};
doses(4,:) = {'25 J/m2','50 J/m2','75 J/m2','100 J/m2'};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/ORIG130404_Begley2001raw.xlsx', 'Sheet1');

% Find the column with the systematic ORF name
ind_orf = find(strcmp('ORF', data.raw(1,:)));

% Find the column with the name of the condition
ind_condition = find(strcmp('Treatment', data.raw(1,:)));

% Find the columns with the relevant data
% This is the ratio of the growth of this spot on the untreated plate to
% the growth on a treated plate divided by the same ratio for the highest
% growing WT colony. In the paper, a value of greater than 1.5 was
% considered sensitive; a value of less than .66 was considered resistant.

ind_data1 = find(strcmp('exp 1 highest control/dose 1', data.raw(1,:))); ind_data1 = [ind_data1:ind_data1+3];
ind_data2 = find(strcmp('exp 1 lowest control/dose 1', data.raw(1,:))); ind_data2 = [ind_data2:ind_data2+3];

ind_data3 = find(strcmp('exp 2 highest control/dose 1', data.raw(1,:))); ind_data3 = [ind_data3:ind_data3+3];
ind_data4 = find(strcmp('exp 2 lowest control/dose 1', data.raw(1,:))); ind_data4 = [ind_data4:ind_data4+3];

ind_data5 = find(strcmp('exp 3 highest control/dose 1', data.raw(1,:))); ind_data5 = [ind_data5:ind_data5+3];
ind_data6 = find(strcmp('exp 3 lowest control/dose 1', data.raw(1,:))); ind_data6 = [ind_data6:ind_data6+3];

% Extract the most relevant information, ignore the rest.
all_strings = cellfun(@num2str, data.raw(:,ind_orf),'UniformOutput',0);
ind_not_empty = setdiff(1:size(data.raw,1), find(strcmp('NaN', all_strings)))';

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


a = mfilename('fullpath');
a = a(1:end-4);
save([a,'begley_samson_2002.mat'],'begley_samson_2002');
return;


% Save data into database
dt = begley_samson_2002;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

