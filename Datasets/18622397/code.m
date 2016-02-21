%% Breslow~Weissman, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

breslow_weissman_2008.source = {'http://www.nature.com/nmeth/journal/v5/n8/extref/nmeth.1234-S4.xls'};
breslow_weissman_2008.downloaddate = {'2013-03-12'};
breslow_weissman_2008.pmid = 18622397;

phenotypes = {'Growth, flow cytometry'};
treatments = {''};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/nmeth.1234-S5.xlsx', 'Deletion growth rates');

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

save('./breslow_weissman_2008.mat','breslow_weissman_2008');

fid = fopen('./breslow_weissman_2008.txt','w');
write_matrix_file(fid, breslow_weissman_2008.orfs, breslow_weissman_2008.ph, breslow_weissman_2008.data);
fclose(fid);

end
