%% Galvan Marquez~Smith, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};


galvan_marquez_smith_2013.source = {'Imelda Galvan Marquez'};
galvan_marquez_smith_2013.downloaddate = {'2013-05-19'};
galvan_marquez_smith_2013.pmid = 23624539;

phenotypes = {'Growth, colony size'};
treatments = {'Chitosan'};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Chitosan Effect on GDA (raw data). Exp. 1 to 3. Imelda Galvan, 2013-1.xlsx', 'Sheet1');

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

save('./galvan_marquez_smith_2013.mat','galvan_marquez_smith_2013');

fid = fopen('./galvan_marquez_smith_2013.txt','w');
write_matrix_file(fid, galvan_marquez_smith_2013.orfs, galvan_marquez_smith_2013.ph, galvan_marquez_smith_2013.data);
fclose(fid);

end
