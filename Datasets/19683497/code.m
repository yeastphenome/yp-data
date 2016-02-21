%% Fillingham~Andrews, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};


fillingham_andrews_2009.source = {'http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S1097276509004614/1-s2.0-S1097276509004614-mmc3.xls/272198/FULL/S1097276509004614/c47432a0acc5d8c66fa4ae36e559c4f9/mmc3.xls'};
fillingham_andrews_2009.downloaddate = {'2013-03-07'};
fillingham_andrews_2009.pmid = 19683497;

phenotypes = {'HTA1 expression, z-score of log2 GFP:RFP'};
treatments = {''};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/mmc3.xlsx', 'Sheet1');

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

save('./fillingham_andrews_2009.mat','fillingham_andrews_2009');

fid = fopen('./fillingham_andrews_2009.txt','w');
write_matrix_file(fid, fillingham_andrews_2009.orfs, fillingham_andrews_2009.ph, fillingham_andrews_2009.data);
fclose(fid);

end
