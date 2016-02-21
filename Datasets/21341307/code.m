%% Yoshikawa~Shimizu, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};


yoshikawa_shimizu_2011.source = {'http://onlinelibrary.wiley.com/store/10.1002/yea.1843/asset/supinfo/yea_1843_supportinginforTS1.xls?v=1&s=ebd728c627c24d863139e83ed47ca3e94b22c39e'};
yoshikawa_shimizu_2011.downloaddate = {'2013-03-08'};
yoshikawa_shimizu_2011.pmid = 21341307;

phenotypes = {'Growth, exponential growth rate (h^-1)'};
treatments = {''};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/yea_1843_supportinginforTS1.xlsx', 'Deletion');

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

save('./yoshikawa_shimizu_2011.mat','yoshikawa_shimizu_2011');

fid = fopen('./yoshikawa_shimizu_2011.txt','w');
write_matrix_file(fid, yoshikawa_shimizu_2011.orfs, yoshikawa_shimizu_2011.ph, yoshikawa_shimizu_2011.data);
fclose(fid);

end
