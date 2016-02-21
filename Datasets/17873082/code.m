%% Botet~Santos, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
% NOTES = data unnormalized; potentially batch, row/col normalization needed?


botet_santos_2007.source = {'Javier Botet'};
botet_santos_2007.downloaddate = {'2013-04-18'};
botet_santos_2007.pmid = 17873082;

phenotypes = {'Growth, OD'};
treatments = {'SMM, 77 h'; 'SMM, 120 h'; 'Sulfanilamide, 0.1 mg/ml, 77 h';'Sulfanilamide, 0.1 mg/ml, 120 h'};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/1ScreenSULFA&MS&MS+PABA.xlsx', 'DATA');

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

save('./botet_santos_2007.mat','botet_santos_2007');

fid = fopen('./botet_santos_2007.txt','w');
write_matrix_file(fid, botet_santos_2007.orfs, botet_santos_2007.ph, botet_santos_2007.data);
fclose(fid);

end
