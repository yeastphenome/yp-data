%% Warringer~Blomberg, 2003
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

warringer_blomberg_2003.pmid = 14676322;

phenotypes = {'Growth, lag phase';'Growth, exponential growth rate';'Growth, saturation level'};
treatments = {'NaCl, 0.85 M'};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/ORIG130305_LPI NaCl.xlsx', 'LPI');

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

save('./warringer_blomberg_2003.mat','warringer_blomberg_2003');

fid = fopen('./warringer_blomberg_2003.txt','w');
write_matrix_file(fid, warringer_blomberg_2003.orfs, warringer_blomberg_2003.ph, warringer_blomberg_2003.data);
fclose(fid);

end

