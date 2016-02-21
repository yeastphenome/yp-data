%% Hartman~Tippery, 2004
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

hartman_tippery_2004.pmid = 15239834;
hartman_tippery_2004.desc = {'The loaded values are Z-scores with respect to the wild-type in a given condition (UNT or HU). The final values (normalized to UNT) are the difference between treated and untreated z-scores.'};

phenotypes = {'Growth, AUGC'};
treatments = {'UNT';'HU, 50 mM';'HU, 150 mM'};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/gb-2004-5-7-r49-s7.xlsx', 'data');

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

save('./hartman_tippery_2004.mat','hartman_tippery_2004');

fid = fopen('./hartman_tippery_2004.txt','w');
write_matrix_file(fid, hartman_tippery_2004.orfs, hartman_tippery_2004.ph, hartman_tippery_2004.data);
fclose(fid);

end
