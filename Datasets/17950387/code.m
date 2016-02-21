%% Galvan~Smith, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

galvan_smith_2008.pmid = 17950387;

phenotypes = {'growth'};
treatments = {'chimaphilin'};

% Load plate maps
[FILENAMES{end+1}, map.raw] = read_data('xlsread','./raw_data/yGDA-Master_Plate_list_Combined(New).xlsx', 'Sheet1');
map.raw(1,:) = [];
map.platerowcol = cell2mat(map.raw(:,4:6));

% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/IMELDA 08Feb2006GDAraw data.xls', 'Sheet1');
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

save('./galvan_smith_2008.mat','galvan_smith_2008');

fid = fopen('./galvan_smith_2008.txt','w');
write_matrix_file(fid, galvan_smith_2008.orfs, galvan_smith_2008.ph, galvan_smith_2008.data);
fclose(fid);

end
