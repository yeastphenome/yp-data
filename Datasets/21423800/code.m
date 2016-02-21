%% Brett~Rao, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

brett_rao_2011.source = {'http://www.plosone.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pone.0017619.s003'};
brett_rao_2011.downloaddate = {'2013-07-18'};
brett_rao_2011.pmid = 21423800;

phenotypes = {'Growth, OD600';'Vacuolar pH (pHv)'};
treatments = {'external pH 2.7';'external pH 4.0';'external pH 7.0'};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/journal.pone.0017619.s003.xlsx', 'Unsorted Data');

% Get indices of the data columns
ind_data = [4:6 8:10];

% Eliminate anything that doesn't have an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,2) = cellfun(@strtrim, data.raw(:,2),'UniformOutput',0);

% A couple of manual fixes
data.raw(strcmp(data.raw(:,2),'YLR228 C'),2) = {'YLR228C'};
data.raw(strcmp(data.raw(:,2), 'YML009c'),2) = {'YML009C'};
data.raw(strcmp(data.raw(:,2), 'YMR062 C'),2) = {'YMR062C'};

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isempty, regexp(data.raw(:,2), 'Y[A-P][RL][0-9]{3}[CW](-[ABC])*')));
data.raw(inds,:) = [];

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

save('./brett_rao_2011.mat','brett_rao_2011');

fid = fopen('./brett_rao_2011.txt','w');
write_matrix_file(fid, brett_rao_2011.orfs, brett_rao_2011.ph, brett_rao_2011.data);
fclose(fid);

end
