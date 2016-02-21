%% Peyroche~Plateau, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};


peyroche_plateau_2012.source = {'http://www.plosone.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pone.0036343.s004'};
peyroche_plateau_2012.downloaddate = {'2013-03-08'};
peyroche_plateau_2012.pmid = 22586468;
peyroche_plateau_2012.desc = {'Relative fitness defect: rf = log2(wt(Se)/mut(Se) - wt/mut + 1), where wt and mut are the generation times of the WT and mutant strains with and without Na2Se.'};

phenotypes = {'Growth, log2 ratio'};
treatments = {'Sodium selenide, 1 uM 16 h'; 'Sodium selenide, 2 uM 16 h'; 'Sodium selenide, 1 uM 27 h'; 'Sodium selenide, 2 uM 27 h'};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/journal.pone.0036343.s004.xlsx', 'data');

% Get indices of the data columns
ind_data = 5:8;

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,1)),strmatch('Y', data.raw(:,1)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,1) = cellfun(@strtrim, data.raw(:,1),'UniformOutput',0);

data2.orfs = upper(data.raw(:,1));
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
peyroche_plateau_2012.orfs = t;
peyroche_plateau_2012.data = t2;
peyroche_plateau_2012.ph = strcat(phenotypes, {'; '}, treatments);

save('./peyroche_plateau_2012.mat','peyroche_plateau_2012');

fid = fopen('./peyroche_plateau_2012.txt','w');
write_matrix_file(fid, peyroche_plateau_2012.orfs, peyroche_plateau_2012.ph, peyroche_plateau_2012.data);
fclose(fid);

end
