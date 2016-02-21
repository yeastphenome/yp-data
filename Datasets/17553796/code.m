%% Thevissen~Francois, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

thevissen_francois_2007.pmid = 17553796;

phenotypes = {'growth [MIC]'};
treatments = {'miconazole [0.025-12.5 ug/ml]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Euroscarf library.xlsx', 'Tabelle1');
tested_orfs = tested.raw(2:end,2);


inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
fid = fopen('raw_data/hits_orfs_scores.txt');
hits = textscan(fid,'%s %d');
hits_orfs = upper(hits{1});
hits_scores = -hits{2};
fclose(fid);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustment
tested_orfs = [tested_orfs; missing];   % 2 ORFs added

thevissen_francois_2007.orfs = tested_orfs;
thevissen_francois_2007.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
thevissen_francois_2007.data(ind2,1) = hits_scores(ind1);

thevissen_francois_2007.ph = [strcat(phenotypes, '; ', treatments)];

save('./thevissen_francois_2007.mat','thevissen_francois_2007');

fid = fopen('./thevissen_francois_2007.txt','w');
write_matrix_file(fid, thevissen_francois_2007.orfs, thevissen_francois_2007.ph, thevissen_francois_2007.data);
fclose(fid);

end
