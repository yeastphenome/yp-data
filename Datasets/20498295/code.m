%% Dakshinamurthy~Garfinkel, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
dakshinamurthy_garfinkel_2010.pmid = 20498295;

phenotypes = {'number of cells (His+) to monitor the frequency of Ty1 transposition'};
treatments = {''};

%% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Matalphakos counted.xlsx', 'Sheet2');
tested_orfs = tested.raw(6:end,2);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = clean_orf(tested_orfs);

% Fix typos
tested_orfs(find(strcmp('YLR228', tested_orfs))) = {'YLR228C'};
tested_orfs(find(strcmp('YMR062', tested_orfs))) = {'YMR062C'};
tested_orfs(find(strcmp('YYKL138C', tested_orfs))) = {'YKL138C'};

inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));

tested_orfs(inds) = [];

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/TableS1-2.xlsx', 'Sheet2');
hits_orfs = data.raw(:,2);
hits_data = data.raw(:,3);

inds = find(cellfun(@isempty, hits_orfs) | cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

hits_orfs = clean_orf(hits_orfs);

inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));

hits_orfs(inds) = [];
hits_data(inds) = [];

hits_data = cellfun(@length, hits_data);
hits_data = -hits_data; % these data describes a decrease in the number if His+ cells, thus a decrease in the frequency of Ty1 transposition

missing = setdiff(hits_orfs, tested_orfs);

% Adding 4 ORFs to the list of tested
dakshinamurthy_garfinkel_2010.orfs = [tested_orfs; missing];

[~,ind1,ind2] = intersect(dakshinamurthy_garfinkel_2010.orfs, hits_orfs);
dakshinamurthy_garfinkel_2010.data(ind1,1) = hits_data(ind2);

dakshinamurthy_garfinkel_2010.ph = strcat(phenotypes, '; ', treatments);

save('./dakshinamurthy_garfinkel_2010.mat','dakshinamurthy_garfinkel_2010');

fid = fopen('./dakshinamurthy_garfinkel_2010.txt','w');
write_matrix_file(fid, dakshinamurthy_garfinkel_2010.orfs, dakshinamurthy_garfinkel_2010.ph, dakshinamurthy_garfinkel_2010.data);
fclose(fid);

end
