%% Parsons~Boone, 2004
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
% TESTED = not available

parsons_boone_2004.source = {'http://www.nature.com/nbt/journal/v22/n1/extref/nbt919-S2.xls'};
parsons_boone_2004.downloaddate = {'2014-02-20'};
parsons_boone_2004.pmid = 14661025;

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/nbt919-S2.xlsx', 'Sheet1');

inds = find(strcmp('ORF', data.raw(:,1)));

phenotypes = {'growth (colony size)'};
treatments = data.raw(inds,3:end)';

% Eliminate white spaces before/after ORF
data.raw(inds+1:end,1) = cellfun(@strtrim, data.raw(inds+1:end,1),'UniformOutput',0);

% Replace 'Inf' with Inf
t = data.raw(inds+1:end,3:end);
t(cellfun(@isnan, t)) = {0};

parsons_boone_2004.orfs = upper(data.raw(inds+1:end,1));
parsons_boone_2004.data = cell2mat(t);

% Flip the sign of the values, such that negative = sensitive, positive =
% resistant
parsons_boone_2004.data = -parsons_boone_2004.data;

parsons_boone_2004.ph = [strcat(phenotypes{1}, '; ', treatments)];


% Eliminate the essential genes
% load essential_genes_100908;
% [t,ind1,ind2] = intersect(parsons_boone_2004.orfs, essential_genes);
% parsons_boone_2004.orfs(ind1) = [];
% parsons_boone_2004.data(ind1,:) = [];

save('./parsons_boone_2004.mat','parsons_boone_2004');

fid = fopen('./parsons_boone_2004.txt','w');
write_matrix_file(fid, parsons_boone_2004.orfs, parsons_boone_2004.ph, parsons_boone_2004.data);
fclose(fid);

end

