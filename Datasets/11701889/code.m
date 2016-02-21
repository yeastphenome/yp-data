%% Ooi~Boeke, 2001
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

ooi_boeke_2001.source = {'Main PDF'};
ooi_boeke_2001.downloaddate = {'2014-02-03'};
ooi_boeke_2001.pmid = 11701889;

[FILENAMES{end+1}, hits] = read_data('textread','./raw_data/ooi_boeke_2001.txt', '%s');
hits = lower(hits);
hits(strcmp('lig4', hits)) = {'dnl4'};
hits(strcmp('gpe2', hits)) = {'yal056w'};

phenotypes = {'NHEJ defect'};
treatments = {''};

% Transform gene names into ORFs
hits_orf = translate(upper(hits));

ooi_boeke_2001.orfs = hits_orf;
ooi_boeke_2001.data = -ones(size(hits_orf));
ooi_boeke_2001.ph = [strcat(phenotypes, '; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(ooi_boeke_2001.orfs, essential_genes);
ooi_boeke_2001.orfs(ind1) = [];
ooi_boeke_2001.data(ind1,:) = [];

save('./ooi_boeke_2001.mat','ooi_boeke_2001');

fid = fopen('./ooi_boeke_2001.txt','w');
write_matrix_file(fid, ooi_boeke_2001.orfs, ooi_boeke_2001.ph, ooi_boeke_2001.data);
fclose(fid);

end
