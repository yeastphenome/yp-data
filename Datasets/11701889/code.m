%% Ooi~Boeke, 2001
function FILENAMES = code()
FILENAMES = {};
% NEED = tested genes

ooi_boeke_2001.source = {'Main PDF'};
ooi_boeke_2001.downloaddate = {'2014-02-03'};
ooi_boeke_2001.pmid = 11701889;

[FILENAMES{end+1}, hits] = dataread('textread','./raw_data/ooi_boeke_2001.txt', '%s');
hits = lower(hits);
hits(strcmp('lig4', hits)) = {'dnl4'};
hits(strcmp('gpe2', hits)) = {'yal056w'};

phenotypes = {'NHEJ defect'};
treatments = {''};

% Transform gene names into ORFs
hits_orf = genename2orf(hits);
tmp = split_by_delimiter('_',hits_orf);
hits_orf = tmp(:,1);

ooi_boeke_2001.orfs = upper(hits_orf);
ooi_boeke_2001.data = -ones(size(hits_orf));
ooi_boeke_2001.ph = [strcat(phenotypes, '; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(ooi_boeke_2001.orfs, essential_genes);
ooi_boeke_2001.orfs(ind1) = [];
ooi_boeke_2001.data(ind1,:) = [];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'ooi_boeke_2001.mat'],'ooi_boeke_2001');
return;

% % Save data into database
% dt = ooi_boeke_2001;
% datasets = get_datasets_for_paper(dt);
%
% [~,database_ix] = sortrows(datasets.names,[1 2 3]);
% [~,ph_ix] = sort(dt.ph);
%
% % Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
% datasets.names(database_ix,:)
% dt.ph(ph_ix)
%
% insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

