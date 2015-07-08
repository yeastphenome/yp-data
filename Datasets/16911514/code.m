%% Kawahata~Iefuji, 2006
% DATA = kawahata_iefuji_2006
function FILENAMES = code()
FILENAMES = {};
kawahata_iefuji_2006.pmid = 16911514;

phenotypes = {'growth [spot assay]'};
treatments = {'lactic acid [5.1% w/v], pH [2.7]'; 'lactic acid [3.1% w/v], pH [2.9]'; 'HCl [0.28% w/v], pH [2.4]';'HCl [0.24% w/v], pH [2.6]'; 'acetic acid [0.5% w/v], pH [4.2]'; 'acetic acid [0.4% w/v], pH [4.3]'};

% Load resistant
[FILENAMES{end+1}, hits_resistant.raw] = dataread('xlsread','raw_data/hits.xlsx', 'Resistant');

hits_resistant_orfs = [hits_resistant.raw(:,1); hits_resistant.raw(:,6)];

inds = cellfun(@isnumeric, hits_resistant_orfs);
hits_resistant_orfs(inds) = [];
hits_resistant_orfs = cellfun(@strtrim, hits_resistant_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_resistant_orfs,1);
hits_resistant_orfs(inds) = [];
hits_resistant_orfs = unique(upper(hits_resistant_orfs));

hits_resistant_scores = zeros(length(hits_resistant_orfs),3);

% Scores
tmp_orfs = [hits_resistant.raw(strcmp('+', hits_resistant.raw(:,3)),1); hits_resistant.raw(strcmp('+', hits_resistant.raw(:,8)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_resistant_orfs);
hits_resistant_scores(ind2,1) = 1;

tmp_orfs = [hits_resistant.raw(strcmp('+', hits_resistant.raw(:,4)),1); hits_resistant.raw(strcmp('+', hits_resistant.raw(:,9)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_resistant_orfs);
hits_resistant_scores(ind2,2) = 1;

tmp_orfs = [hits_resistant.raw(strcmp('+', hits_resistant.raw(:,5)),1); hits_resistant.raw(strcmp('+', hits_resistant.raw(:,10)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_resistant_orfs);
hits_resistant_scores(ind2,3) = 1;

% Load sensitive
[FILENAMES{end+1}, hits_sensitive.raw] = dataread('xlsread','raw_data/hits.xlsx', 'Sensitive');

hits_sensitive_orfs = [hits_sensitive.raw(:,1); hits_sensitive.raw(:,6)];

inds = cellfun(@isnumeric, hits_sensitive_orfs);
hits_sensitive_orfs(inds) = [];
hits_sensitive_orfs = cellfun(@strtrim, hits_sensitive_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_sensitive_orfs,1);
hits_sensitive_orfs(inds) = [];
hits_sensitive_orfs = unique(upper(hits_sensitive_orfs));

hits_sensitive_scores = zeros(length(hits_sensitive_orfs),3);

% Scores
tmp_orfs = [hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,3)),1); hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,8)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_sensitive_orfs);
hits_sensitive_scores(ind2,1) = -1;

tmp_orfs = [hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,4)),1); hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,9)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_sensitive_orfs);
hits_sensitive_scores(ind2,2) = -1;

tmp_orfs = [hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,5)),1); hits_sensitive.raw(strcmp('-', hits_sensitive.raw(:,10)),6)];
[~,~,ind2] = intersect(upper(tmp_orfs), hits_sensitive_orfs);
hits_sensitive_scores(ind2,3) = -1;


% Check overlap between resistant and sensitive
length(intersect(hits_resistant_orfs(hits_resistant_scores(:,3)>0), hits_sensitive_orfs(hits_sensitive_scores(:,3)<0)));

kawahata_iefuji_2006.orfs = unique([hits_resistant_orfs; hits_sensitive_orfs]);
kawahata_iefuji_2006.data = zeros(length(kawahata_iefuji_2006.orfs), length(phenotypes));
[~,ind1,ind2] = intersect(hits_resistant_orfs, kawahata_iefuji_2006.orfs);
kawahata_iefuji_2006.data(ind2,[1 3 5]) = hits_resistant_scores(ind1,:);
[~,ind1,ind2] = intersect(hits_sensitive_orfs, kawahata_iefuji_2006.orfs);
kawahata_iefuji_2006.data(ind2,[2 4 6]) = hits_sensitive_scores(ind1,:);


kawahata_iefuji_2006.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'kawahata_iefuji_2006.mat'],'kawahata_iefuji_2006');
return;

% Save data into database
dt = kawahata_iefuji_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

