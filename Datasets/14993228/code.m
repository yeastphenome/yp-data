%% Serrano~Arino, 2004
% DATA = serrano_arino_2004
function FILENAMES = code()
FILENAMES = {};

serrano_arino_2004.source = {'main PDF'};
serrano_arino_2004.downloaddate = {'2014-03-03'};
serrano_arino_2004.pmid = 14993228;

[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/serrano_arino_2004.xlsx', 'Sheet1');

phenotypes = {'growth (colony size)'};
treatments = {'pH 6-2.-7.5'};

% Eliminate white spaces before/after ORF
data.raw(:,1) = cellfun(@strtrim, data.raw(:,1),'UniformOutput',0);

% Translate genenames to ORF
data.raw(:,4) = genename2orf(data.raw(:,1),'noannot');

% Manually adjust the genenames that couldn't be matched
data.raw(strcmpi('cwh36', data.raw(:,4)),4) = {'YCL007C'};
data.raw(strcmpi('lys7', data.raw(:,4)),4) = {'YMR038C'};
data.raw(strcmpi('rcs1', data.raw(:,4)),4) = {'YGL071W'};

hits_orfs = data.raw(:,4);
scores = cell2mat(data.raw(:,2));

% Flip the scores such that -5 is the most sensitive and -1 is the least
% sensitive
scores = scores - 6;

% Load tested genes
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','raw_data/BY4741.xlsx', 'Tabelle1');
tested_orfs = unique(upper(tested.raw(2:end,2)));

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
serrano_arino_2004.orfs = tested_orfs;
serrano_arino_2004.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(serrano_arino_2004.orfs, hits_orfs);
serrano_arino_2004.data(ind1,:) = scores(ind2,:);

serrano_arino_2004.ph = [strcat(phenotypes{1}, '; ', treatments)];


a = mfilename('fullpath');
a = a(1:end-4);
save([a,'serrano_arino_2004.mat'],'serrano_arino_2004');
return;

% Save data into database
dt = serrano_arino_2004;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

