%% Dunn~Jensen, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};


dunn_jensen_2006.source = {'http://www.molbiolcell.org/content/suppl/2005/11/02/E05-06-0585.DC1/Supp_Table3.xls'};
dunn_jensen_2006.downloaddate = {'2013-03-06'};
dunn_jensen_2006.pmid = 16267274;

phenotypes = {'Growth, log2 hybridization ratio'};
treatments = {'EtBr, 25 ug/ml'};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Supp_Table3.xlsx', 'Sheet1');

ind_orf = strmatch('Systemic name', data.raw(12,:));
data2.orfs = upper(data.raw(14:end, ind_orf));

ind_data1 = strmatch('(-EtBr/+EtBr)', data.raw(12,:));

data2.data = data.raw(14:end, ind_data1);

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

inds = setdiff(1:length(data2.orfs),strmatch('Y', data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);


% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
dunn_jensen_2006.orfs = t;
dunn_jensen_2006.data = 1./t2;  % reverse the ratio so that the lower the value, the sicker the mutant.
dunn_jensen_2006.ph = strcat(phenotypes, {'; '}, treatments);

save('./dunn_jensen_2006.mat','dunn_jensen_2006');
return;

% Save data into database
dt = dunn_jensen_2006;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
datasets_ids(i,1) = datasets(i).id;
datasets_names{i,1} = datasets(i).name;
datasets_names{i,2} = datasets(i).shortname;
datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

fid = fopen('./dunn_jensen_2006.txt','w');
write_matrix_file(fid, dunn_jensen_2006.orfs, dunn_jensen_2006.ph, dunn_jensen_2006.data);
fclose(fid);

end

