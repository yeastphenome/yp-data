%% Lee~Giaever, 2014
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
lee_giaever_2014.pmid = 24723613;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(lee_giaever_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load text file
[FILENAMES{end+1}, data] = read_data('read_matrix_file', './raw_data/fitness_defect_matrix_hom.txt',1,1);
data.labels_col = strtrim(data.labels_col);

[FILENAMES{end+1}, data_het] = read_data('read_matrix_file', './raw_data/fitness_defect_matrix_het.txt',1,1);
data_het.labels_col = strtrim(data_het.labels_col);

% Get compound names
[FILENAMES{end+1}, compounds] = read_data('readtable','./raw_data/1250217s1.xlsx', 'Sheet','compound library');
mu = native2unicode([206 188],'UTF-8');

inds = find(strcmp([mu 'M'], compounds.unit));
compounds.unit(inds) = {'uM'};

inds = find(strcmp('-', compounds.PCID));
compounds.PCID(inds) = {'NaN'};

compounds.PCID = cellfun(@str2num, compounds.PCID);

% Match screenIDs with dataset ids
[FILENAMES{end+1}, screen2dataset] = read_data('readtable', './extras/screenid_datasetid.txt','delimiter','\t','ReadVariableNames', false);
screen2dataset.Properties.VariableNames = {'screen','hom_dataset_id'};

[FILENAMES{end+1}, hom_dataset2conditionset] = read_data('readtable', './private/hom_datasetid_conditionsetid.txt','delimiter','\t','ReadVariableNames', false);
hom_dataset2conditionset.Properties.VariableNames = {'hom_dataset_id','conditionset_id'};

[FILENAMES{end+1}, het_dataset2conditionset] = read_data('readtable', './private/het_datasetid_conditionsetid.txt','delimiter','\t','ReadVariableNames', false);
het_dataset2conditionset.Properties.VariableNames = {'het_dataset_id','conditionset_id'};

% Keep only one dataset id whenever the same conditionset occurs

% Hom
[C, ia, ic] = unique(hom_dataset2conditionset.conditionset_id);

counts = zeros(length(C),1);
for i = 1 : length(C)
    counts(i) = length(find(hom_dataset2conditionset.conditionset_id == C(i)));
end
dups = C(counts>1);
sortrows(hom_dataset2conditionset(ismember(hom_dataset2conditionset.conditionset_id, dups),:), 'conditionset_id')

datasets_to_remove = [];
for i = 1 :length(dups)
    datasets_to_remove = [datasets_to_remove; 
        min(hom_dataset2conditionset{hom_dataset2conditionset.conditionset_id == dups(i),'hom_dataset_id'})];
end

hom_dataset2conditionset(ismember(hom_dataset2conditionset.hom_dataset_id, datasets_to_remove),:) = [];

% Het
[C, ia, ic] = unique(het_dataset2conditionset.conditionset_id);

counts = zeros(length(C),1);
for i = 1 : length(C)
    counts(i) = length(find(het_dataset2conditionset.conditionset_id == C(i)));
end
dups = C(counts>1);
sortrows(het_dataset2conditionset(ismember(het_dataset2conditionset.conditionset_id, dups),:), 'conditionset_id')

datasets_to_remove = [];
for i = 1 :length(dups)
    datasets_to_remove = [datasets_to_remove; 
        min(het_dataset2conditionset{het_dataset2conditionset.conditionset_id == dups(i),'het_dataset_id'})];
end

het_dataset2conditionset(ismember(het_dataset2conditionset.het_dataset_id, datasets_to_remove),:) = [];

screen2dataset = innerjoin(screen2dataset, hom_dataset2conditionset, 'Keys', {'hom_dataset_id'});
screen2dataset = innerjoin(screen2dataset, het_dataset2conditionset, 'Keys', {'conditionset_id'});


%% HOM 
% Match data with dataset ids
[~,ind1,ind2] = intersect(screen2dataset.screen, data.labels_col);
data.datasetid = nan(size(data.labels_col));
data.datasetid(ind2) = screen2dataset.hom_dataset_id(ind1);

% Clean up the ORFs
data.labels_row = clean_orf(data.labels_row);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(data.labels_row));
disp(data.labels_row(inds));

% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
[data.datasetid, data.data] = grpstats(data.data', data.datasetid, {'gname','mean'});
data.data = data.data';
data.datasetid = cellfun(@str2num, data.datasetid);

%% HET 
% Match data with dataset ids
[~,ind1,ind2] = intersect(screen2dataset.screen, data_het.labels_col);
data_het.datasetid = nan(size(data_het.labels_col));
data_het.datasetid(ind2) = screen2dataset.het_dataset_id(ind1);

% Clean up the ORFs
data_het.labels_row = clean_orf(data_het.labels_row);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(data_het.labels_row));
disp(data_het.labels_row(inds));

% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
[data_het.datasetid, data_het.data] = grpstats(data_het.data', data_het.datasetid, {'gname','mean'});
data_het.data = data_het.data';
data_het.datasetid = cellfun(@str2num, data_het.datasetid);

%% Prepare final dataset

strains_all = unique([data.labels_row; data_het.labels_row]);
datasetids_all = [data.datasetid; data_het.datasetid];
data_all = nan(length(strains_all), length(datasetids_all));

[~,ind1,ind2] = intersect(strains_all, data.labels_row);
data_all(ind1,1:length(data.datasetid)) = data.data(ind2,:);

[~,ind1,ind2] = intersect(strains_all, data_het.labels_row);
data_all(ind1,length(data.datasetid)+1:end) = data_het.data(ind2,:);

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, datasetids_all);
datanames_all = cell(size(datasetids_all));
datanames_all(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
lee_giaever_2014.orfs = strains_all;
lee_giaever_2014.ph = datanames_all;
lee_giaever_2014.data = -data_all;         % taking the opposite because the data are fitness scores = log2(ctrl/sample)
lee_giaever_2014.dataset_ids = datasetids_all;


%% Save

save('./lee_giaever_2014.mat','lee_giaever_2014');

fid = fopen('./lee_giaever_2014.txt','w');
write_matrix_file(fid, lee_giaever_2014.orfs, lee_giaever_2014.ph, lee_giaever_2014.data);
fclose(fid);

% fid = fopen('./private/hom_data_to_import.txt','w');
% t = data.datasetid;
% for i = 1 : length(t)
%     inds = find(lee_giaever_2014.dataset_ids == t(i));
%     for j = 1 : length(lee_giaever_2014.orfs)
%         if ~isnan(lee_giaever_2014.data(j, inds))
%             fprintf(fid, '%d\t%s\t%.3f\n', t(i), lee_giaever_2014.orfs{j}, lee_giaever_2014.data(j,inds));
%         end
%     end
%     if mod(i, 100) == 0
%         fprintf('%d of %d\n', i, length(t));
%     end
% end
% fclose(fid);
            
end