%% Giaever~Johnston, 2002
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
giaever_johnston_2002.pmid = 12140549;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(giaever_johnston_2002.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% First set of data: growth in various conditions

% Load the phenotype mapping file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/phenotype_mapping.xlsx','Sheet1');
pm.experiment = cell2mat(data(2:end,1));
pm.datasetid = cell2mat(data(2:end,2));

hit_data_files = dir('./raw_data/');
hit_data_files_names = {hit_data_files(:).name}';
inds = find(cellfun(@isempty, regexp(hit_data_files_names,'^[0-9]{1,2}\_')));
hit_data_files_names(inds) = [];

hit_dataset_ids = unique(pm.datasetid);

phenotypes_hit_strains = cell(length(hit_dataset_ids),1);
phenotypes_hit_data = cell(length(hit_dataset_ids),1);

for i = 1 : length(hit_data_files_names)
            
        % Get the dataset id
        t = regexp(hit_data_files_names{i}, '(\d)+(?=\_)','match');
        t = str2num(t{1});
        
        dt = pm.datasetid(pm.experiment==t);

        % Load hit strains
        [FILENAMES{end+1}, data] = read_data('textscan',['./raw_data/' hit_data_files_names{i}], '%s %f %*[^\n]');

        % Get the list of ORFs and the correponding data 
        hit_strains = data{1};

        % Get the data itself
        hit_data = data{2};
        if ~isempty(regexp(hit_data_files_names{i},'sen'))
            hit_data = -hit_data;
        end

        % Eliminate all white spaces & capitalize
        hit_strains = clean_orf(hit_strains);
        
        % If in gene name form, transform into ORF name
        hit_strains = translate(hit_strains);

        % Find anything that doesn't look like an ORF
        inds = find(~is_orf(hit_strains));
        disp(hit_strains(inds)); 
        
        phenotypes_hit_strains{hit_dataset_ids==dt} = [phenotypes_hit_strains{hit_dataset_ids==dt}; hit_strains];
        phenotypes_hit_data{hit_dataset_ids==dt} = [phenotypes_hit_data{hit_dataset_ids==dt}; hit_data];

end

for i = 1 : length(hit_dataset_ids)

    % If the same strain is present more than once, average its values
    [all_hit_strains, all_hit_data] = grpstats(phenotypes_hit_data{i}, phenotypes_hit_strains{i}, {'gname','mean'});
    
    phenotypes_hit_strains{i} = all_hit_strains;
    phenotypes_hit_data{i} = all_hit_data;
    
end

% Combine all data together into a single matrix
all_orfs = unique(vertcat(phenotypes_hit_strains{:}));
all_data = nan(length(all_orfs),length(hit_dataset_ids));
for i = 1 : length(hit_dataset_ids)
    [~,ind1,ind2] = intersect(all_orfs, phenotypes_hit_strains{i});
    all_data(ind1,i) = phenotypes_hit_data{i}(ind2);
end

%% Second set of data: YPD

[FILENAMES{end+1}, data_ypd] = read_data('readtable','./raw_data/ypd.txt','delimiter','\t');
all_orfs_ypd = data_ypd.ORF;
all_data_ypd = -data_ypd.AverageRatio;

all_orfs_ypd = clean_orf(all_orfs_ypd);

% If in gene name form, transform into ORF name
all_orfs_ypd = translate(all_orfs_ypd);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(all_orfs_ypd));
disp(all_orfs_ypd(inds)); 

% If the same strain is present more than once, average its values
[all_orfs_ypd, all_data_ypd] = grpstats(all_data_ypd, all_orfs_ypd, {'gname','mean'});

hit_dataset_ids_ypd = [16187];

all_orfs2 = unique([all_orfs; all_orfs_ypd]);
all_datasets2 = [hit_dataset_ids; hit_dataset_ids_ypd];
all_data2 = zeros(length(all_orfs2), length(all_datasets2));

[~,ind1,ind2] = intersect(all_orfs, all_orfs2);
all_data2(ind2,1:length(hit_dataset_ids)) = all_data(ind1,:);
[~,ind1,ind2] = intersect(all_orfs_ypd, all_orfs2);
all_data2(ind2,length(hit_dataset_ids)+1:end) = all_data_ypd(ind1,:);


%% Third set of data: morphology

[FILENAMES{end+1}, data] = read_data('textscan','./raw_data/Cell_Morph_Screen_Table.txt','%s %s %s %s','delimiter','\t');

hit_strains = data{1}(3:end);
hit_data = data{4}(3:end);

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds)); 

% Fix errors
hit_strains(strcmp('YELOO1C', hit_strains)) = {'YEL001C'};
inds = find(strcmp('YMR41W', hit_strains));    % ambiguos typo, can't be fixed
hit_strains(inds) = [];
hit_data(inds) = [];

hit_strains = translate(hit_strains);


% Transform data
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/phenotype_mapping2.xlsx','Sheet1');

pm2.datasetid = cell2mat(data(2:end,8));
pm2.coeff = cell2mat(data(2:end,7));
pm2.orig = data(2:end,1);

hit_dataset_ids2 = unique(pm2.datasetid);
hit_data2 = zeros(length(hit_strains), length(hit_dataset_ids2));

for i = 1 : length(hit_strains)
    
    tmp = regexp(hit_data{i}, ';', 'split');
    for j = 1 : length(tmp)
        tmp{j} = strtrim(tmp{j});
        tmp2 = regexp(tmp{j}, ' ', 'split');
        
        if strcmp(tmp2{1},'WT')
            hit_data2(i,:) = 0;
        elseif strcmp(tmp2{1},'ND')
            hit_data2(i,:) = NaN;
        else
            ind1 = find(strcmp(tmp2{1}, pm2.orig));
            
            if ~isempty(ind1)   % if isempty, the phenotype is some of the minor ones and we don't record it.
                dataset_id = pm2.datasetid(ind1);
                coeff = pm2.coeff(ind1);
                dataset_indx = find(hit_dataset_ids2 == dataset_id);
            
                hit_data2(i,dataset_indx) = coeff * str2num(tmp2{2});
            end
        end        
    end
    
end

% If the same strain is present more than once, average its values
[all_hit_strains2, all_hit_data2] = grpstats(hit_data2, hit_strains, {'gname','mean'});

%% Merge the growth data and the morphology data

all_dataset_ids3 = [all_datasets2; hit_dataset_ids2];
all_orfs3 = unique([all_orfs2; all_hit_strains2]);
all_data3 = nan(length(all_orfs3), length(all_dataset_ids3));

[~,ind1,ind2] = intersect(all_orfs2, all_orfs3);
all_data3(ind2,1:length(all_datasets2)) = all_data2(ind1,:);

[~,ind1,ind2] = intersect(all_hit_strains2, all_orfs3);
all_data3(ind2,length(all_datasets2)+1:end) = all_hit_data2(ind1,:);


%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, all_dataset_ids3);
hit_data_names = cell(size(all_dataset_ids3));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
giaever_johnston_2002.orfs = all_orfs3;
giaever_johnston_2002.ph = hit_data_names;
giaever_johnston_2002.data = all_data3;
giaever_johnston_2002.dataset_ids = all_dataset_ids3;

%% Save

save('./giaever_johnston_2002.mat','giaever_johnston_2002');

%% Print out

fid = fopen('./giaever_johnston_2002.txt','w');
write_matrix_file(fid, giaever_johnston_2002.orfs, giaever_johnston_2002.ph, giaever_johnston_2002.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(giaever_johnston_2002)
end

end

