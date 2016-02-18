%% Giaever~Johnston, 2002
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
giaever_johnston_2002.pmid = 12140549;

%% First set of data: growth in various conditions

% Load the phenotype mapping file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/phenotype_mapping.xlsx','Sheet1');
pm.datasetid0 = data(2:end,1);
pm.who = data(2:end,7);
pm.what = data(2:end,8);
pm.when = data(2:end,6);
pm.how = data(2:end,9);
pm.condition = data(2:end,3);

pm.datasetid0 = cellfun(@num2str, pm.datasetid0,'UniformOutput',0);

for i = 1 : length(pm.condition)
    inds = regexp(pm.condition{i}, ' ');
    pm.condition{i} = pm.condition{i}(1:inds(end-2)-1);
end

pm = build_phenotype_name(pm);

[phenotypes, ia, ic] = unique(pm.ph);

hit_data_files = dir('./raw_data/');
hit_data_files_names = {hit_data_files(:).name}';
inds = find(cellfun(@isempty, regexp(hit_data_files_names,'^[0-9]{1,2}\_')));
hit_data_files_names(inds) = [];

phenotypes_hit_strains = cell(size(phenotypes));
phenotypes_hit_data = cell(size(phenotypes));

for i = 1 : length(phenotypes)
    
        % Find dataset_ids
        inds = find(ic == i);
        
        % Find the files
        e = strjoin(pm.datasetid0(inds),'|');
        indfiles = find(~cellfun(@isempty, regexp(hit_data_files_names, ['^(' e ')\_'])));

        all_hit_strains = [];
        all_hit_data = [];
        
        for j = 1 : length(indfiles)
            % Load hit strains
            [FILENAMES{end+1}, data] = read_data('textscan',['./raw_data/' hit_data_files_names{indfiles(j)}], '%s %f %*[^\n]');

            % Get the list of ORFs and the correponding data 
            hit_strains = data{1};

            % Get the data itself
            hit_data = exp(-data{2});   % this step is necessary to combine sensitivity and resistance data together (otherwise, they both have positive & negative values)
            if ~isempty(regexp(hit_data_files_names{indfiles(j)},'sen'))
                hit_data = -hit_data;
            end
    
            % Eliminate all white spaces & capitalize
            hit_strains = clean_orf(hit_strains);

            % Find anything that doesn't look like an ORF
            inds = find(~is_orf(hit_strains));
            disp(hit_strains(inds)); 
            
            all_hit_strains = [all_hit_strains; hit_strains];
            all_hit_data = [all_hit_data; hit_data];
        end

    % If the same strain is present more than once, average its values
    [all_hit_strains, all_hit_data] = grpstats(all_hit_data, all_hit_strains, {'gname','mean'});
    
    phenotypes_hit_strains{i} = all_hit_strains;
    phenotypes_hit_data{i} = all_hit_data;
end

% Combine all data together into a single matrix
all_orfs = unique(vertcat(phenotypes_hit_strains{:}));
all_data = nan(length(all_orfs),length(phenotypes));
for i = 1 : length(phenotypes)
    [~,ind1,ind2] = intersect(all_orfs, phenotypes_hit_strains{i});
    all_data(ind1,i) = phenotypes_hit_data{i}(ind2);
end

%% Second set of data: morphology

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

% Transform data
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/phenotype_mapping2.xlsx','Sheet1');
pm2.who = data(2:end,2);
pm2.what = data(2:end,3);
pm2.where = data(2:end,4);
pm2.when = data(2:end,5);
pm2.how = data(2:end,6);
pm2.coeff = cell2mat(data(2:end,7));
pm2.orig = data(2:end,1);

pm2 = build_phenotype_name(pm2);

morphology_phenotypes = unique(pm2.ph);

hit_data2 = zeros(length(hit_strains), length(morphology_phenotypes));

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
                ind2 = find(strcmp(pm2.ph{ind1}, morphology_phenotypes));
                hit_data2(i,ind2) = str2num(tmp2{2});
                if strcmp(pm2.what{ind1}, 'size')
                    hit_data2(i,ind2) = pm2.coeff(ind1) .* hit_data2(i,ind2);
                end
            end
        end        
    end
    
end

% If the same strain is present more than once, average its values
[all_hit_strains2, all_hit_data2] = grpstats(hit_data2, hit_strains, {'gname','mean'});

%% Merge the growth data and the morphology data

all_phenotypes3 = [phenotypes; morphology_phenotypes];
all_orfs3 = unique([all_orfs; all_hit_strains2]);
all_data3 = nan(length(all_orfs3), length(all_phenotypes3));

[~,ind1,ind2] = intersect(all_orfs, all_orfs3);
all_data3(ind2,1:length(phenotypes)) = all_data(ind1,:);

[~,ind1,ind2] = intersect(all_hit_strains2, all_orfs3);
all_data3(ind2,length(phenotypes)+1:length(phenotypes)+length(morphology_phenotypes)) = all_hit_data2(ind1,:);


%% Prepare final dataset

% If the dataset is quantitative:
giaever_johnston_2002.orfs = all_orfs3;
giaever_johnston_2002.ph = all_phenotypes3;
giaever_johnston_2002.data = all_data3;

%% Save

save('./giaever_johnston_2002.mat','giaever_johnston_2002');

end

