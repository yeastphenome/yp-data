%% Szymanski~Goodman, 2007
function FILENAMES = code()
FILENAMES = {};
szymanski_goodman_2007.pmid = 18093937;

treatments = {'log phase'; 'stationary phase'};


%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table2.xlsx', 'Sheet1');
data(1,:) = []; % remove column headers

% Get the list of ORFs and the correponding data 
hit_strains = data(:,1);

[FILENAMES{end+1}, data2] = read_data('xlsread', './raw_data/phenotype_mapping.xlsx','Sheet1');
ph_mapping.ph_orig = data2(2:end,1);
ph_mapping.mat = cell2mat(data2(2:end,2:end));
ph_mapping.ph_new = data2(1,2:end)';

for k = 1 : length(treatments)
    hit_data{k} = zeros(length(hit_strains),length(ph_mapping.ph_new));
    for i = 1 : length(hit_strains)
        ind = find(strcmp(hit_strains{i}, data(:,1)));
        phs = regexp(data{ind,k+1},',','split');
        phs = strtrim(lower(phs));
        for j = 1 : length(phs)
            ind2 = find(strcmp(phs{j}, ph_mapping.ph_orig));
            ind3 = find(abs(ph_mapping.mat(ind2,:))>0);
            hit_data{k}(ind,ind3) = ph_mapping.mat(ind2,ind3);
        end
    end
end

 % Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

%% Prepare final dataset

phenotypes = strcat({'lipid droplets '}, ph_mapping.ph_new);
conditions = [repmat({'log phase'}, length(phenotypes),1); repmat({'stationary phase'}, length(phenotypes),1)];
phenotypes = [phenotypes; phenotypes];

% If the dataset is quantitative:
szymanski_goodman_2007.orfs = hit_strains;
szymanski_goodman_2007.ph = strcat(phenotypes, '; ', conditions);
szymanski_goodman_2007.data = [hit_data{1} hit_data{2}];


%% Save

save('./szymanski_goodman_2007.mat','szymanski_goodman_2007');

end

