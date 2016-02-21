%% Corbacho~Hernandez, 2005
function FILENAMES = code()
FILENAMES = {};
corbacho_hernandez_2005.pmid = 15993632;

phenotypes = {'alcian blue staining intensity'};
treatments = {''};


%% Hit strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('textread','./raw_data/Corbacho_tables.txt','%s %s %s %s %*[^\n]','delimiter','\t');

% Get the list of ORFs and the correponding data 
hit_strains = data{1};

% Get the data itself
hit_data = data{4};
    
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

% One phenotype observed only in Mat-alpha strain. Remove the annotation
inds = find(strcmp('1 (alpha)', hit_data));
hit_data(inds) = {'1'};
hit_data = cellfun(@str2num, hit_data);

% Transform data from: 1 = large decrease of staining, 3 = small decrease of staining
% to: -3 = super-low staining, -1 = low staining
hit_data = hit_data-4;

%% Prepare final dataset

corbacho_hernandez_2005.orfs = hit_strains;
corbacho_hernandez_2005.ph = strcat(phenotypes, '; ', treatments);
corbacho_hernandez_2005.data = hit_data;

%% Save

save('./corbacho_hernandez_2005.mat','corbacho_hernandez_2005');

fid = fopen('./corbacho_hernandez_2005.txt','w');
write_matrix_file(fid, corbacho_hernandez_2005.orfs, corbacho_hernandez_2005.ph, corbacho_hernandez_2005.data);
fclose(fid);

end

