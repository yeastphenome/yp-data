%% Jonikas, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
jonikas_schuldiner_2009.pmid = 19325107;

treatments = {'standard'};
phenotypes = {'UPR response'};

% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/NIHMS201195-supplement-st1.xlsx','Table_S1');

orfs = data.raw(3:end,1);
raw_data = data.raw(3:end,11:end);
headings = data.raw(1:2, 11:19);

orfs = clean_orf(orfs);

orfs(strcmp('YOLO62C', orfs)) = {'YOL062C'};
orfs(strcmp('YKLO72W', orfs)) = {'YKL072W'};
orfs(strcmp('YOLO57W', orfs)) = {'YOL057W'};
orfs(strcmp('YLR287-A', orfs)) = {'YLR287C-A'};

inds = find(~is_orf(orfs));
disp(orfs(inds));

orfs(inds) = [];
raw_data(inds,:) = [];

raw_data = raw_data(:,7);
headings = headings(:,7);

raw_data = cell2mat(raw_data);

[orfs, raw_data] = grpstats(raw_data, orfs, {'gname','mean'});

jonikas_schuldiner_2009.orfs = orfs;
jonikas_schuldiner_2009.ph = strcat(phenotypes, '; ', treatments);
jonikas_schuldiner_2009.data = raw_data;

%% Save

save('./jonikas_schuldiner_2009.mat','jonikas_schuldiner_2009');

%% Print out

fid = fopen('./jonikas_schuldiner_2009.txt','w');
write_matrix_file(fid, jonikas_schuldiner_2009.orfs, jonikas_schuldiner_2009.ph, jonikas_schuldiner_2009.data);
fclose(fid);



end