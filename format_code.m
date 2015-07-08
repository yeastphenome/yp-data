function format_code(parse_dir, ext)
% format all files of format EXT
% Input:
%   parse_dir : root dir of code to be formated, must be full path
%   ext       : code file type, i.e. '*.m', '*.c'

pat = '\*\.\w{1,4}';
mat_res = regexp(ext, pat, 'match');
assert(~isempty(mat_res), sprintf('File extension ''%s'' not right.', ext));

dir_c = regexp(parse_dir, ';', 'split')';
cellfun(@(x) format_dir(x, ext), dir_c);

disp('done.');

function format_dir(dir_c, ext)
if isempty(dir_c), return; end

disp(['[+]' dir_c]);
f_list = dir(fullfile(dir_c, ext));
f_list = {f_list.name}';
f_list = cellfun(@(x) fullfile(dir_c, x), f_list, 'UniformOutput', false);
cellfun(@format_file, f_list);

function format_file(f_path)
disp(f_path);
doc = matlab.desktop.editor.openDocument(f_path);
assert(~isempty(doc), sprintf('! %s :open failed.', f_path));
doc.smartIndentContents();
doc.save();
doc.close();
