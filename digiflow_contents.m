function digiflow_contents(category)

m_path = mpath;
string_chars = ' \n ';
if nargin == 0
    contents = ls([m_path]);
    category = '';
else
    contents = ls([m_path category]);
end
fprintf(['\n ------------------------ \n ' category ' \n ------------------------ \n'])
for i = 3:size(contents, 1)
    part = contents(i, :);
    [~, part] = fileparts(part);
    string_chars = [string_chars part];
    if mod(i-1, 4) == 1
        string_chars = [string_chars ' \n '];
    else
        string_chars = [string_chars ' \t '];
    end
end
string_chars = [string_chars ' \n '];

fprintf(string_chars);


end
