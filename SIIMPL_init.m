SIIMPL_path = cd;

add_path = [SIIMPL_path, '\plot'];
if isempty(findstr(path, add_path))
   path(path, add_path);
end
