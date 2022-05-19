
function [Z1, M1, E0] = implantdir(varargin);

if nargin == 1
   filename = varargin{1};
else
   filename = 'siimplLog.txt';
end

fid=fopen(filename);

if fid==-1
   disp(['Could not open file:' ,filename]);
   for i=1:nargout, varargout{i}=[]; end
   return    
end

tag = 'Z1';

i = 0;
tries = 1000;

while (i < tries)
   
      i = i + 1;
      string=fgetl(fid);
      if ~isstr(string), break, end
      
      % find tag 
      if findstr(string, tag)
         
         [trsh, string] = strtok(string);
         Z1 = str2num(strtok(string));

         string=fgetl(fid);
         [trsh, string] = strtok(string);
         M1 = str2num(strtok(string));

         string=fgetl(fid);
         [trsh, string] = strtok(string);
         E0 = str2num(strtok(string));

         string=fgetl(fid);
         [trsh, string] = strtok(string);
         dose = str2num(strtok(string));

        fclose(fid);
         return
      end

end % while
   
disp('Import direction failed!!!')  
fclose(fid);

