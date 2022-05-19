
function [Efermi, A1, A2, p2, A3, A4, A5] = implantdir(varargin);

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

tag = 'Efermi';

i = 0;
tries = 1000;

while (i < tries)
   
      i = i + 1;
      string=fgetl(fid);
      if ~isstr(string), break, end
      
      % find tag 
      if findstr(string, tag)
         
         [trsh, string] = strtok(string);
         Efermi = str2num(strtok(string));

         string=fgetl(fid);

         string=fgetl(fid);
         [trsh, string] = strtok(string);
         A1 = str2num(strtok(string));

         string=fgetl(fid);
         [trsh, string] = strtok(string);
         A2 = str2num(strtok(string));

         string=fgetl(fid);
         [trsh, string] = strtok(string);
         p2 = str2num(strtok(string));

         string=fgetl(fid);
         [trsh, string] = strtok(string);
         A3 = str2num(strtok(string));

         string=fgetl(fid);
         [trsh, string] = strtok(string);
         A4 = str2num(strtok(string));

         string=fgetl(fid);
         [trsh, string] = strtok(string);
         A5 = str2num(strtok(string));

        fclose(fid);
         return
      end

end % while
   
disp('Import direction failed!!!')  
fclose(fid);

