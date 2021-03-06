
function [Z2, M2, Natom] = implantdir(varargin);

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


tag = 'Z2';

tries = 1000;
i = 0;

while (i < tries)
   
      i = i + 1;
      string=fgetl(fid);
      if ~isstr(string), break, end
      
      % find tag 
      if findstr(string, tag)
         
         [trsh, string] = strtok(string);
         
         Z2 = str2num(string);
                  
         M2 = [];
         string=fgetl(fid);
         [trsh, string] = strtok(string);
         for j = 1:length(Z2)
            [M2i, string] = strtok(string);
            M2 = [M2 str2num(M2i)];
         end
        i = tries;

      end

end % while
   

tag = 'nAtom';

i = 0;

while (i < tries)
   
      i = i + 1;
      string=fgetl(fid);
      if ~isstr(string), break, end
      
      % find tag 
      if findstr(string, tag)
         
         [trsh, string] = strtok(string);
         
         Natom = str2num(strtok(string));

         fclose(fid);
         return
     end   
end % while

disp('Import direction failed!!!')  
fclose(fid);

