% function [varargout]=import(filename);

function A = import(varargin);

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

tag = 'Implanted rectangle area (Å):';

i = 0;
tries = 1000;

while (i < tries)
   
      i = i + 1;
      string=fgetl(fid);
      if ~isstr(string), break, end
      
      % find tag 
      if findstr(string, tag)
         
    		string=fgetl(fid);
         A = str2num(string);
    		string=fgetl(fid);
         A = [A; str2num(string)];
    		string=fgetl(fid);
         A = [A; str2num(string)];
			string=fgetl(fid);
         A = [A; str2num(string)];
         A = [A; A(1,:)];
         fclose(fid);
         return
      end

end % while
   
disp('Import rectangle failed!!!')  
fclose(fid);

