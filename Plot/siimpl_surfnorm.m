
function A = implantdir(varargin);

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

tag = 'Surface Normal:';

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

         fclose(fid);
         return
      end

end % while
   
disp('Import direction failed!!!')  
fclose(fid);

