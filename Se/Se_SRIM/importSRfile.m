function [E, Se, Sn] = importSRfile(filename);

fid = fopen(filename,'r');

tries = 100;

i = 1;
found1 = 0;

while ((i < tries) & (~found1))
      i = i + 1;
      
      string = fgetl(fid);
      if ~isstr(string), break, end
            
      % Find  data (after Line that start with 'Energy')
      if (findstr(string,'Energy') & ~found1)
         found1 = 1;
         data = fscanf(fid,'%f%f%f',[3 inf]);         
      end 
end % while
data = data';

E = data(:,1);
Se = data(:,2);
Sn = data(:,3);

