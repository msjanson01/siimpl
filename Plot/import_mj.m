% function [varargout]=import(filename);

function varargout=import(filename);

fid=fopen(filename);

if fid==-1
   disp(['Could not open file:' ,filename]);
   for i=1:nargout, varargout{i}=[]; end
   return    
end
fclose(fid);


A = load(filename);

if isempty(A)
   disp('Could not read file into matrix!');
   for i=1:nargout, varargout{i}=[]; end
   return    
end

if (nargout==1)
   varargout{1}=A;
	return   
end

[varLen, Nvar]=size(A);

for i=1:min(nargout,Nvar)
   varargout{i}=A(:,i)';
end

if (nargout>Nvar);
   for i=(Nvar+1):nargout;
      varargout{i}=[];
   end
end