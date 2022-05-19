% function plotatomR3(R, id);

function plotatomsR3(R,id);

[m, n]=size(R);
if prod(size(id)) == 1, id = ones(m,1)*id; end


% primary ion
i = find(id == 0);
if ~isempty(i)
   lh = line(R(i,1),R(i,2),R(i,3));
   if length(R) > 100
       marker='.';
       color='k';
       siz = 3;
   else
    marker='o';
    color='g';
    siz=3;
end
   set(lh,'lineStyle','none');
   set(lh,'marker',marker,'markerSize',siz);
   set(lh,'markerFaceColor',color,'markerEdgeColor',color,'tag','ion');
end

% primary ion, dots
i = find(id == -100);
if ~isempty(i)
	lh = line(R(i,1),R(i,2),R(i,3));
   set(lh,'lineStyle','none');
   set(lh,'marker',marker,'markerSize',siz);
   set(lh,'markerFaceColor',color,'markerEdgeColor',color,'tag','ion');
end


% interstitlal 1
i = find(id == 1);
if ~isempty(i)
	lh = line(R(i,1),R(i,2),R(i,3));
   marker='o';
   color='r';
   siz=2;
   set(lh,'lineStyle','none');
   set(lh,'marker',marker,'markerSize',siz);
   set(lh,'markerFaceColor',color,'markerEdgeColor',color,'Tag','Int1');
end
% interstitlal 2
i = find(id == 2);
if ~isempty(i)
	lh = line(R(i,1),R(i,2),R(i,3));
   marker='o';
   color='b';
   siz=2;
   set(lh,'lineStyle','none');
   set(lh,'marker',marker,'markerSize',siz);
   set(lh,'markerFaceColor',color,'markerEdgeColor',color,'Tag','Int2');
end
% interstitlal 3,4,...
i = find(id >= 3);
if ~isempty(i)
	lh = line(R(i,1),R(i,2),R(i,3));
   marker='o';
   color='c';
   siz=2;
   set(lh,'lineStyle','none');
   set(lh,'marker',marker,'markerSize',siz);
   set(lh,'markerFaceColor',color,'markerEdgeColor',color,'Tag','Int3-');
end
% ------------------------

% vacancy 1
i = find(id == -1);
if ~isempty(i)
	lh = line(R(i,1),R(i,2),R(i,3));
   marker='o';
   color='k';
   siz=2;
   set(lh,'lineStyle','none');
   set(lh,'marker',marker,'markerSize',siz);
   set(lh,'markerFaceColor',color,'markerEdgeColor',color,'Tag','Vac1');
end
% vacancy 2
i = find(id == -2);
if ~isempty(i)
	lh = line(R(i,1),R(i,2),R(i,3));
   marker='o';
   color=[0.6 0.6 0.6];
   siz=2;
   set(lh,'lineStyle','none');
   set(lh,'marker',marker,'markerSize',siz);
   set(lh,'markerFaceColor',color,'markerEdgeColor',color,'Tag','Vac2');
end
% vacancy 3,4,...
i = find(id <= -3);
if ~isempty(i)
	lh = line(R(i,1),R(i,2),R(i,3));
   marker='o';
   color='y';
   siz=2;
   set(lh,'lineStyle','none');
   set(lh,'marker',marker,'markerSize',siz);
   set(lh,'markerFaceColor',color,'markerEdgeColor',color,'Tag','Vac3-');
end
