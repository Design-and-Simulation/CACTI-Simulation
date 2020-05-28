function origin2center(ax)
%ORIGIN2CENTER  put the axis's origin to the center of the figure
% 
%   Input: axes's handle


% extend the x,y axis lengths
xL=xlim;
yL=ylim;
extend_x = ( xL(2)-xL(1) ) * 0.1 ;
extend_y = ( yL(2)-yL(1) ) * 0.1 ;
xxL = xL + [ -extend_x extend_x] ;
yyL = yL + [ -extend_y extend_y] ;
set(gca,'xlim', xxL) ;
set(gca,'ylim', yyL) ;

% put origin to center
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.Box = 'off';

pos = get(gca,'Position');
x_Lim = get(gca,'Xlim');
y_Lim = get(gca,'Ylim');
if prod(y_Lim)>0
   position_x = [pos(1), pos(2)+pos(4)/2, pos(3), eps];
else
   position_x = [pos(1), pos(2)-y_Lim(1)/diff(y_Lim)*pos(4), pos(3), eps];
end
if prod(x_Lim)>0
   position_y = [pos(1)+pos(3)/2, pos(2), eps, pos(4)];
else
   position_y = [pos(1)-x_Lim(1)/diff(x_Lim)*pos(3), pos(2), eps, pos(4)];
end
annotation('arrow', [pos(1)-0.065*pos(3), pos(1)+pos(3)+0.065*pos(3)], ...
[position_x(2)-0.001,position_x(2)-0.001],'HeadLength',6,'HeadWidth',6);
annotation('arrow', [position_y(1)+0.001, position_y(1)+0.001],...
[pos(2)-0.065*pos(4),pos(2)+pos(4)+0.065*pos(4)],...
'HeadLength',6,'HeadWidth',6);
end