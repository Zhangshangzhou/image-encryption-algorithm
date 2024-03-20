 function [x,y]=henon(x0,y0,num)

x=zeros(1);y=zeros(1);
x(1)=x0;
y(1)=y0;
a=0.6;
k=0.8;
for n=1:num
	x(n+1)=sin(21./a*(y(n)+3)*k*x(n)*(1-k*x(n)));
	y(n+1)=sin(21./(a*(k*x(n+1)+3)*y(n)*(1-y(n))));
end
% figure;
% H=plot(x(1000:end),y(1000:end),'b');%1000
% set(H,'linestyle','none','marker','.','markersize',1)
% xlabel('x');ylabel('y');
 end

