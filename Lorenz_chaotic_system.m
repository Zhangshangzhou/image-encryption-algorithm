function [x,y,z]=Lorenz_chaotic_system(x_0,y_0,z_0,num)
%设定求解的时间范围和初�?
tspan=(0:500/(num+3000):500);
vec0=[x_0,y_0,z_0];
%用库函数求解(默认精度�?
[t,v]=ode45('lorenz',tspan,vec0);
%拆分结果变量
x=v(:,1);
y=v(:,2);
z=v(:,3);

%  plot(x,z);
%  xlabel('x');
%  ylabel('y');
%  zlabel('z');
end

