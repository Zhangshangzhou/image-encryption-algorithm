function [xi,yi,zi,wi]=siweihundun(x0,y0,z0,w0,a,b,c,d,e,r,f,num)
tspan=(0:20/(num+3000):20);
K0=[x0,y0,z0,w0];
funn=@(t,X)[a*(X(2)-X(1)-X(4))+b*X(2)*X(3);
    c*(4*X(1)+X(2))-X(1)*X(3);
    d*X(1)-e*X(3)+X(1)*X(2);
    r*X(1)+f*(3*X(2)*X(3)+X(2)*X(2))];

[~,X]=ode45(funn,tspan,K0);
xi=X(:,1);
yi=X(:,2);
zi=X(:,3);
wi=X(:,4);

end