function jiema
clc; clear
Q7=imread('Q7.bmp');
y_01=0.0742;
y_02=0.1133;
y_03=0.1016;
y_04=0.0742;
y_05=0.0820;
y_06=0.1094;
y_07=0.0586;
y_08=0.0703;
key=2;
%% zigzag混沌
[M,N]=size(Q7);
SUM=M*N;
Q6=reshape(Q7,1,SUM);
zigzag1=Q6(1:SUM-key+1);
zigzag2=Q6(SUM-key+2:SUM);
zigzag3=[zigzag2 zigzag1];
Q5=izigzagScan(zigzag3,M,N);

%% DNA
t=4;
block=N/t;  %e表示每一行可以分为多少块
r_block=(M/t)*(N/t);
u=3.9999;     %Logistic参数μ，自定为3.99
logi2=zeros(1,M*N+1000);        %预分配内存
logi2(1)=y_02;
for i=1:M*N+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    logi2(i+1)=u*logi2(i)*(1-logi2(i));
end
logi2=logi2(1001:length(logi2)); 
p=mod(round(logi2*10^4),256); %（round四舍五入）（mod除后的余数）
bian_R=reshape(p,M,N);  %reshape转成M行N列的随机矩阵R

%[xi,yi,zi,wi]=siweihundun(y_05,y_06,y_07,y_08,80,45,22,5,21,100,8,r_block);
%xi=xi(3002:length(xi));        %去除前3001项，获得更好的随机性（求解陈氏系统的子函数多计算了3000点）
%yi=yi(3002:length(yi));
%zi=zi(3002:length(zi));
%wi=wi(3002:length(wi));

%xi=mod(round(xi*10^4),8)+1;  %（round四舍五入）(mod求余数)
%yi=mod(round(yi*10^4),8)+1;
%zi=mod(round(zi*10^4),4);
%wi=mod(round(wi*10^4),8)+1;

xi=[4;8;8;2;5;8;6;7;5;5;8;6;8;2;7;8;6;2;8;7;8;3;7;2;8;7;7;2;6;5;5;2;7;6;1;4;1;8;8;8;5;7;3;2;8;8;8;4;6;6;1;1;7;4;3;6;8;3;8;2;5;8;1;8;6;2;7;5;4;4;7;1;5;7;6;2;3;6;4;2;8;2;6;6;5;5;4;6;1;2;4;6;7;8;5;3;6;3;4;8;8;1;6;7;1;5;8;3;7;7;8;2;2;8;4;5;6;2;7;7;3;2;8;7;2;1;2;4;8;8;6;5;3;8;8;3;6;8;8;4;5;6;8;7;1;8;2;7;4;6;2;5;1;1;2;5;5;3;8;2;5;4;7;3;1;5;5;6;3;7;2;8;5;6;8;8;2;7;3;8;6;7;4;7;7;4;3;5;4;3;5;4;4;7;1;4;3;7;4;6;4;1;1;4;5;7;2;2;3;7;6;8;6;7;6;2;3;4;4;7;2;2;5;6;8;3;4;6;1;7;6;6;5;1;6;1;1;5;1;3;2;2;4;2;1;7;6;4;3;1;1;6;8;4;2;2;5;1;4;1;8;5;3;3;3;7;3;6;1;3;7;6;6;1;5;3;1;3;1;5;3;1;2;5;5;5;6;2;5;2;2;4;3;4;3;8;5;5;2;4;2;4;6;5;7;2;1;7;7;2;2;1;6;2;1;5;6;5;2;4;4;2;3;2;6;6;7;7;8;4;4;2;1;2;4;3;8;3;1;3;7;2;6;3;2;3;7;2;6;4;6;4;3;6;1;4;3;4;8;3;7;8;2;7;7;4;4;2;2;6;5;2;5;4;4;6;4;8;2;8;8;7;5;2;2;2;3;8;3;1;8;5;2;2;7;4;7;7;5;5;8;1;6;4;7;1;3;4;6;6;5;3;8;3;8;7;1;6;1;8;6;8;6;8;5;7;2;5;3;5;7;7;3;5;6;3;6;3;3;6;3;5;6;3;2;8;2;6;4;5;5;6;5;4;1;8;4;1;4;2;4;2;6;5;8;3;8;6;8;7;4;1;3;4;8;5;2;1;2;3;8;2;5;7;8;1;1;1;7;1;6;2;4;1;5;2;1;1;2;6;8;6;3;2;3;5;1;3;1;4;6;4;5;7;6;2;4;6;5;4;8;2;6;6;3;4;6;4;3;8;7;5;8;2;5;4;6;3;3;4;6;1;6;8;6;2;8;6;4;7;1;8;8;2;2;3;5;7;6;8;7;4;4;1;7;8;6;3;8;1;4;4;7;4;2;1;6;3;6;2;6;8;1;6;5;4;1;6;2;3;1;2;2;3;5;4;6;2;8;1;4;1;1;3;6;8;8;2;1;2;4;8;7;3;7;5;1;4;5;6;5;6;1;8;6;6;3;5;4;2;6;6;3;8;4;3;6;1;1;8;3;8;5;6;4;8;5;7;2;2;4;3;6;1;3;4;4;1;2;2;2;1;6;6;4;8;3;4;1;8;4;2;2;6;1;4;6;3;7;5;5;2;3;3;1;5;5;4;8;6;2;2;8;6;2;5;2;8;8;5;7;4;8;1;1;1;6;6;6;2;5;6;2;2;8;7;2;6;4;3;2;2;2;1;5;7;8;1;4;8;6;5;2;1;5;1;7;6;5;8;7;4;2;6;8;2;8;2;5;8;6;3;2;3;7;2;1;1;7;3;3;4;2;7;3;2;2;4;5;6;8;8;6;4;2;6;4;4;5;4;3;4;5;5;1;8;7;6;5;8;7;6;8;7;7;6;1;1;2;3;5;2;4;7;3;3;6;1;1;6;5;6;1;6;3;3;5;1;5;2;1;3;2;2;7;1;3;3;2;6;2;6;7;6;1;2;1;4;2;8;4;6;3;1;2;1;2;5;2;1;4;7;6;4;2;4;1;2;5;4;8;5;8;1;1;7;5;1;4;5;1;5;7;2;5;3;5;3;7;1;6;7;5;5;3;8;2;6;2;5;7;2;2;6;1;7;4;2;3;6;7;8;5;8;3;8;1;8;8;2;6;2;7;4;6;6;5;4;7;6;6;4;2;7;1;5;2;7;5;4;4;4;5;6;3;5;8;1;8;3;6;1;4;4;2;1;2;7;2;7;2;5;1;2;6;8;7;7;4;2;8;2;5;8;3;7;2;2;7;1;2;4;8;6;8;4;6;3;3;6;1;4;8;8;1;6;6;8;2;2;7;3;7;7;6;5;1;6;3;3;4;4;3;3;7;2;1;1;3;5;5;5;4;7;3;4;6;1;3;5;2;8;7;4;7;5;6;6;7;6;3;2;2;3;6;7;1;7;2;2;1;2;3;7;1;7;8;2;2;7;6;7;7;1;1;6;1;2;1;8;2;2;1;2;3;5;5;2;5;1;4;6;7;2;7;1;8;3;5;8;6;6;3;1;3;5;6;5;6;2;1;7;4;3;8;6;3;1;7;7;5;7;2;1;5;8;7;7;5;3;1;2;8;3;3;7;2;3;7;8;5;1;2;5;8;6;4;1;5;3;3;7;3;5;1;8;6;4;1;3;1;4;6;3;7;5;6;5;8;2;5;4;1;1;8;2;3;4;8;1;5;3;8;4;5;2;6;8;8;4;3;8;6;8;7;4;4;1;7;4;4;1;2;2;7;7;3;4;1;4;7;1;5;1;7;7;2;6;3;5;4;5;6;2;7;7;1;7;5;7;7;4;4;5;6;4;5;1;5;5;2;4;8;5;2;3;5;5;4;8;2;6;2;7;3;7;5;3;6;8;7;5;8;8;6;7;7;3;7;6;2;5;4;1;7;7;7;8;8;1;5;8;6;6;5;2;2;4;7;5;3;8;6;5;1;8;2;5;7;5;4;2;6;6;4;1;6;4;3;4;6;8;8;6;8;7;5;5;3;8;5;3;2;8;3;6;8;7;6;8;6;8;2;3;8;1;1;2;4;7;5;4;6;5;4;2;2;7;5;8;5;3;3;6;7;5;2;6;7;6;2;4;1;1;8;1;5;8;7;1;4;3;3;8;3;2;3;2;5;6;8;3;6;3;6;5;7;6;1;7;2;5;1;8;4;6;6;8;1;6;7;4;5;4;7;5;7;3;6;8;2;2;1;1;3;2;6;8;1;2;2;3;5;7;5;2;5;5;2;3;3;3;5;4;7;5;6;6;8;2;6;6;1;8;3;4;3;2;2;3;6;8;8;8;3;6;8;3;4;7;5;8;1;4;7;4;6;2;7;7;3;1;8;6;7;6;6;3;2;1;4;2;3;2;7;6;6;8;6;3;7;5;6;1;4;8;1;8;4;8;5;7;5;5;5;3;8;8;4;7;6;2;8;5;3;4;2;3;2;6;7;7;8;1;5;1;6;2;1;6;3;4;3;8;5;3;1;5;7;5;2;2;7;3;5;7;7;8;7;2;3;1;2;5;1;2;5;1;2;7;1;4;3;1;1;2;2;3;4;8;6;5;4;2;7;2;6;6;2;3;4;4;1;6;5;3;4;5;8;8;3;3;8;6;2;5;8;1;2;4;6;7;8;5;4;1;3;3;3;4;3;3;3;2;4;1;2;7;7;4;7;4;8;5;3;1;7;3;1;8;6;3;8;1;7;1;3;1;6;1;7;4;2;7;3;6;2;3;3;2;4;4;8;7;7;7;7;3;8;3;4;6;3;3;3;1;5;5;2;4;4;5;4;2;3;2;8;3;1;8;6;8;6;3;5;2;8;2;5;4;3;6;1;7;7;4;4;5;3;7;7;5;7;5;8;8;2;7;5;7;6;7;5;7;5;6;8;7;6;6;7;3;5;8;3;8;1;5;4;6;7;3;8;3;2;5;4;2;6;2;6;2;8;5;2;2;7;8;4;5;5;8;6;3;6;7;7;8;2;8;2;2;5;1;4;1;3;1;4;2;3;7;2;6;3;3;5;8;3;8;2;5;5;4;2;6;8;4;7;3;8;6;1;2;1;4;3;3;2;5;1;8;6;6;4;8;2;1;8;4;2;1;4;1;7;2;7;2;6;8;4;3;1;3;2;3;3;4;5;8;2;6;7;1;4;2;8;4;1;6;3;2;1;3;8;8;6;5;7;4;7;4;7;7;2;8;7;5;2;7;5;5;5;5;1;3;1;5;2;3;4;5;3;5;3;7;3;6;7;7;1;5;8;1;1;2;3;2;8;8;7;3;4;7;5;8;5;3;6;4;3;3;6;1;2;3;8;3;1;3;5;7;6;5;4;2;6;1;3;6;2;2;3;3;3;2;5;2;4;8;3;4;5;5;6;3;6;8;3;6;7;4;1;1;6;4;5;4;8;3;2;1;4;6;7;8;8;1;4;8;6;1;6;6;7;8;1;5;1;6;7;7;6;4;4;8;8;1;8;2;8;6;7;6;3;5;1;1;5;7;3;7;5;8;3;5;8;6;3;3;2;3;7;2;8;4;1;5;8;5;7;2;7;1;8;1;7;2;2;7;2;2;5;1;4;3;6;5;7;4;5;6;7;7;7;3;1;7;8;1;1;3;4;5;7;6;4;4;7;8;7;6;5;6;1;8;1;7;8;8;2;3;2;7;5;7;8;7;2;3;4;5;3;8;5;7;3;2;6;1;4;4;4;7;4;4;7;1;1;8;2;8;5;2;6;5;6;8;5;4;1;7;2;5;6;2;3;6;4;3;8;7;4;1;6;4;1;8;4;2;1;8;3;8;1;4;3;6;2;5;2;2;8;8;2;3;1;3;4;1;2;8;3;7;2;6;6;6;8;6;6;2;5;7;5;5;6;3;3;6;7;4;2;3;2;2;3;3;5;5;1;5;2;1;3;4;2;7;6;7;5;8;4;7;1;6;8;4;5;7;6;1;6;4;1;3;1;7;2;3;2;5;6;3;1;4;7;1;4;8;2;3;5;6;1;4;8;2;6;6;4;2;1;1;2;6;3;7;8;6;6;1;6;4;2;4;2;7;2;7;1;3;8;3;1;4;8;8;5;4;8;1;6;5;4;3;6;7;2;6;4;8;2;7;1;5;3;6;8;3;1;8;6;6;5;2;8;2;1;8;5;1;5;7;8;1;1;5;4;3;3;2;2;8;2;7;6;8;4;7;7;5;2;4;5;1;1;1;7;4;4;1;7;4;5;6;1;1;7;2;2;3;1;5;6;2;2;7;5;5;2;4;3;7;2;8;1;2;5;2;4;7;3;8;2;4;4;8;8;6;8;8;6;8;2;3;7;8;2;7;2;6;3;6;7;4;4;1;3;6;5;1;4;4;1;4;2;7;2;6;8;3;3;4;5;2;6;6;5;1;3;1;1;7;4;5;2;1;7;3;5;7;6;6;1;5;5;4;6;2;3;4;6;1;5;3;5;8;5;7;4;7;2;5;8;5;4;4;8;2;8;4;4;5;7;1;7;4;1;3;2;6;8;2;6;8;4;4;7;6;7;8;3;8;5;8;1;4;1;4;7;7;5;8;7;6;8;5;3;2;5;2;2;2;2;6;8;5;4;5;4;4;2;5;5;1;3;7;4;3;7;4;3;7;3;7;6;2;2;2;1;4;1;6;1;3;3;2;5;1;3;1;7;3;4;5;1;8;7;5;8;3;7;3;5;3;7;6;4;5;5;7;8;4;1;3;6;8;2;1;4;6;4;3;7;6;8;2;8;4;7;4;5;6;2;4;1;7;4;4;1;2;2;6;7;8;3;1;4;6;2;8;4;1;7;6;6;7;8;6;2;4;2;6;6;8;4;1;7;3;4;2;8;5;7;6;7;4;8;7;5;3;5;8;3;5;4;5;6;7;2;2;3;3;4;1;2;8;2;5;3;4;5;4;1;8;8;1;5;3;2;4;7;4;2;3;6;5;2;4;3;3;1;3;1;6;6;4;1;4;4;6;1;5;2;5;5;1;3;3;6;4;2;2;2;2;1;4;8;5;4;6;3;4;8;5;2;6;8;3;4;7;1;7;7;3;5;2;2;4;1;1;6;7;8;6;2;3;5;2;4;8;4;4;8;8;2;6;6;4;8;1;7;2;4;1;8;4;1;1;3;7;2;7;2;2;2;5;2;1;5;2;7;6;8;5;1;2;7;1;4;2;4;3;4;3;2;2;1;4;1;7;7;8;8;3;3;5;8;5;7;6;6;1;5;1;7;7;1;4;1;5;2;2;3;5;5;7;6;7;2;2;4;5;3;8;4;4;3;8;2;5;5;4;1;5;5;7;1;2;5;5;2;1;4;8;2;1;4;3;4;3;3;2;6;7;5;2;2;3;6;5;7;1;6;2;6;8;1;1;2;3;5;3;6;3;5;3;6;4;5;4;6;5;4;2;7;1;3;4;5;6;1;7;1;2;5;1;1;6;3;2;1;8;3;4;7;7;5;3;1;7;4;1;5;5;6;7;5;7;3;7;6;3;6;6;7;2;6;8;4;3;2;1;1;5;3;6;6;2;7;2;7;7;4;4;5;5;3;6;3;5;1;4;5;6;5;4;7;2;5;6;5;8;8;3;6;1;4;3;3;8;3;6;8;7;5;3;4;6;1;6;1;2;7;5;5;1;8;4;3;5;7;3;8;3;6;6;1;2;1;8;8;3;6;8;3;6;5;3;7;3;2;7;6;3;2;2;4;8;4;1;7;1;6;5;2;6;3;6;6;5;6;1;3;1;8;4;2;1;8;4;5;2;4;1;4;7;1;2;4;1;8;8;4;5;8;7;4;5;4;3;2;5;3;1;5;4;4;1;5;7;6;7;5;7;1;1;1;5;3;7;1;6;1;8;2;8;1;3;2;8;1;7;8;6;4;5;3;3;8;4;8;5;3;2;4;1;3;8;2;1;8;8;2;6;8;1;5;7;6;8;1;2;4;5;6;6;5;3;3;2;6;6;1;5;6;5;3;1;2;8;2;8;1;1;1;6;4;7;4;2;2;7;7;1;6;6;2;7;2;3;7;3;7;4;6;6;2;6;3;1;4;5;2;6;8;4;5;6;4;2;8;3;1;8;2;3;1;3;8;6;7;6;3;2;2;5;7;4;2;7;2;5;2;7;7;2;4;4;3;6;3;2;8;5;5;2;3;7;6;7;8;6;8;8;4;1;5;2;3;3;6;7;3;7;7;3;5;3;3;2;1;2;3;1;5;5;2;5;5;8;3;5;4;7;2;7;3;7;6;7;6;5;5;4;8;6;4;3;5;2;7;8;8;5;8;6;8;7;5;5;5;8;6;8;1;3;7;4;2;5;8;7;5;2;4;6;2;4;5;6;2;8;2;3;4;6;2;1;4;3;8;4;1;6;7;1;5;6;4;1;4;6;4;5;2;2;5;1;8;6;6;5;5;8;5;2;7;4;7;5;7;3;2;1;5;4;5;2;5;6;2;3;5;1;5;5;7;7;2;8;7;5;2;8;6;4;1;8;7;6;6;8;8;5;3;3;2;2;1;5;8;6;3;8;7;7;2;4;5;7;5;7;7;3;3;3;7;6;8;2;1;4;8;1;2;8;4;8;3;2;8;5;2;7;6;3;6;5;8;7;8;6;1;8;7;4;3;7;8;7;4;2;8;8;7;1;2;1;5;6;7;1;5;7;4;3;2;8;5;6;2;4;3;1;5;4;2;1;3;1;7;4;3;3;3;4;7;5;1;4;2;5;2;3;5;2;4;6;1;6;2;4;7;3;3;6;3;7;7;7;8;6;7;4;8;1;2;8;4;1;7;2;2;7;4;6;4;5;5;3;8;5;6;7;3;6;7;3;7;2;8;6;1;5;1;3;4;7;4;1;2;4;2;8;6;4;4;3;6;2;5;4;6;4;4;8;4;5;1;5;2;1;6;1;3;7;8;7;7;6;6;1;4;4;6;7;4;3;5;5;3;4;5;7;5;7;4;2;6;2;2;7;4;7;8;8;7;4;7;6;2;3;1;4;3;8;6;6;7;8;3;1;6;8;7;5;3;1;1;8;2;3;7;8;6;4;4;5;1;8;5;5;8;7;4;5;8;2;8;5;1;3;4;5;7;4;7;2;8;1;2;3;7;1;4;6;3;2;8;3;3;2;1;2;2;2;3;2;5;6;5;7;2;5;8;3;3;8;1;1;1;8;2;2;3;5;5;5;5;6;5;8;8;4;3;4;7;5;7;5;4;8;5;5;2;2;6;1;7;7;5;5;4;7;3;8;3;5;2;4;1;4;2;6;8;5;2;3;7;3;8;6;2;3;7;7;8;4;3;2;8;4;3;1;7;5;4;7;1;4;6;7;4;3;6;1;4;5;7;6;8;8;7;5;2;4;5;1;6;3;8;1;4;8;5;6;4;2;2;1;3;6;2;5;4;1;2;2;8;3;2;2;7;3;7;8;2;3;3;2;7;8;2;8;3;5;4;5;5;5;2;4;6;1;7;2;3;5;7;7;4;8;5;2;3;6;4;6;5;3;4;7;7;2;7;8;5;5;7;3;6;7;7;2;2;2;7;3;2;6;3;5;1;6;1;2;3;8;7;8;6;7;7;3;2;5;3;4;8;1;6;3;4;5;6;3;6;4;7;1;3;5;7;3;3;6;8;2;1;3;5;5;4;8;2;6;4;3;6;4;3;1;7;2;6;7;2;5;6;5;6;8;6;6;2;7;4;6;1;3;2;7;2;7;2;1;6;4;4;8;6;8;2;7;2;7;6;5;8;7;5;5;1;2;6;8;6;6;7;7;1;6;6;1;7;7;4;4;3;6;6;7;2;5;7;5;5;8;7;2;2;4;2;3;8;3;3;7;6;3;8;2;8;6;8;1;7;3;1;8;7;5;1;7;2;2;6;5;2;2;7;7;7;8;5;7;7;7;6;2;6;3;7;3;1;4;2;3;5;7;1;6;7;3;8;8;3;3;8;7;7;2;8;8;3;6;2;8;6;1;5;1;5;2;5;5;5;8;8;6;5;1;5;3;8;2;2;4;2;1;2;3;7;3;3;8;7;7;7;5;7;6;4;4;2;3;2;3;2;3;7;2;8;4;5;6;5;4;8;4;2;7;2;7;5;8;5;2;1;7;2;1;7;6;1;5;7;7;1;7;1;5;2;8;2;3;3;2;1;8;6;4;2;8;2;5;4;6;2;5;7;4;3;5;1;6;7;7;8;8;1;7;8;2;7;4;1;1;2;1;6;3;1;8;3;5;6;8;8;6;5;2;8;8]
yi=[3;7;3;1;4;2;2;4;2;7;6;3;6;4;7;3;3;8;7;2;4;8;7;8;3;8;3;7;1;3;5;1;2;5;4;3;7;5;8;5;5;6;4;2;2;7;8;4;3;5;6;3;2;5;4;6;3;5;2;3;6;7;6;2;5;7;4;7;1;1;6;7;1;5;7;3;1;5;2;3;3;4;1;4;2;3;6;2;1;5;2;1;8;5;4;5;1;7;6;8;6;2;5;3;6;2;6;8;1;2;6;1;2;1;3;7;2;6;2;6;4;7;6;8;5;3;3;2;1;4;5;2;4;7;3;6;4;8;1;7;1;7;3;8;6;6;4;4;7;3;2;2;7;1;5;5;6;8;6;4;2;1;1;6;1;1;6;1;1;4;4;1;1;8;5;1;1;3;7;1;8;3;8;2;4;6;7;4;8;7;3;1;4;7;6;3;2;7;1;8;7;7;8;8;6;6;6;3;2;7;3;3;4;7;5;7;1;1;1;3;6;8;5;5;3;1;4;2;1;8;2;5;5;6;5;2;5;8;3;4;2;6;8;1;3;7;6;7;5;1;2;6;2;4;3;5;6;4;5;3;1;3;7;4;1;5;3;5;4;2;7;4;2;7;6;7;3;6;2;8;4;6;6;8;2;6;4;4;4;7;3;8;8;4;8;3;6;3;5;6;8;5;8;2;2;5;2;7;5;7;6;5;6;5;5;5;6;3;3;6;6;2;4;3;1;7;3;4;6;5;2;1;4;1;6;1;8;6;7;6;5;1;1;7;3;2;2;8;5;1;7;2;8;4;2;1;5;3;3;5;2;2;8;2;3;3;8;7;7;4;7;4;5;1;5;1;3;2;3;2;6;3;4;4;2;5;5;1;5;1;7;8;3;7;1;4;3;4;2;1;1;5;6;2;1;4;2;8;7;3;3;5;5;2;3;3;2;1;8;4;5;2;7;3;1;7;7;2;3;2;4;6;6;5;8;8;6;1;1;8;6;3;7;2;7;6;7;6;2;6;1;3;2;1;8;2;8;4;6;8;3;1;6;2;7;6;6;8;1;2;4;5;3;1;3;2;4;1;4;3;2;5;7;7;3;2;4;1;3;8;1;6;1;8;3;2;1;5;4;5;7;6;7;8;3;7;6;3;5;4;7;5;5;6;4;7;1;2;3;2;5;3;5;3;7;7;4;6;1;4;5;8;5;5;3;2;1;2;2;4;5;2;8;1;4;7;3;2;4;5;1;1;7;2;7;6;8;4;4;4;5;2;8;1;5;6;2;7;2;4;8;6;8;2;3;2;8;3;1;1;4;4;8;2;4;1;5;5;3;3;5;7;4;2;6;2;7;3;7;2;4;1;8;4;2;2;3;2;2;1;7;5;4;7;2;3;7;4;1;7;3;3;1;1;4;8;7;3;1;5;6;8;1;5;4;3;8;6;6;1;3;4;7;5;7;1;2;1;4;2;7;6;7;2;5;5;2;8;6;1;8;5;7;1;4;1;6;8;3;1;2;8;6;4;7;5;6;3;5;4;5;6;3;8;7;3;7;7;8;7;6;3;1;3;4;3;3;3;5;7;3;8;5;3;3;6;8;4;6;6;4;2;3;4;3;4;5;7;3;1;2;5;3;2;8;2;6;5;8;1;8;4;3;5;8;3;2;7;8;3;1;5;4;2;3;7;4;2;8;3;1;1;8;3;5;6;2;2;1;1;1;3;2;7;3;1;4;5;7;8;2;3;1;2;6;7;1;1;4;4;7;4;8;5;6;4;1;3;2;6;6;6;1;7;2;8;3;8;6;8;4;3;8;5;2;2;2;7;4;5;7;7;3;3;4;1;1;4;3;6;3;1;4;7;1;1;5;6;6;1;4;7;2;5;2;4;7;7;4;8;8;3;6;7;5;2;7;2;4;7;7;7;8;4;8;5;6;8;8;6;1;1;3;8;1;5;5;3;3;5;6;5;8;7;2;1;2;7;6;2;7;1;5;1;3;7;7;6;6;6;5;2;3;6;2;2;6;4;1;8;4;5;3;2;6;4;2;4;6;6;7;2;6;3;1;1;8;8;6;6;4;3;8;5;6;5;7;6;1;7;2;6;6;3;5;2;8;4;3;1;3;5;8;4;8;6;7;4;3;8;1;3;4;7;1;8;2;8;2;4;2;6;2;4;5;4;1;3;4;4;6;8;1;8;1;1;2;3;3;3;2;8;6;3;6;8;8;1;5;5;8;6;7;2;2;4;8;6;4;3;6;2;5;2;1;4;6;3;6;5;7;8;2;1;4;8;7;7;4;7;5;8;4;4;7;4;5;7;5;7;5;4;5;5;2;6;8;6;6;2;3;7;5;7;5;2;8;5;6;6;2;5;6;8;2;8;4;7;7;7;6;5;4;6;5;1;6;6;3;6;3;1;7;7;6;2;5;2;6;2;1;1;3;3;6;5;3;6;4;5;2;3;1;1;1;8;8;8;6;5;4;8;2;2;8;5;7;6;4;2;2;4;1;5;1;1;2;4;1;2;7;4;4;6;8;6;8;3;5;2;5;2;3;2;8;3;3;3;7;6;3;3;2;3;1;7;1;3;8;7;8;7;2;7;1;8;4;4;6;6;8;3;5;6;8;7;1;4;2;8;8;6;7;5;3;2;5;1;6;5;4;5;6;7;5;5;6;6;4;1;2;2;5;4;6;3;3;8;1;1;2;7;3;8;1;6;4;7;4;3;6;7;1;5;4;7;5;2;3;4;6;1;4;3;7;1;4;5;2;8;6;1;7;5;8;7;8;6;6;5;7;1;1;3;6;6;4;1;2;7;7;8;7;5;7;3;8;7;6;3;1;4;1;3;2;3;3;2;8;3;6;1;3;1;2;7;7;6;2;7;6;7;1;5;7;7;2;6;6;2;8;8;7;5;7;3;2;5;1;1;5;6;6;7;1;8;3;7;2;6;4;7;6;7;7;8;2;2;5;1;7;1;2;5;2;4;5;3;8;4;3;4;6;6;3;3;1;3;4;1;8;2;5;1;7;6;5;6;7;6;6;8;3;4;5;4;4;3;2;5;2;4;3;8;7;3;7;4;5;2;5;4;1;6;8;7;5;2;5;8;8;1;8;3;6;7;2;1;3;6;1;8;3;6;1;7;7;4;4;8;2;3;3;4;2;3;6;7;1;3;2;2;6;4;5;4;4;5;8;4;8;5;8;7;4;1;5;6;3;5;6;5;3;5;4;1;8;3;5;2;2;4;3;4;8;7;6;1;1;6;7;7;5;5;4;3;6;2;2;7;8;3;4;2;7;1;3;7;4;8;7;1;7;5;3;4;2;5;1;3;1;6;5;5;8;2;5;3;4;8;3;5;8;4;7;1;7;4;4;1;6;7;5;4;7;8;7;2;1;3;4;3;6;7;7;8;6;8;8;5;2;5;6;3;7;5;3;6;3;1;7;6;8;4;7;8;6;3;5;3;8;5;5;7;1;2;8;6;6;1;2;5;1;7;3;7;7;2;6;2;2;6;6;8;5;5;7;8;8;3;7;4;1;3;5;7;2;1;7;8;8;8;8;7;1;7;4;2;7;6;5;1;6;8;1;3;5;5;5;3;1;6;7;1;5;4;3;1;2;3;7;4;8;6;2;6;8;5;1;3;1;7;5;8;4;2;5;3;1;3;5;7;4;6;6;6;5;7;3;1;7;1;5;4;2;1;4;3;3;2;2;5;1;2;1;5;7;7;6;6;5;5;4;4;4;5;2;2;5;7;4;7;5;2;8;4;7;5;2;1;2;7;8;2;8;2;6;6;7;1;4;5;7;7;3;6;2;4;8;6;3;2;2;4;2;4;6;1;4;8;6;7;3;7;7;4;5;7;8;3;7;3;6;5;8;6;8;4;5;4;2;4;2;6;2;8;5;6;6;8;7;5;2;1;1;4;3;8;5;5;4;5;8;7;1;2;2;4;4;4;5;7;5;1;6;1;4;2;6;4;8;2;7;3;6;7;7;6;2;8;3;6;5;8;7;1;5;5;2;8;8;5;2;5;2;2;7;6;7;1;1;7;7;5;3;1;8;6;1;1;1;3;1;4;1;8;4;1;6;6;7;7;5;8;5;4;8;1;1;8;1;5;5;3;6;1;2;7;7;3;6;7;6;4;5;1;7;4;7;8;1;7;7;4;7;5;6;6;8;7;5;3;1;2;2;1;7;1;6;2;2;5;8;3;1;1;2;7;1;8;8;3;4;2;3;2;2;6;3;8;3;1;3;2;3;5;3;3;7;3;1;2;3;8;3;7;6;7;6;7;6;2;3;5;3;5;6;7;8;1;1;3;2;5;4;4;4;2;8;2;7;1;6;5;3;4;4;6;1;2;6;4;7;4;2;2;1;8;3;5;3;5;2;2;3;7;7;5;3;1;8;2;6;7;2;7;4;2;8;5;3;1;6;4;6;1;6;6;5;1;2;8;8;6;2;7;3;5;4;5;6;8;3;6;4;4;2;8;2;7;8;5;2;2;7;2;5;6;8;5;5;8;8;6;8;7;4;7;2;5;4;6;6;1;4;7;4;8;5;5;1;4;6;8;8;3;4;2;2;4;3;3;1;6;4;2;3;3;7;5;8;1;2;4;7;4;5;2;6;2;3;6;1;2;2;3;2;5;8;2;8;1;6;5;3;8;4;7;6;4;7;1;4;4;4;7;4;4;5;8;8;7;7;4;5;5;6;1;1;6;4;4;1;5;4;4;3;2;1;5;1;8;6;4;4;2;1;5;5;4;4;7;2;7;4;5;6;4;7;2;6;4;2;8;6;2;4;1;7;5;3;8;5;2;1;2;2;2;7;5;7;2;5;4;8;7;3;3;4;4;7;3;5;4;6;4;6;6;6;4;4;2;6;2;8;5;3;7;4;8;4;5;6;7;6;1;6;3;8;7;3;8;5;3;5;3;6;5;2;1;2;7;1;4;7;8;5;8;8;3;2;1;5;2;8;5;1;8;2;4;1;1;7;3;7;8;2;5;8;2;2;5;5;4;1;2;8;8;1;1;3;3;3;7;1;3;6;2;7;8;2;5;7;7;8;5;5;6;6;4;3;7;2;6;6;2;2;4;2;2;8;8;3;1;3;4;4;3;3;1;4;5;4;8;7;3;3;2;1;7;4;6;1;1;1;7;2;3;2;4;6;8;5;2;3;6;7;5;2;6;2;3;8;3;7;1;5;2;2;6;4;4;8;1;8;4;1;8;8;5;1;5;1;8;7;3;3;3;5;8;7;4;8;5;4;4;3;8;7;2;2;5;1;2;7;4;5;3;6;6;3;4;6;2;4;6;5;1;5;7;5;7;2;8;3;1;5;3;1;5;2;1;1;6;5;1;1;6;8;5;3;2;1;2;2;4;6;4;5;4;5;3;6;1;4;1;3;8;1;7;7;4;3;2;8;5;8;5;1;8;1;4;1;5;5;3;2;4;6;3;3;8;3;7;5;5;5;4;5;2;1;7;8;1;6;6;3;8;6;1;7;3;3;3;2;2;4;7;3;8;8;6;7;6;1;7;3;5;7;5;4;3;3;7;7;3;6;8;3;2;1;3;6;6;6;6;3;3;5;3;3;4;1;5;7;4;7;1;8;3;4;1;3;4;2;4;5;3;6;3;2;6;1;8;3;7;7;2;2;8;7;2;8;6;6;1;3;8;3;3;8;2;8;6;8;6;2;7;2;2;4;6;7;1;1;6;5;8;1;8;6;2;2;3;5;3;4;8;7;6;6;3;3;2;4;6;2;7;6;3;6;1;3;8;1;1;7;5;4;3;1;3;3;2;3;3;5;8;3;8;6;5;1;6;7;4;4;7;5;1;4;6;5;5;4;1;2;7;2;1;4;2;7;3;2;8;1;6;5;1;8;6;6;8;5;1;4;6;1;7;5;1;2;4;7;6;8;1;3;1;1;1;2;5;8;8;2;5;2;5;6;2;7;5;8;4;4;2;4;7;7;3;7;3;8;4;6;6;2;1;2;1;7;7;1;3;3;5;1;6;6;7;1;7;1;6;1;1;1;4;6;7;6;4;2;1;8;5;6;5;3;4;1;3;1;8;7;7;3;6;3;6;1;6;5;2;3;7;6;1;3;8;5;7;2;2;3;4;4;4;5;3;2;8;1;1;8;3;4;6;7;6;7;7;5;6;8;5;3;5;1;7;5;5;6;7;5;1;6;4;4;7;8;7;4;8;6;3;5;4;8;5;1;5;5;6;4;5;6;8;8;4;6;6;2;6;4;2;7;3;4;8;6;2;3;3;5;1;6;2;5;2;3;8;3;5;7;2;5;5;4;1;1;4;5;4;7;6;3;3;3;1;6;5;6;8;6;1;7;1;7;1;5;1;5;1;3;7;1;5;8;1;7;1;7;7;2;4;1;3;1;4;4;3;6;8;5;7;4;6;4;7;7;8;4;1;2;3;1;2;8;4;2;8;5;2;5;7;2;1;2;2;1;8;4;2;8;8;5;1;6;5;8;2;8;8;6;4;1;3;7;4;7;2;2;6;1;7;1;8;3;8;5;8;3;4;7;7;5;6;2;3;4;8;7;3;7;1;8;8;6;6;6;2;3;2;8;7;8;1;5;6;2;1;1;7;5;1;6;3;5;1;5;8;3;3;7;2;4;7;6;5;5;5;2;8;7;7;8;2;7;1;3;1;4;8;5;8;2;2;2;3;5;7;5;8;4;8;3;7;4;4;6;6;5;3;4;2;6;3;5;6;6;2;7;1;3;7;4;4;2;7;1;1;1;1;3;2;4;3;5;2;3;1;8;1;2;5;8;3;8;7;4;1;3;5;7;1;5;6;2;6;3;8;7;5;6;3;1;1;5;6;4;6;1;5;3;6;3;7;2;5;7;8;1;8;2;1;7;1;3;7;1;1;1;4;2;4;4;7;1;4;3;7;1;3;5;5;8;2;8;3;7;7;1;4;1;7;4;6;4;2;8;5;1;2;5;6;5;8;6;7;5;6;2;3;1;3;4;8;2;8;7;4;1;5;1;3;5;2;3;8;5;3;3;1;6;6;1;5;4;5;7;4;5;7;8;3;1;2;1;8;1;6;1;2;8;2;1;8;3;3;5;1;2;4;7;5;3;3;1;5;3;2;5;6;4;2;4;2;6;4;1;8;5;3;5;4;1;6;7;6;3;8;5;7;5;1;5;7;2;7;2;4;8;2;8;5;5;2;6;2;2;7;7;2;6;3;2;8;7;2;1;6;6;8;1;4;1;7;3;8;5;8;5;1;4;5;1;5;1;1;2;2;7;7;1;2;7;1;3;4;8;6;6;5;4;6;1;1;7;6;6;5;3;4;1;1;8;7;5;3;3;8;4;2;6;1;5;3;6;3;6;5;5;6;5;8;1;8;7;3;6;3;3;2;1;7;8;2;1;4;8;1;3;5;7;5;7;5;4;5;2;7;7;4;3;1;5;1;3;2;7;1;3;7;2;2;3;3;4;4;2;4;1;8;1;3;5;8;7;3;6;1;6;2;5;6;2;8;5;8;2;2;6;6;6;2;5;6;7;2;5;7;1;8;6;5;4;4;5;6;6;6;5;4;6;5;3;6;3;5;3;3;8;5;1;3;5;4;6;5;1;5;5;1;1;1;6;7;1;8;5;6;3;4;4;8;3;3;1;2;4;5;2;5;3;2;8;7;4;5;2;4;3;1;7;7;4;5;7;2;7;1;7;5;2;3;1;3;4;1;2;5;4;6;4;8;6;1;8;8;5;8;1;7;1;7;4;5;6;7;7;1;2;1;5;4;5;8;5;3;7;2;6;8;4;4;3;5;3;1;4;2;4;1;4;8;5;8;6;1;6;8;2;3;2;2;7;1;5;3;3;4;8;4;7;2;8;3;3;4;2;4;8;5;7;2;8;8;5;5;4;2;8;7;8;5;3;6;2;7;1;6;5;1;5;4;1;8;3;2;3;6;2;4;7;6;8;6;4;1;5;3;7;2;7;2;4;1;5;5;6;7;3;8;5;6;7;3;7;6;1;4;2;6;2;8;8;6;2;2;1;1;1;4;6;5;4;1;2;8;5;3;1;5;2;2;2;1;1;7;5;7;7;3;5;3;4;4;6;3;4;8;1;5;5;6;1;6;7;1;2;6;7;3;6;2;7;7;3;8;3;1;5;2;8;2;1;5;5;4;3;7;1;4;7;5;7;7;7;6;6;4;7;4;7;1;5;1;1;7;5;6;2;2;6;6;8;7;4;4;7;3;6;3;4;6;4;2;7;5;1;1;4;7;1;4;3;6;5;1;2;4;8;7;1;5;6;8;4;6;6;5;1;4;7;8;6;7;7;1;4;1;2;4;8;7;8;2;7;7;4;6;5;8;5;1;5;3;8;8;8;6;7;3;5;1;8;4;5;6;7;5;6;8;6;3;1;8;2;4;3;7;4;1;5;5;1;5;2;6;7;3;4;6;3;2;7;2;6;7;2;2;1;1;3;2;6;4;3;8;8;1;5;5;8;1;2;6;4;8;8;8;2;2;3;6;3;6;2;1;2;1;7;5;3;4;7;3;3;2;8;2;3;3;3;5;7;3;8;1;6;8;5;1;4;5;7;6;4;7;8;1;4;7;8;3;4;6;4;3;3;7;5;5;5;3;6;1;3;4;8;4;2;8;4;8;2;2;3;3;4;1;4;5;7;4;2;8;4;3;5;4;1;3;3;8;1;5;2;1;4;3;6;3;5;3;4;5;3;7;7;6;3;3;1;6;6;5;8;8;3;6;3;6;8;1;2;4;5;5;7;5;2;6;5;3;2;8;4;5;2;7;4;3;3;4;1;3;8;8;7;1;8;2;8;8;5;5;6;1;2;6;8;1;4;3;7;3;1;3;5;3;3;5;2;5;1;4;8;4;8;4;4;3;7;3;2;4;2;2;8;2;7;1;3;6;3;6;6;5;8;1;3;2;1;2;4;7;1;6;3;5;2;6;2;8;8;3;1;8;4;1;4;7;3;5;3;4;5;1;2;2;8;1;3;5;7;2;1;2;1;5;6;4;2;4;6;2;8;8;1;8;6;7;1;4;4;2;1;5;8;7;1;6;3;6;3;6;4;5;7;3;8;1;1;4;8;2;6;6;7;6;8;8;5;1;5;4;2;4;8;5;4;3;5;6;1;6;4;7;1;8;8;1;1;4;7;6;3;3;6;6;5;8]
zi=[1;0;3;2;1;1;1;3;1;0;0;2;1;0;2;2;1;0;3;0;1;0;3;0;2;1;2;3;1;3;1;0;0;0;3;1;0;2;1;1;2;0;2;2;0;1;2;1;2;3;2;1;1;0;2;1;2;2;1;1;3;0;1;2;2;1;2;3;1;0;2;2;1;0;0;1;3;2;3;1;0;0;1;1;0;1;1;0;0;3;1;3;3;1;2;0;0;2;2;3;3;2;0;3;2;2;3;0;3;1;1;0;2;1;0;3;2;3;3;1;3;0;2;2;3;0;0;1;1;0;2;3;0;3;3;3;3;3;1;1;0;2;3;1;1;2;2;3;1;0;0;1;1;0;2;0;0;1;2;3;3;3;3;0;1;0;0;2;3;2;0;1;2;1;1;2;2;0;3;3;1;0;1;0;0;2;1;3;2;2;1;0;3;1;3;2;3;3;1;0;1;2;3;1;1;1;3;3;0;1;1;0;0;2;2;1;1;0;3;1;0;0;0;0;1;0;0;2;2;0;0;0;3;2;3;3;2;2;2;3;0;3;3;0;0;2;0;3;1;0;1;2;3;0;1;2;0;0;3;2;3;2;3;1;0;0;3;0;3;3;1;1;3;0;3;1;1;2;2;0;3;3;0;3;0;0;0;3;3;0;2;0;3;1;0;3;3;1;3;3;0;2;0;2;3;1;1;2;3;3;1;3;0;0;0;2;3;3;3;0;3;1;1;1;3;0;0;0;0;2;3;0;2;0;3;2;1;2;3;3;1;3;1;2;0;1;3;1;3;2;0;1;2;1;0;3;1;0;1;3;2;0;1;2;1;0;2;0;1;2;0;2;2;2;2;0;1;0;3;1;2;1;2;1;0;2;0;0;2;1;2;0;0;3;3;3;1;3;0;3;1;2;3;2;3;1;0;2;2;2;2;0;2;3;3;0;3;0;0;2;0;3;1;2;1;0;0;1;1;1;2;0;1;1;2;1;1;0;0;0;0;2;2;0;0;2;3;3;0;0;0;2;1;2;0;1;0;1;0;0;0;2;2;2;0;1;3;1;1;0;2;2;0;0;1;1;3;2;1;3;0;0;0;0;3;1;2;0;1;0;0;0;1;3;2;3;2;1;0;1;3;2;3;0;3;2;2;0;3;2;0;3;0;3;2;2;1;1;3;0;0;0;2;3;1;3;1;3;3;0;1;2;3;3;0;1;0;3;1;1;0;1;2;2;1;1;2;3;3;2;1;0;1;0;2;3;0;3;0;3;1;0;2;0;1;1;2;2;1;3;1;0;3;1;2;0;3;3;2;3;0;2;2;1;0;2;0;3;2;3;2;1;2;2;3;3;1;3;0;0;3;0;2;3;0;1;1;1;3;3;1;2;1;1;0;3;1;0;2;3;3;1;3;0;2;0;0;1;3;0;0;0;3;2;1;1;1;2;1;1;0;0;1;3;1;1;0;2;3;1;0;2;3;2;2;0;0;3;1;2;3;0;2;1;3;1;3;3;1;3;0;0;1;3;0;2;3;3;1;3;0;0;1;1;3;1;2;1;1;0;1;2;2;3;2;0;3;1;0;0;3;1;0;0;3;0;3;1;3;0;3;1;0;0;1;3;3;1;2;2;3;0;3;0;0;3;2;2;1;0;2;1;2;1;1;0;0;2;1;1;1;1;3;2;2;1;3;1;3;3;3;1;0;0;2;3;2;2;0;1;1;1;3;1;0;0;1;3;0;3;1;0;3;3;0;2;3;0;2;1;1;3;2;0;0;2;1;3;0;1;2;1;1;0;1;0;1;0;1;2;0;1;3;1;2;3;1;1;1;1;1;3;1;0;0;3;1;0;2;2;3;0;0;0;1;0;2;2;1;2;0;2;1;3;2;2;0;2;2;3;3;2;0;3;0;2;2;1;1;2;1;2;2;0;1;3;0;1;2;1;0;0;2;3;2;2;0;1;1;0;1;1;2;1;0;3;3;1;1;1;3;1;3;3;1;1;1;1;3;3;3;1;0;3;0;3;2;1;2;3;2;3;0;0;3;0;2;2;1;0;1;0;2;0;1;2;2;2;3;1;0;0;2;0;1;3;2;3;0;0;0;1;1;2;3;1;0;0;2;1;1;3;1;3;0;3;3;1;3;2;0;3;2;1;1;1;3;1;1;1;0;3;3;3;2;2;3;0;0;2;1;0;2;1;2;1;2;3;3;3;0;0;0;1;2;3;3;1;2;0;2;2;3;3;3;3;2;1;1;2;0;3;0;3;3;3;0;1;1;0;2;0;0;3;2;0;2;1;2;3;0;2;0;3;3;3;3;0;3;1;0;0;0;1;3;3;2;1;2;0;3;1;1;1;2;1;2;2;1;0;2;3;0;2;3;1;1;2;0;1;3;0;0;1;2;2;0;2;3;2;1;0;0;2;0;0;1;3;3;1;0;0;0;0;1;1;0;0;0;3;3;0;1;3;1;1;0;0;1;3;0;3;2;2;3;3;2;2;3;1;3;0;1;2;0;0;1;2;2;0;3;1;0;2;0;3;1;3;3;2;3;1;1;1;3;0;0;0;1;3;3;2;1;3;1;1;1;0;1;2;0;0;2;1;3;3;3;2;3;2;3;0;1;2;0;3;1;3;1;0;3;1;3;3;0;2;3;1;3;3;3;0;2;1;1;3;1;0;1;0;3;2;1;3;3;3;3;3;1;3;3;1;1;0;2;2;3;2;2;2;1;2;0;1;2;1;0;1;1;3;1;1;2;1;3;3;3;3;2;0;3;1;1;0;2;0;0;2;3;2;3;1;3;1;0;3;3;0;0;0;1;2;0;2;1;1;1;2;3;3;0;1;3;0;3;0;2;1;0;1;0;2;1;0;3;1;3;3;3;3;3;3;3;3;2;1;2;0;3;1;2;1;0;2;1;2;3;1;0;3;2;0;0;0;2;3;0;0;3;3;0;1;0;3;3;0;1;1;0;0;0;2;3;1;3;0;3;2;0;1;1;3;1;2;0;2;1;3;3;1;3;2;0;0;0;3;0;0;3;1;0;2;1;2;1;2;2;0;2;1;1;2;2;0;0;3;3;2;0;2;0;3;2;3;2;0;0;2;0;2;0;1;3;1;2;2;3;1;3;2;1;1;0;1;3;1;1;2;0;2;3;3;2;2;3;1;1;0;2;2;1;3;1;0;3;1;0;0;1;1;2;0;0;2;0;2;1;0;2;3;0;3;3;1;2;2;0;3;0;0;3;1;3;2;0;1;2;0;0;1;2;2;3;3;1;2;1;1;1;1;1;2;2;2;2;0;2;1;3;3;2;2;2;2;2;2;1;2;2;0;3;1;1;0;2;3;3;3;3;2;3;1;1;0;0;1;3;1;1;2;0;1;0;1;1;2;0;2;3;2;1;3;3;3;0;1;1;2;2;1;2;0;0;0;3;0;1;0;0;0;1;1;3;0;2;3;3;0;3;3;1;0;3;0;3;0;1;2;1;1;3;2;2;2;0;1;1;0;2;3;0;3;2;3;0;1;3;1;0;3;1;0;1;3;2;1;2;3;0;2;0;3;3;2;2;1;1;3;3;2;3;0;0;2;2;3;3;2;1;1;2;0;1;3;3;2;1;0;0;1;0;3;3;2;3;1;1;1;1;3;3;0;3;1;1;0;1;3;0;1;0;2;2;1;1;1;2;2;1;1;3;3;0;2;3;0;2;2;3;0;3;0;2;0;1;3;2;1;1;0;3;2;2;1;3;2;0;3;3;0;3;0;3;3;1;2;1;3;0;3;3;3;0;3;2;0;1;2;2;2;3;0;0;1;0;3;1;0;2;3;1;2;0;2;3;0;1;0;0;0;0;3;0;1;3;0;0;0;3;2;0;2;3;1;2;3;1;0;3;1;2;1;1;3;2;1;0;1;3;1;0;1;1;2;2;2;0;3;2;3;1;3;3;3;3;1;3;1;3;1;0;2;3;1;0;3;3;1;1;0;2;2;2;0;3;1;1;2;2;1;0;0;3;1;1;1;2;0;3;2;3;2;1;1;2;2;3;3;1;3;0;0;0;0;1;1;2;1;3;3;2;3;0;0;3;0;1;3;3;1;1;0;2;0;0;1;3;2;2;2;0;2;2;1;1;2;1;0;3;0;3;3;1;2;3;2;2;0;2;0;3;1;0;1;1;2;2;2;0;1;2;1;1;2;1;3;0;1;1;0;1;3;1;3;0;3;2;0;1;0;2;2;2;1;1;1;1;1;1;0;0;0;0;1;2;1;2;3;2;3;3;3;2;3;0;1;3;2;1;0;3;0;3;0;2;0;0;1;0;3;2;2;3;2;2;0;1;1;0;3;1;1;3;2;0;1;0;0;0;0;3;0;3;0;0;0;0;3;3;0;0;1;1;3;2;3;2;1;2;2;2;3;0;2;2;0;1;0;1;1;2;0;1;2;3;0;3;1;0;2;0;0;0;1;2;0;1;0;2;3;1;1;3;1;2;2;3;1;0;0;1;1;1;0;2;3;2;0;2;1;3;3;2;2;2;2;1;3;1;1;0;0;0;3;2;2;2;1;2;1;0;3;1;1;2;3;1;1;1;3;0;0;1;1;3;2;3;0;2;2;0;2;1;2;3;2;1;2;0;3;2;3;0;1;0;3;2;1;2;2;1;1;1;0;1;3;2;0;2;0;1;2;2;1;2;3;0;1;0;0;2;0;0;2;1;1;2;0;1;2;0;1;1;2;3;3;0;1;0;1;1;1;0;3;3;2;3;2;3;2;3;0;0;2;0;2;0;1;0;0;2;2;3;1;1;0;2;1;3;2;1;3;2;0;3;0;2;0;1;1;2;3;1;2;3;2;2;0;2;0;3;2;3;0;3;0;0;1;2;3;1;3;0;3;1;1;3;2;2;2;1;0;3;0;2;1;0;3;1;1;2;3;3;0;0;0;3;0;0;0;0;2;0;1;1;0;1;2;3;1;2;2;2;2;0;3;2;3;0;1;1;1;0;1;3;2;0;3;2;1;2;0;0;2;1;0;2;1;3;3;1;1;1;2;3;1;0;2;2;0;3;2;1;1;2;0;0;1;3;1;2;1;3;3;2;0;0;2;3;1;2;2;2;1;3;0;0;2;0;0;2;0;0;3;2;2;2;3;0;3;2;2;2;2;3;0;0;0;3;2;3;1;2;3;3;3;3;1;1;0;1;0;2;0;3;0;3;3;3;3;2;3;3;1;0;3;3;0;1;2;2;0;3;2;0;2;1;3;2;1;0;0;0;2;2;2;0;2;0;0;1;0;3;1;0;1;0;3;2;3;0;0;2;0;2;2;1;0;3;2;3;0;2;1;1;1;3;1;1;1;3;3;1;1;2;0;2;1;3;2;2;0;1;2;0;1;1;3;1;0;3;3;0;3;0;3;3;1;3;3;3;2;2;0;3;0;3;2;0;1;2;3;1;3;0;2;1;1;2;0;1;1;0;1;2;2;0;2;0;0;3;1;1;1;2;2;3;1;3;0;0;1;2;3;1;0;1;2;1;1;0;1;0;1;1;2;0;3;2;3;1;1;2;1;2;0;0;1;1;1;1;2;1;2;2;2;3;3;2;1;3;1;2;0;0;2;0;0;0;3;0;1;1;3;3;3;3;3;0;2;2;3;3;3;2;2;2;3;0;3;2;2;2;3;1;2;1;0;3;1;0;2;0;0;3;1;2;1;2;2;3;3;0;1;1;2;2;2;1;3;0;1;0;2;2;2;2;2;3;1;1;1;0;1;1;2;2;3;3;0;0;3;2;0;2;2;1;3;1;2;1;3;2;3;3;3;0;2;2;1;3;2;1;0;2;2;2;3;0;3;3;1;2;0;0;1;1;2;0;3;0;0;0;1;1;0;1;3;3;1;1;1;0;1;3;0;1;2;2;0;0;3;1;0;3;0;3;3;0;1;3;3;3;1;3;3;3;2;1;1;2;3;1;3;2;3;3;0;1;0;2;3;1;1;3;0;1;2;1;2;3;0;1;0;3;3;3;1;3;2;1;1;0;2;2;1;1;2;1;3;1;1;2;3;1;1;0;1;2;3;3;0;2;1;3;0;1;2;1;2;2;3;0;2;1;1;0;3;2;0;0;2;0;1;2;3;2;1;2;2;3;0;1;0;1;1;1;2;2;1;3;1;1;0;2;2;1;3;3;1;3;3;3;3;1;1;2;1;3;1;1;1;1;1;2;2;2;2;2;2;0;2;0;0;3;0;1;3;2;1;3;2;0;1;1;2;0;1;1;2;0;2;3;0;0;1;1;0;1;1;0;0;3;1;1;2;0;3;2;1;3;2;3;1;3;2;2;0;3;0;3;2;3;2;2;1;1;1;0;1;2;0;3;2;1;3;0;2;2;0;1;2;3;2;1;0;0;0;2;0;1;0;2;2;2;1;0;0;3;2;0;2;0;2;3;1;0;2;0;3;0;3;3;3;0;3;2;1;0;2;2;3;0;1;0;3;2;0;0;3;1;2;3;1;2;0;0;2;3;2;1;1;3;1;2;1;1;2;3;3;1;3;2;1;1;2;1;3;3;0;2;2;2;0;3;0;2;1;3;0;2;3;3;1;0;3;2;3;3;1;0;3;2;2;0;3;1;0;3;1;2;1;1;1;2;2;1;3;0;1;0;3;3;1;0;1;0;3;0;1;2;1;1;3;3;1;1;3;1;0;3;2;3;2;0;1;1;3;1;0;1;0;2;2;2;2;0;1;0;2;1;2;1;1;0;2;2;1;3;3;3;3;0;2;0;1;2;3;1;3;2;0;2;3;3;1;1;1;0;0;0;3;2;0;2;2;0;3;0;2;0;2;1;2;2;0;0;3;3;2;3;0;2;2;1;0;2;2;3;3;3;3;3;0;1;1;1;2;3;2;1;2;2;0;2;0;1;1;0;0;2;2;0;3;3;3;1;1;2;3;0;2;1;2;2;1;0;2;0;2;1;0;3;0;2;3;3;3;2;3;0;0;3;3;0;1;3;0;3;2;3;1;0;0;0;2;0;1;2;0;1;3;2;2;3;0;3;2;2;2;1;2;1;1;1;2;1;0;2;1;2;2;2;3;2;0;3;0;2;3;2;2;1;2;3;0;1;2;3;0;2;3;1;3;2;1;1;3;3;1;1;0;3;1;2;0;1;0;2;0;2;0;2;2;3;3;1;2;3;2;2;0;0;1;1;3;3;1;3;1;0;2;1;1;1;2;1;0;0;0;0;2;3;0;0;2;3;3;1;0;1;3;1;1;2;3;2;3;2;3;0;1;2;0;2;1;0;2;1;0;2;1;2;1;1;1;1;2;1;3;3;2;0;3;2;3;3;1;1;2;0;3;3;3;2;2;2;3;2;3;3;1;2;0;2;0;2;3;3;2;1;1;1;2;1;1;0;0;1;2;1;3;2;1;3;0;2;0;2;3;1;3;0;1;1;0;1;3;2;0;1;1;0;1;1;2;0;3;1;3;2;0;2;0;1;2;3;1;0;0;3;3;0;1;2;3;2;2;1;1;3;3;3;2;0;2;1;0;2;0;2;0;2;1;1;1;1;1;1;2;1;2;1;2;3;0;2;2;1;2;1;2;0;2;1;3;3;2;3;0;3;0;2;1;2;2;2;2;3;0;0;3;3;1;0;0;1;3;0;3;1;1;0;0;0;3;2;2;1;3;2;3;3;1;3;3;3;3;2;2;0;3;0;0;2;2;0;2;1;1;3;1;2;1;2;2;0;3;0;2;0;2;1;3;3;1;3;0;3;2;2;2;0;3;0;3;0;3;3;2;3;3;1;2;1;2;3;0;0;1;1;2;1;2;2;1;3;1;2;2;1;0;3;3;3;3;0;3;2;0;3;1;3;0;0;1;1;1;0;2;2;2;0;1;2;0;2;1;2;0;1;1;0;2;2;0;2;3;1;1;1;0;0;0;2;0;1;3;3;3;1;1;2;1;2;3;3;0;1;1;0;2;1;0;1;3;1;3;0;1;2;2;2;1;2;2;0;1;0;1;0;1;2;0;2;2;1;3;3;1;1;2;0;0;0;0;2;1;2;1;2;1;1;2;1;1;3;1;1;3;3;1;0;2;2;0;2;3;1;1;1;2;1;1;0;3;1;0;2;0;1;3;1;1;1;3;1;3;3;0;2;0;1;1;3;3;3;0;0;1;2;3;1;3;0;0;1;1;0;2;3;2;0;3;1;3;2;0;0;1;1;3;2;0;2;0;1;1;2;3;1;3;3;2;0;2;1;0;2;2;0;2;3;1;2;3;1;1;3;2;0;3;3;2;3;2;2;1;0;1;1;2;2;3;2;1;3;3;2;3;1;3;0;2;2;1;0;1;1;2;0;1;0;2;0;0;3;2;1;1;0;3;3;3;1;2;3;1;2;0;2;0;3;1;0;2;1;1;2;3;0;0;0;1;1;0;0;2;3;3;2;3;3;3;0;1;2;0;1;2;0;2;0;1;0;1;0;3;1;1;3;0;2;0;2;3;3;1;3;0;0;3;0;1;0;2;2;2;1;2;1;2;1;2;3;0;3;3;0;0;0;1;2;3;2;1;2;2;0;2;2;2;2;2;2;2;3;3;2;3;3;1;2;3;0;1;0;0;3;0;3;0;1;3;3;0;2;0;1;2;2;1;2;0;1;0;3;3;1;0;0;0;1;0;0;2;3;2;1;0;1;3;1;2;1;3;0;2;1;3;2;1;0;1;2;3;2;0;2;2;0;0;0;1;3;0;1;0;2;1;0;2;0;0;0;1;2;2;1;2;1;3;0;0;0;1;0;3;2;1;3;1;1;0;1;2;3;2;0;2;0;1;3;1;1;2;3;2;0;0;2;0;1;1;1;3;3;1;0;3;0;0;1;2;1;2;0;1;0;0;3;0;2;0;0;3;1;1;1;0;2;3;1;3;0;0;2;1;0;3;2;2;3;3;0;2;0;1;1;0;3;0;1;3;3;3;3;0;1;3;2;3;2;3;3;0;0;2;1;2;1;0;3;3;0;0;0;0;0;3;2;2;1;3;2;1;3;0;0;2;3;1;1;2;0;2;3;3;3;1;3;0;0;2;0;3;1;3;0;2;2;1;0;0;0;0;2;0;2;2;1;2;1;1;3;1;2;3;1;0;0;0;3;1;0;1;1;2;2;3;0;1;1;2;3;1;0;1;2;0;0;3;1;0;3;1]
wi=[1;1;4;8;7;1;3;4;3;5;4;5;7;8;4;2;5;6;8;7;4;6;3;7;2;6;2;4;1;8;2;7;5;6;6;5;6;3;8;3;6;8;4;1;3;5;1;7;3;6;8;2;2;1;4;6;3;2;8;6;2;4;6;8;5;6;2;8;2;2;3;5;3;8;1;4;2;7;5;6;8;2;5;6;6;6;5;8;3;8;6;8;7;5;5;7;6;1;8;1;5;8;8;6;6;1;2;6;5;3;1;4;5;5;1;8;5;1;5;2;2;3;2;7;6;6;3;1;8;3;4;2;1;6;3;1;4;1;7;1;7;3;5;3;6;2;6;3;1;2;4;4;5;7;4;8;8;5;6;4;5;8;2;8;7;5;7;6;3;3;6;7;3;1;4;3;3;1;4;6;6;1;1;1;4;4;5;6;2;3;5;1;7;1;7;1;1;8;3;2;5;4;5;1;3;5;5;4;3;8;4;4;7;3;5;6;5;2;2;4;2;5;5;7;8;4;6;5;7;8;8;2;8;8;3;6;8;7;4;1;4;1;3;7;1;4;1;4;7;6;5;1;7;3;3;2;1;3;4;2;2;1;8;5;7;7;5;8;5;6;8;5;3;3;2;8;4;6;3;1;1;6;3;4;7;5;6;7;8;1;8;7;1;3;5;8;2;3;1;4;7;8;5;6;7;5;7;1;8;7;6;6;6;4;2;8;8;8;8;8;6;7;4;1;8;7;8;7;5;4;3;3;6;8;3;1;2;1;2;1;1;7;7;4;7;7;5;8;2;6;5;6;5;8;7;5;4;8;7;5;3;1;2;7;1;5;3;6;7;7;3;6;2;3;6;3;5;2;8;7;3;6;6;7;4;4;3;2;7;3;3;5;8;2;4;7;7;4;2;2;3;1;7;8;3;4;7;3;6;2;7;8;8;2;4;7;8;3;5;5;7;6;2;4;3;7;8;3;1;7;4;7;4;8;5;2;7;8;3;2;5;8;6;3;4;8;5;4;2;2;7;7;7;1;7;1;3;2;5;8;1;3;4;6;8;5;4;3;3;1;3;5;7;7;1;6;6;2;2;5;4;1;7;3;2;3;6;3;1;4;3;1;2;5;2;2;6;7;5;1;2;4;5;7;3;1;3;2;1;6;3;6;7;1;5;3;3;6;5;3;5;8;1;3;6;5;3;4;8;4;2;4;7;4;7;6;1;4;5;7;2;1;7;8;4;2;3;1;2;4;4;3;5;7;3;7;1;7;2;3;2;2;1;3;6;8;3;3;5;4;1;1;1;8;4;7;3;8;7;5;2;3;5;5;6;8;5;6;6;6;7;6;8;2;3;3;5;6;5;6;7;2;3;6;4;4;5;5;8;5;6;2;6;7;8;1;8;4;3;2;8;8;5;8;1;3;7;2;8;5;3;3;8;5;4;6;7;5;4;4;2;3;1;3;5;1;1;7;8;3;2;4;1;3;1;1;4;1;4;5;6;6;6;3;7;7;4;8;6;2;5;4;1;1;5;8;5;2;8;7;7;4;6;7;6;2;1;4;2;2;2;8;3;5;5;1;2;6;2;6;2;4;4;5;7;1;1;1;8;8;8;6;2;3;2;3;7;6;3;5;7;6;4;2;4;8;8;8;7;6;7;6;6;3;4;1;5;3;5;3;1;3;8;5;8;6;1;1;6;5;6;1;8;6;4;1;6;3;7;3;2;3;4;1;4;3;8;1;6;4;3;8;4;2;4;4;5;6;8;4;2;1;7;3;3;3;7;5;8;6;8;7;8;8;4;6;1;2;1;8;7;8;8;7;5;1;5;4;7;4;1;5;4;1;6;7;2;4;5;3;6;8;3;7;6;1;3;1;3;7;7;8;8;3;6;6;8;8;8;5;5;5;7;4;1;3;2;8;7;2;4;3;6;4;5;8;4;3;7;1;2;8;4;1;7;3;3;1;1;4;6;6;2;7;4;1;3;3;5;8;6;4;5;4;3;8;4;3;2;7;1;2;5;3;6;7;1;2;7;6;8;7;5;1;2;1;6;2;3;4;1;6;3;2;8;4;4;7;2;4;1;6;7;6;3;3;3;3;7;5;2;4;3;7;5;1;7;6;6;4;5;6;6;2;4;4;4;6;1;2;3;8;8;2;4;1;6;2;5;4;3;2;1;5;7;6;5;5;2;3;8;7;3;3;7;8;1;6;5;7;3;8;7;1;8;1;8;4;6;6;6;1;2;8;1;7;8;3;4;1;4;7;7;3;6;6;3;4;2;3;3;2;5;6;6;2;8;6;2;6;6;4;3;5;5;3;7;4;7;5;8;2;3;4;7;6;2;5;8;2;3;3;3;3;5;4;2;7;5;1;1;7;2;5;6;2;4;5;5;1;5;3;7;2;4;3;2;7;2;6;3;3;3;4;6;1;2;5;4;1;2;6;3;4;2;2;2;2;6;1;5;7;7;1;5;8;4;3;7;6;3;5;3;6;2;2;7;1;8;4;6;8;5;5;6;6;7;4;3;7;6;5;6;5;4;7;1;1;3;1;4;5;7;6;5;5;8;8;7;7;5;2;2;1;2;6;8;8;4;1;5;3;6;1;2;5;8;6;4;5;2;2;4;8;8;5;5;2;6;4;7;7;5;6;1;5;6;4;8;8;4;8;5;4;6;4;7;4;8;4;7;7;5;5;3;2;7;5;1;8;1;5;2;2;6;7;7;8;1;6;3;2;7;3;3;6;5;1;3;5;3;1;5;7;1;1;6;7;3;1;7;7;1;5;3;6;3;2;3;4;3;2;1;5;8;7;2;6;1;8;1;8;4;1;8;3;5;7;5;5;7;7;3;3;5;3;4;1;4;3;8;4;4;4;2;5;1;5;7;3;1;1;3;2;7;1;2;6;4;5;1;4;1;2;3;5;6;7;4;6;6;4;5;7;1;7;2;3;4;2;1;8;3;8;4;6;6;4;1;7;8;1;2;6;6;2;4;4;1;5;2;7;5;3;1;4;7;1;2;1;6;4;3;5;7;2;5;3;3;6;6;6;1;8;8;5;2;3;7;2;1;7;5;7;2;2;4;3;8;1;8;5;3;5;4;2;8;5;7;2;8;4;1;6;2;2;1;6;2;1;7;5;7;5;2;3;2;2;5;6;5;6;7;2;1;7;7;7;2;7;2;2;6;4;7;1;8;3;4;5;7;7;1;6;4;7;2;6;7;8;6;2;6;3;1;7;8;3;2;3;8;5;2;5;1;2;5;4;5;4;5;6;7;4;7;5;1;8;6;5;5;4;1;7;3;1;6;3;8;1;8;1;2;6;2;6;2;7;4;2;8;4;1;3;6;5;3;3;5;4;2;4;4;5;6;7;8;7;1;1;5;8;7;2;5;1;7;1;2;7;5;6;6;5;8;6;5;7;1;4;4;2;6;6;2;2;5;1;6;4;1;4;1;2;2;3;6;6;6;4;5;2;5;7;8;8;4;8;3;8;2;5;3;8;5;3;7;5;3;1;7;6;8;5;7;2;8;6;7;4;2;8;4;6;6;8;8;8;4;6;7;7;5;4;8;1;2;5;1;8;7;2;7;3;2;8;2;6;8;6;2;2;7;3;7;2;3;3;1;6;4;8;7;4;5;4;3;8;8;1;3;1;2;5;4;4;1;4;1;5;4;4;8;6;8;5;5;6;4;8;7;3;2;2;7;6;5;1;7;6;3;3;8;4;3;6;4;5;1;4;8;7;1;6;3;7;4;4;7;3;7;1;3;7;3;6;5;4;4;4;6;1;7;1;3;8;6;8;3;3;3;2;1;6;2;7;6;4;2;4;5;4;2;1;4;7;6;6;2;6;6;6;3;8;1;8;4;5;2;5;2;3;5;1;1;3;4;8;5;3;5;4;7;6;5;7;4;7;6;5;7;7;5;4;7;1;1;6;4;1;3;1;3;8;1;7;7;7;3;7;4;1;6;4;1;2;5;5;6;4;6;7;7;8;2;4;3;2;6;5;4;5;6;5;1;6;3;4;6;5;7;2;1;5;5;5;1;2;4;5;2;5;5;3;4;8;7;4;2;1;8;8;3;3;1;6;3;8;1;2;2;6;2;2;3;3;7;2;7;1;8;3;7;6;4;6;5;4;2;8;2;8;3;7;2;3;6;5;2;5;3;2;1;7;6;5;3;4;2;1;1;7;6;3;5;1;8;7;5;8;4;6;1;7;6;8;7;2;7;5;5;2;4;2;6;5;6;5;6;2;7;3;6;5;4;8;1;6;1;2;1;6;8;7;6;4;2;2;4;7;5;8;7;4;4;1;5;8;1;1;5;3;7;5;4;1;5;1;4;3;6;8;5;8;4;4;4;8;2;6;4;3;5;2;7;2;6;6;1;3;8;1;7;8;8;7;5;6;3;3;4;4;7;7;6;4;5;4;2;1;8;1;8;6;6;1;1;4;7;6;7;6;8;7;8;1;6;3;2;5;6;3;8;3;6;5;6;7;3;5;8;6;2;5;6;4;1;5;5;5;2;4;2;8;4;2;4;5;8;3;8;3;6;7;5;1;3;8;3;2;7;7;6;7;5;5;4;3;6;5;4;7;5;4;3;7;8;3;3;2;2;8;2;2;4;6;8;6;4;7;8;6;6;6;2;3;7;1;8;6;5;3;4;6;5;1;4;3;7;8;6;7;3;2;1;3;4;2;3;8;5;4;4;5;2;7;6;7;8;4;8;2;1;4;7;2;4;3;6;2;6;8;7;7;5;6;5;8;5;4;6;1;8;2;4;3;7;4;4;3;8;7;5;3;4;7;6;1;4;2;5;1;3;7;7;7;2;4;5;3;1;3;2;1;8;3;5;4;5;1;8;2;7;6;4;8;2;4;4;5;8;2;3;1;3;1;4;4;3;1;1;3;4;8;4;8;3;3;4;2;6;4;7;2;2;8;4;2;2;7;8;2;8;3;1;6;8;4;1;7;7;5;2;4;4;2;3;3;6;6;2;7;2;8;3;6;8;7;6;8;2;1;6;1;1;1;2;5;8;3;1;6;3;7;6;2;6;8;2;2;5;2;5;4;8;5;1;8;5;8;3;8;6;2;8;8;2;8;7;5;4;8;5;7;7;3;1;4;8;1;4;7;3;8;5;4;3;7;3;5;6;6;1;1;8;6;1;7;8;6;3;6;3;8;7;5;2;8;6;3;4;3;6;5;6;7;1;2;2;8;4;8;4;2;8;7;7;7;6;4;6;1;1;8;4;6;7;3;1;1;6;8;2;5;3;5;3;1;8;7;5;6;3;8;8;4;3;7;6;6;4;8;3;5;5;2;5;2;7;6;2;4;7;2;7;7;8;5;4;7;1;1;3;2;6;6;3;3;5;6;2;2;7;6;2;1;7;7;4;4;4;6;5;4;7;8;5;4;6;4;4;7;8;7;6;7;2;1;6;4;8;6;2;7;7;2;2;2;8;5;4;7;8;6;1;2;8;1;6;3;6;4;3;5;8;5;6;4;6;6;4;3;8;6;8;3;8;3;2;5;3;6;4;8;8;1;8;3;1;8;5;8;1;2;2;8;7;3;6;5;2;3;2;6;3;8;5;6;8;3;4;1;8;2;7;1;7;7;1;6;1;5;1;8;2;3;4;3;1;4;8;7;8;2;3;1;5;5;4;8;8;6;8;6;8;6;5;2;2;1;1;7;4;2;7;3;5;5;1;5;2;4;1;7;5;5;4;2;4;7;3;4;7;5;6;2;6;2;6;5;6;7;6;7;5;7;1;2;1;6;8;5;1;8;3;5;2;7;1;4;2;3;7;1;6;5;8;8;1;3;8;7;2;5;8;7;7;1;1;1;7;4;5;3;3;2;3;5;5;4;5;8;1;6;2;2;5;6;5;7;6;5;2;3;6;5;3;4;8;6;3;3;7;4;3;2;6;2;3;8;1;6;3;2;6;7;5;6;2;7;3;3;1;2;4;3;8;4;6;7;4;2;3;6;2;1;1;8;3;5;5;5;6;2;3;8;1;1;1;6;6;1;4;5;2;1;7;3;4;4;7;6;7;3;5;4;4;3;4;5;3;1;8;8;3;7;6;3;8;7;8;8;5;2;5;2;8;4;2;6;4;5;5;6;6;1;4;3;3;2;4;3;6;8;1;1;3;2;3;8;1;2;3;8;5;4;2;7;4;1;8;4;2;8;2;2;2;5;8;4;1;1;4;8;3;7;2;6;1;4;1;3;5;5;1;4;2;4;6;3;1;5;4;7;5;7;5;5;1;2;8;5;7;1;4;8;7;3;5;1;1;2;4;5;1;8;4;5;5;5;8;1;5;8;7;2;2;3;2;6;2;1;8;8;1;6;1;1;6;5;4;2;1;2;1;4;4;3;5;8;8;1;2;2;2;2;7;3;4;6;1;4;4;4;5;8;1;3;4;4;1;3;5;6;7;2;4;6;3;8;1;7;5;8;7;3;7;8;4;7;7;5;2;4;4;1;2;8;5;3;2;1;6;3;6;8;4;3;4;3;6;7;3;1;3;8;5;6;8;6;2;2;5;5;1;3;3;1;8;8;8;2;3;3;8;1;8;3;7;3;1;5;6;2;5;6;5;5;2;5;8;3;5;7;3;3;2;3;6;4;3;2;5;1;8;3;8;5;2;5;7;6;4;3;3;5;7;1;5;6;3;5;1;2;3;3;4;2;5;2;5;2;5;8;4;7;6;8;8;8;1;4;5;3;8;3;4;6;3;5;6;8;1;6;6;8;1;8;5;6;7;3;2;1;7;8;7;1;4;6;7;7;6;1;7;6;1;5;3;3;2;5;8;7;3;8;4;4;5;4;6;7;7;4;3;2;7;2;4;3;5;7;2;1;1;8;1;3;5;5;2;4;4;3;2;8;1;4;3;6;8;5;6;1;2;7;5;4;4;4;3;1;8;6;3;8;5;5;5;2;2;6;7;6;3;6;1;5;6;8;8;2;1;3;1;1;7;3;5;3;7;8;4;3;4;7;5;8;2;1;4;5;7;8;4;4;1;2;7;7;1;6;4;6;1;7;5;7;8;7;2;6;2;1;7;6;7;1;1;7;4;2;2;5;7;8;2;5;8;7;5;8;3;5;2;5;5;6;1;7;2;4;8;2;5;5;3;3;4;4;2;1;7;6;6;2;6;8;3;4;5;1;3;4;5;2;7;5;5;1;3;8;4;5;6;2;4;4;8;8;5;8;3;5;3;4;5;4;7;7;8;5;6;6;7;5;6;3;7;7;8;1;7;7;4;7;8;4;2;7;5;7;8;4;6;6;1;1;6;1;4;2;6;2;3;7;5;3;6;3;4;4;6;1;8;4;6;3;4;7;7;4;2;2;8;6;5;6;7;8;5;2;2;7;8;1;6;7;4;7;1;7;3;5;4;7;4;1;6;5;5;1;8;5;2;5;7;3;7;3;7;2;5;6;2;4;6;2;3;7;6;5;1;6;3;1;3;2;6;6;6;6;8;2;4;2;2;4;1;3;4;3;6;8;6;3;5;5;2;5;8;1;2;8;2;3;3;7;7;2;8;6;3;2;1;2;5;3;5;8;2;8;5;3;6;1;3;4;6;8;8;7;5;1;6;2;1;5;8;2;4;8;8;7;3;2;1;6;1;3;8;3;2;8;3;1;6;5;1;3;8;3;6;2;6;5;2;4;5;4;2;6;2;3;6;5;4;6;4;5;7;2;7;1;8;2;6;1;3;6;6;1;5;6;2;4;1;7;6;1;6;4;5;2;8;4;1;7;2;2;6;6;4;4;2;5;2;4;2;1;1;2;1;5;7;3;6;8;8;3;3;4;4;5;7;2;8;2;5;3;6;2;5;2;4;4;6;1;6;7;4;5;8;3;7;6;7;5;1;1;5;7;7;7;6;7;6;3;7;1;3;5;5;7;7;5;8;2;4;3;5;2;8;5;2;3;1;3;2;2;8;6;4;3;2;6;1;7;5;8;6;7;3;5;2;4;6;3;3;5;3;1;6;1;4;7;5;2;5;7;3;2;7;5;7;6;4;7;3;3;6;2;6;7;7;5;7;4;3;3;1;6;3;6;5;8;2;8;2;3;3;6;7;3;2;1;6;7;1;1;6;1;3;5;1;8;3;5;7;3;3;7;7;5;2;1;5;2;3;8;7;2;1;8;7;5;3;6;6;7;1;6;2;7;2;6;2;3;1;4;3;5;4;8;1;7;4;1;8;5;6;3;1;1;3;6;7;5;3;1;8;3;5;7;6;6;8;6;3;1;4;4;7;2;7;2;1;6;2;5;1;7;3;2;3;6;2;8;2;1;4;6;2;1;3;4;2;8;6;3;7;4;6;2;2;3;1;1;8;4;2;6;7;3;6;3;7;2;1;7;7;6;3;8;4;2;8;7;3;5;5;8;8;4;1;3;3;8;4;7;2;5;5;3;1;2;8;3;4;5;1;8;5;6;6;4;4;6;7;7;7;5;4;2;3;3;2;7;1;1;4;6;1;7;1;3;7;3;1;3;2;1;5;5;2;4;6;3;6;8;3;1;5;1;8;2;3;7;3;4;2;5;6;5;5;4;2;8;8;6;8;4;8;3;7;4;6;2;4;7;3;6;1;2;5;2;5;7;2;6;2;8;8;2;6;7;3;1;1;4;1;3;2;3;5;8;6;5;8;2;8;6;8;4;1;1;8;7;5;5;3;8;1;8;4;2;8;1;4;6;2;5;8;4;6;8;8;5;8;1;1;4;6;8;2;7;4;6;1;1;3;2;7;3;6;6;4;3;8;4;1;2;6;7;2;4;2;3;7;5;1;4;4;8;5;4;4;8;1;4;5;6;8;8;8;1;8;1;6;5;2;3;3;3;8;2;3;1;2;3;8;3;4;8;8;6;6;5;3;2;6;3;3;8;4;8;6;2;5;8;6;1;5;5;3;6;8;1;7;3;3;6;5;2;6;4;7;6;3;7;6;8;4;4;1;1;1;6;6;4;3;2;1;1;1;3;7;3;5;7;7;4;8;6;8;2;2;2;3;1;6;6;5;1;2;4;6;2;1;1;1;1;1;8;5;3;8;6;5;6;8;2;3;3;1;5;7;7;8;7;7;8;7;6;8;8;8;2;7;2;3;6;7;1;8;7;7;2;5;3;1;8;5;4;6;7;5;1;2;7;8;8;1;7;7;2;8;8;6;5;5;1;1;3;1;6]

zi(zi==0)=4;      %加减法互换
zi(zi==1)=0;
zi(zi==4)=1;

for i=1:r_block
    Q4=DNA_bian(fenkuai(t,Q5,i),wi(i));

    Q3=DNA_bian(fenkuai(t,bian_R,i),yi(i));
    
    Q2=DNA_yunsuan(Q4,Q3,zi(i));
   
    
    xx=floor(i/block)+1;    %（floor朝负无穷方向取整）
    yy=mod(i,block);        %（mod求余数）
    if yy==0
        xx=xx-1;
        yy=block;
    end
   
    Y((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q2,xi(i));
end

%% 6、去除加密时补的零
M1=0;   %加密时补零的参数，M1=mod(M,t);作为密钥
N1=0;   %加密时补零的参数，N1=mod(N,t);作为密钥
if M1~=0
    Y=Y(1:M-t+M1,:,:);
end
if N1~=0
    Y=Y(:,1:N-t+N1,:);
end

[m,n]=size(Y);
max_1=221365.353693676;
min_1=-118683.475787179;

for i=1:m
    for j=1:n    
       Y(i,j)= Y(i,j)*(max_1-min_1)/256+min_1;
    end
end

Y1=Y(1:m/2,:);
Y2=Y(m/2+1:m,:);

[lorenze_x,lorenze_y,lorenze_z]=Lorenz_chaotic_system(y_02,y_03,y_04,m*n);
lorenze_x=lorenze_x(3002:length(lorenze_x));
lorenze_y=lorenze_y(3002:length(lorenze_y)); 
lorenze_z=lorenze_z(3002:length(lorenze_z)); 
 
lorenze_z=mod(round(lorenze_z*10^4),2);
  measure=zeros((m/2)*n,1);
for Q=1:((m/2)*n) 
    if(lorenze_z==0)
        measure(Q)=lorenze_x(Q);
    else
        measure(Q)=lorenze_y(Q);
    end
end
measure=reshape(measure,(m/2),n);

X3=zeros(m,n);  %  恢复矩阵
for i=1:n  %  列循环       
    rec=omp(Y1(:,i),measure,m);
    X3(:,i)=rec;
end

X4=zeros(m,n);  %  恢复矩阵
for i=1:n  %  列循环       
    rec=omp(Y2(:,i),measure,m);
    X4(:,i)=rec;
end


ww=DWT(m);
[c,d]=size(ww);

u=3.9999;     %Logistic参数μ，自定为3.99
logid=3;
logix=zeros(1,logid*c+1000);        %预分配内存
logiy=zeros(1,logid*d+1000);
logix(1)=y_01;
logiy(1)=(y_01+y_02)/2;

for i=1:logid*c+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    logix(i+1)=u*logix(i)*(1-logix(i));
end
for i=1:logid*d+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    logiy(i+1)=u*logiy(i)*(1-logiy(i));
end
logix=logix(1001:logid:length(logix));            %去除前1000点，获得更好的随机性（length数组长度）
logiy=logiy(1001:logid:length(logiy));  

[~,Ux]=sort(logix,'descend');      %（对矩阵进行降序排列，Ux为排序后元素在原矩阵的列位置）
[~,Uy]=sort(logiy,'descend');  

for i=1:c   %行置换
    temp = ww(i,:);
    ww(i,:) = ww(Ux(i),:);
    ww(Ux(i),:) = temp;    
end
for i=1:d   %列置换
    temp = ww(i,:);
    ww(i,:) = ww(Uy(i),:);
    ww(Uy(i),:) = temp;    
end
    
X2=ww'*sparse(X3)*ww;  %  小波反变换
X2=full(X2); 
X1=ww'*sparse(X4)*ww;  %  小波反变换
X1=full(X4); 
figure(1);
imshow(uint8(X1));
figure(2);
imshow(uint8(X2));
figure(3);



%  OMP的函数
%  s-测量；T-观测矩阵；N-向量大小
function hat_y=omp(s,T,N)

Size=size(T);                                     %  观测矩阵大小
M=Size(1);                                        %  测量
hat_y=zeros(1,N);                                 %  待重构的谱域(变换域)向量                     
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_n=s;                                            %  残差值

for times=1:M/4;                                  %  迭代次数(稀疏度是测量的1/4)
    for col=1:N;                                  %  恢复矩阵的所有列向量
        product(col)=abs(T(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置
    Aug_t=[Aug_t,T(:,pos)];                       %  矩阵扩充
    T(:,pos)=zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零）
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  最小二乘,使残差最小
    r_n=s-Aug_t*aug_y;                            %  残差
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
    
    if (norm(r_n)<9)                              %  残差足够小  2范数
        break;
    end
end
hat_y(pos_array)=aug_y;                           %  重构的向量




