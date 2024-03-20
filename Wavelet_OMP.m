%  本程序实现图像LENA的压缩传感
%  程序作者：沙威，香港大学电气电子工程学系，wsha@eee.hku.hk
%  算法采用正交匹配法，参考文献 Joel A. Tropp and Anna C. Gilbert 
%  Signal Recovery From Random Measurements Via Orthogonal Matching
%  Pursuit，IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 53, NO. 12,
%  DECEMBER 2007.
%  该程序没有经过任何优化

function Wavelet_OMP

clc;clear

%  读文件
X=imread('Lena512.bmp');
% X=X(:,:,1);        %R
%imwrite(X,'../Wavelet_OMP/walkbridge.bmp','bmp');    %将图片以bmp形式保存  

% I=ind2gray(X,map);
% figure(2);
% imshow(uint8(X));

%X=imread('4.2.06.tiff');
%X=imread('Rice.tif');
%X= rgb2gray(X);



X=double(X);
[a,b]=size(X);
%X=ones(a,b);

%  小波变换矩阵生成
ww=DWT(a);
% wavename='db2';
% [cA2,cB2,cD2,cE2]=dwt2(X,wavename);


%  小波变换让图像稀疏化（注意该步骤会耗费时间，但是会增大稀疏度）
X1=ww*sparse(X);   %   sparse将矩阵x转换为稀疏矩阵，矩阵中去除零元素
X1=full(X1);           %   full 把稀疏矩阵转换成一个全矩阵

%  随机矩阵生成
%  M=256;
%  R=randn(M,a);    % randn生成M*a的标准正态分布的随机矩阵 测量矩阵
%  [lorenze_x,lorenze_y,lorenze_z]=Lorenz_chaotic_system(0.1234,0.2345,0.3456,a*b);
%  lorenze_x=lorenze_x(3002:length(lorenze_x));
%  lorenze_y=lorenze_y(3002:length(lorenze_y)); 
%  lorenze_z=lorenze_z(3002:length(lorenze_z));
%  [max_lorenze,~]=max(lorenze_x);
%  [max_lorenze,~]=max(lorenze_x);
%  [min_lorenze,~]=min(lorenze_x);
%  [min_lorenze,~]=min(lorenze_x);
%  for i=1:length(lorenze_x)
%     lorenze_a=(lorenze_x-min_lorenze)/(max_lorenze-min_lorenze);
%  end
 Y=zeros(a/2,b);
 for i=1:b
    x0=0.2+0.0001*i;y0=0.2; 
   [heon_x,heon_y]=henon(x0,y0,(a/2)*512-1);
   R=0.08838834756*reshape(heon_y,a/2,b);
   Y(:,i)=R*X1(:,i);
 end

 
% u=3.9999;     %Logistic参数μ，自定为3.99
% logid=3;
% c=a;
% logix=zeros(1,logid*c+1000);        %预分配内存
% 
% logix(1)=0.234;
% for i=1:logid*c+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
%     logix(i+1)=u*logix(i)*(1-logix(i));
% end
% logix=logix(1001:logid:length(logix));            %去除前1000点，获得更好的随机性（length数组长度）
% 
% logiy(1)=0.345;
% for i=1:logid*c+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
%     logiy(i+1)=u*logiy(i)*(1-logiy(i));
% end
% logiy=logiy(1001:logid:length(logiy));            %去除前1000点，获得更好的随机性（length数组长度）
% 
% [~,Ux]=sort(logix,'descend');      %（对矩阵进行降序排列，Ux为排序后元素在原矩阵的列位置）
% [~,Uy]=sort(logiy,'descend');  
% 
% for i=1:a   %行置换
%     temp = R(i,:);
%     R(i,:) = R(Ux(i),:);
%     R(Ux(i),:) = temp;    
% end
% for i=1:b   %列置换
%     temp = R(i,:);
%     R(i,:) = R(Uy(i),:);
%     R(Uy(i),:) = temp;    
% end

%R=reshape(logix,a/2,b);


% [max_1,~]=max(Y);
% [max_1,~]=max(max_1);
% [min_1,~]=min(Y);
% [min_1,~]=min(min_1);
% for i=1:(a/2)
%     for j=1:b
%         Y(i,j)=uint8((Y(i,j)-min_1)/(max_1-min_1)*255);
%      %  W(i,j)= Q(i,j)*(max_1-min_1)/256+min_1;
%     end
% end
% for i=1:(a/2)
%     for j=1:b
%         
%        W(i,j)= Y(i,j)*(max_1-min_1)/256+min_1;
%     end
% end
% 
% 
% %  测量
% Y=R*X1;

%  OMP算法
X2=zeros(a,b);  %  恢复矩阵
for i=1:b  %  列循环   
    x0=0.2+0.0001*i;y0=0.2; 
    [heon_x,heon_y]=henon(x0,y0,(512/2)*512-1);
    R=0.08838834756*reshape(heon_y,a/2,b);
   rec=SL0(R,Y(:,i), 0.004);
    % rec=omp(Y(:,i),R,a);
    X2(:,i)=rec;
end


%  原始图像
figure(1);
imshow(uint8(X));
title('原始图像');

%  变换图像
figure(2);
imshow(uint8(X1));
title('小波变换后的图像');

%  压缩传感恢复的图像
figure(3);
X3=ww'*sparse(X2);  %  小波反变换
X3=full(X3);
subplot(2,2,2);imshow(uint8(X3));title('恢复的图像');

subplot(2,2,1);imshow(uint8(X));title('原始的图像');

%  误差(PSNR)
errorx=sum (sum(abs(uint8(X3)-uint8(X)).^2))        %  MSE误差
psnr=10*log10(255*255/(errorx/a/b))   %  PSNR
 

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
    
    if (norm(r_n)<2)                              %  残差足够小  2范数
        break;
    end
end
hat_y(pos_array)=aug_y;                           %  重构的向量
