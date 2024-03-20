function main

clc;clear
tic
%  读文件
% X1=imread('X111.bmp');
% X2=imread('X444.bmp');


X1=imread('DICOM1.png');
X2=imread('DICOM2.png');

%X2=imread('boat.512.tiff');
% X2=imread('jetplane.tif');
% X2=X2(:,:,1);
figure(1);
imshow(uint8(X1));

%X2=imread('4.tiff');
%X= rgb2gray(X2);
%imwrite(X,'../Wavelet_OMP/Plane.bmp','bmp');    %将图片以bmp形式保存    


X1=double(X1);
X2=double(X2);
[m,n]=size(X1);
X12=[X1;X2];

%% sha256
%SUM_P=sum(X1(:))+sum(X2(:));
%SUM_P=num2str(SUM_P);
%OUTPUT = sha256( SUM_P ) ; 

OUTPUT=Hash(X12,'SHA-256');
% a=hex2dec('0d5c811a4a4b7b77f1a65cd6deb6c6149a1ff4dd2e27504d590ba4afa0168489');
% b=dec2bin(a);

digest1=OUTPUT(1:8);
digest2=OUTPUT(9:16); 
digest3=OUTPUT(17:24); 
digest4=OUTPUT(25:32); 
digest5=OUTPUT(33:40);
digest6=OUTPUT(41:48); 
digest7=OUTPUT(49:56); 
digest8=OUTPUT(57:64); 

for i=1:8   
        a(i)=digest1(i);
        XXXX(i)=  cellstr(a(i));
        b(i)=hex2dec(XXXX(i));
end

y_01=(bitxor((bitxor((bitxor(b(1),b(8))),b(2))),b(7))+bitxor((bitxor((bitxor(b(3),b(6))),b(4))),b(5)))/256;
% y_01=0.0313;
%y_01=0.2;
for i=1:8    
a(i)=digest2(i);
XXXX(i)=  cellstr(a(i));
b(i)=hex2dec(XXXX(i));
end

y_02=(bitxor((bitxor((bitxor(b(1),b(8))),b(2))),b(7))+bitxor((bitxor((bitxor(b(3),b(6))),b(4))),b(5)))/256;
%y_02=0.0898;
for i=1:8    
a(i)=digest3(i);
XXXX(i)=  cellstr(a(i));
b(i)=hex2dec(XXXX(i));
end

y_03=(bitxor((bitxor((bitxor(b(1),b(8))),b(2))),b(7))+bitxor((bitxor((bitxor(b(3),b(6))),b(4))),b(5)))/256;
% y_03=0.0234;
%y_03=0.2;
for i=1:8    
a(i)=digest4(i);
XXXX(i)=  cellstr(a(i));
b(i)=hex2dec(XXXX(i));
end

y_04=(bitxor((bitxor((bitxor(b(1),b(8))),b(2))),b(7))+bitxor((bitxor((bitxor(b(3),b(6))),b(4))),b(5)))/256;
%y_04=0.0625;
for i=1:8    
a(i)=digest5(i);
XXXX(i)=  cellstr(a(i));
b(i)=hex2dec(XXXX(i));
end

y_05=(bitxor((bitxor((bitxor(b(1),b(8))),b(2))),b(7))+bitxor((bitxor((bitxor(b(3),b(6))),b(4))),b(5)))/256;
%y_05=0.0977;
for i=1:8    
a(i)=digest6(i);
XXXX(i)=  cellstr(a(i));
b(i)=hex2dec(XXXX(i));
end

y_06=(bitxor((bitxor((bitxor(b(1),b(8))),b(2))),b(7))+bitxor((bitxor((bitxor(b(3),b(6))),b(4))),b(5)))/256;
%y_06=0.0195;
for i=1:8    
a(i)=digest7(i);
XXXX(i)=  cellstr(a(i));
b(i)=hex2dec(XXXX(i));
end

y_07=(bitxor((bitxor((bitxor(b(1),b(8))),b(2))),b(7))+bitxor((bitxor((bitxor(b(3),b(6))),b(4))),b(5)))/256;
%y_07=0.1016;
for i=1:8    
a(i)=digest8(i);
XXXX(i)=  cellstr(a(i));
b(i)=hex2dec(XXXX(i));
end

y_08=(bitxor((bitxor((bitxor(b(1),b(8))),b(2))),b(7))+bitxor((bitxor((bitxor(b(3),b(6))),b(4))),b(5)))/256;
%y_08=0.0898;
%% 像素置乱
   [lorenze_x,lorenze_y,lorenze_z]=Lorenz_chaotic_system(y_02,y_03,y_04,m*n);
 lorenze_x=lorenze_x(3002:length(lorenze_x));
 lorenze_y=lorenze_y(3002:length(lorenze_y)); 
 lorenze_z=lorenze_z(3002:length(lorenze_z)); 
  %         X3=pixelEncryption(X1,lorenze_x);
  %         X4=pixelEncryption(X2,lorenze_y);
   
 %   X3=encrypt(X1,key2);
 %   X4=encrypt(X2,key2);

%%  小波变换矩阵生成
ww=DWT(m);
[c,d]=size(ww);

u=3.9999;     %Logistic参数μ，自定为3.99
logid=3;
logix=zeros(1,logid*c+1000);        %预分配内存

logix(1)=y_01;


for i=1:logid*c+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
    logix(i+1)=u*logix(i)*(1-logix(i));
end

logix=logix(1001:logid:length(logix));            %去除前1000点，获得更好的随机性（length数组长度）


[~,Ux]=sort(logix,'descend');      %（对矩阵进行降序排列，Ux为排序后元素在原矩阵的列位置）
%[~,Uy]=sort(logiy,'descend');  

for i=1:c   %行置换
    temp = ww(i,:);
    ww(i,:) = ww(Ux(i),:);
    ww(Ux(i),:) = temp;    
end
% for i=1:d   %列置换
%     temp = ww(:,i);
%     ww(:,i) = ww(:,Uy(i));
%     ww(:,Uy(i)) = temp;    
% end


%  小波变换让图像稀疏化（注意该步骤会耗费时间，但是会增大稀疏度）
X5=ww*sparse(X1);   %   sparse将矩阵x转换为稀疏矩阵，矩阵中去除零元素
X5=full(X5);           %   full 把稀疏矩阵转换成一个全矩阵
X6=ww*sparse(X2);
X6=full(X6);
%X1=mat2gray(X1);
% imshow(uint8(X1));
%   figure(1);
  %subplot(2,2,1); imshow(uint8(X));
  %subplot(2,2,1); imshow(uint8(Y));
  %figure(2);
 % subplot(1,2,1); imshow(uint8(X5));
 % subplot(1,2,2);imshow(uint8(X6));
% spa=0;
% k=0;
% while(1)
% for mm=1:m
%     for nn=1:n
%         if X5(mm,nn)<1-k && X5(mm,nn)>-1+k   
%               spa=spa+1;
%         end  
%     end
% end
%   k=k+0.1;
% if(spa<m*n*0.01)
%     break;
% else
%     spa=0;
% end
% end
% 
% for mm=1:m
%     for nn=1:n
%        if X5(mm,nn)<1-k && X5(mm,nn)>-1+k
%            X5(mm,nn)=0;
%        end
%     end
% end
% 
% 
% spa=0;
% k=0;
% while(1)
% for mm=1:m
%     for nn=1:n
%         if X6(mm,nn)<1-k && X6(mm,nn)>-1+k   
%               spa=spa+1;
%         end  
%     end
% end
%   k=k+0.1;
% if(spa<m*n*0.01)
%     break;
% else
%     spa=0;
% end
% end
% 
% for mm=1:m
%     for nn=1:n
%        if X6(mm,nn)<1-k && X6(mm,nn)>-1+k
%            X6(mm,nn)=0;
%        end
%     end
% end

% for Q=1:((m/2)*n) 
%     if(lorenze_z==0)
%         measure(Q)=lorenze_x(Q);
%     else
%         measure(Q)=lorenze_y(Q);
%     end
% end
% measure=sqrtm(4/m)*measure;
% measure=reshape(measure,(m/2),n);

  Y1=zeros(m/2,n);
  Y2=zeros(m/2,n);
  for i=1:n
     x0=y_03+0.00001*i;y0=y_01; 
    [heon_x,heon_y]=henon(x0,y0,(m/2)*n-1);
     R1=reshape(heon_x,m/2,n);
     R2=reshape(heon_y,m/2,n);
     Y1(:,i)=R1*X5(:,i);
     Y2(:,i)=R2*X6(:,i);
  end

% R=randn(m/2,n);    % randn生成M*a的标准正态分布的随机矩阵 测量矩阵




%Y1=measure*X3;
% Y1=R*X5;
% Y2=R*X6;
Y3=[Y1;Y2];


[max_1,~]=max(Y3);
[max_1,~]=max(max_1);
[min_1,~]=min(Y3);
[min_1,~]=min(min_1);


for i=1:(m)
    for j=1:n
        Y(i,j)=round((Y3(i,j)-min_1)/(max_1-min_1)*255);
     %  W(i,j)= Q(i,j)*(max_1-min_1)/256+min_1;
    end
end

%Y=mod(round(Y*10^4),256);


%  随机矩阵生成
%M=190;
%R=randn(M,n);    % randn生成M*a的标准正态分布的随机矩阵 测量矩阵
%Y=R*X1;
%%  测量
%figure(1);
%imshow(uint8(Y));



%% 补零
%将图像的行列数都补成可以被t整除的数，t为分块的大小。
t=4;
M1=mod((m),t);    %可作为固定密钥，以便解码时可以去除补上的0(mod求余数）
N1=mod(n,t);    %可作为固定密钥，以便解码时可以去除补上的0(mod求余数）
if M1~=0
    Y((m)+1:(m)+t-M1,:)=0; 
end
if N1~=0
    Y(:,n+1:n+t-N1)=0; 
end
[M,N]=size(Y);  %补零后的行数和列数
SUM=M*N;

% logi2=zeros(1,M*N+1000);        %预分配内存
% logi2(1)=y_02;
% for i=1:M*N+999                 %进行SUM+999次循环，共得到SUM+1000点（包括初值）
%     logi2(i+1)=u*logi2(i)*(1-logi2(i));
% end
lorenze_z=lorenze_z(1:SUM); 
p=mod(round(lorenze_z*10^4),256); %（round四舍五入）（mod除后的余数）
bian_R=reshape(p,M,N);  %reshape转成M行N列的随机矩阵R

block=N/t;  %e表示每一行可以分为多少块
r_block=(M/t)*(N/t);

[xi,yi,zi,wi]=siweihundun(y_05,y_06,y_07,y_08,80,45,22,5,21,100,8,r_block);
xi=xi(3002:length(xi));        %去除前3001项，获得更好的随机性（求解陈氏系统的子函数多计算了3000点）
yi=yi(3002:length(yi));
zi=zi(3002:length(zi));
wi=wi(3002:length(wi));

xi=mod(round(xi*10^4),8)+1;  %（round四舍五入）(mod求余数)
yi=mod(round(yi*10^4),8)+1;
zi=mod(round(zi*10^4),4);
wi=mod(round(wi*10^4),8)+1;

for i=1:r_block
    
   Q1=DNA_bian(fenkuai(t,Y,i),xi(i));
   Q2=DNA_bian(fenkuai(t,bian_R,i),yi(i));
   Q3=DNA_yunsuan(Q1,Q2,zi(i));  
   
    xx=floor(i/block)+1; %（floor朝负无穷方向取整）
    yy=mod(i,block);     %（mod取余）
    if  yy==0
        xx=xx-1;
        yy=block;
    end
    Q4((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q3,wi(i));    %将每一块合并成完整的图Q
end
  
 % figure(2);
%  imshow(uint8(Q4));
%% zigzag混沌
  Q5=reshape(Q4,1,SUM);
  zigzagQ=zigzag(M); 
  key=2;    %   randi([1 M*N]);
  Q6=Q5(zigzagQ);
  zigzag1=Q6(1:key-1);
  zigzag2=Q6(key:SUM);
  zigzag3=[zigzag2 zigzag1];
  Q7=reshape(zigzag3,M,N);
  
 
%% 分块混沌
    
%    figure(3);
 %   imshow(uint8(Q7));
%     Q7=uint8(Q7);
%     imwrite(Q7,'../Wavelet_OMP/Q7.bmp','bmp'); 
  % figure;imhist(uint8(Q8));title('原始图片R通道直方图');%imhist获取图像数据直方图
  Q8=pixelEncryption(Q7,lorenze_x);
  
  L = lorenze_y(1:(M/8)*(N/8));
[L_s,H] = sort(L,'descend'); % H can be used as scrambling array.H 可以用作加扰数组。
         Q8_1 =encrypt(Q8,H);
   %     Q8_1=uint8(Q4);
   %      Q8_1=imnoise(Q8_1,'gaussian',0,0.00001);
 toc 
  %figure;imhist(uint8(Q8_1));%imhist获取图像数据直方图
  %figure;imhist(uint8(X1));%imhist获取图像数据直方图
  %figure;imhist(uint8(X2));%imhist获取图像数据直方图
   figure(2);
  imshow(uint8(Q8_1));
  Q8_1=uint8(Q8_1);
 % Q8_1=imnoise(Q8_1,'salt & pepper',0.0000001);
 % imwrite(uint8(Q8_1),'../Wavelet_OMP/DICOM3.bmp','bmp');    %将图片以bmp形式保存    
 
 % P1=imread('Q8_1.bmp');
 % P2=imread('Q8_4.bmp');
 %  nu=NPCRUACIBACI(P2,Q8_1);
 %% 图像解密
%    [lorenze_x,lorenze_y,lorenze_z]=Lorenz_chaotic_system(y_02+10^-15,y_03,y_04,m*n);
%   lorenze_x=lorenze_x(3002:length(lorenze_x));
%   lorenze_y=lorenze_y(3002:length(lorenze_y)); 
%   lorenze_z=lorenze_z(3002:length(lorenze_z)); 
%   zi=mod(round(zi*10^4),4);
 
 Q8_1=double(Q8_1);
 %Q11=Q8_1;
 Q8_2=decrypt(double(Q8_1),H);

 



%[M,N]=size(Q8);
%SUM=M*N;

Q9=pixelUnEncryption(Q8_2,lorenze_x);


Q10=reshape(Q9,1,SUM);
zigzag1=Q10(1:SUM-key+1);
zigzag2=Q10(SUM-key+2:SUM);
zigzag3=[zigzag2 zigzag1];
Q11=izigzagScan(zigzag3,M,N);
  
zi(zi==0)=4;      %加减法互换
zi(zi==1)=0;
zi(zi==4)=1;

for i=1:r_block
    Q12=DNA_bian(fenkuai(t,Q11,i),wi(i));

    Q13=DNA_bian(fenkuai(t,bian_R,i),yi(i));
    
    Q14=DNA_yunsuan(Q12,Q13,zi(i));
   
    
    xx=floor(i/block)+1;    %（floor朝负无穷方向取整）
    yy=mod(i,block);        %（mod求余数）
    if yy==0
        xx=xx-1;
        yy=block;
    end
   
    Q15((xx-1)*t+1:xx*t,(yy-1)*t+1:yy*t)=DNA_jie(Q14,xi(i));
end

M1=0;   %加密时补零的参数，M1=mod(M,t);作为密钥
N1=0;   %加密时补零的参数，N1=mod(N,t);作为密钥
if M1~=0
    Q15=Q15(1:M-t+M1,:,:);
end
if N1~=0
    Q15=Q15(:,1:N-t+N1,:);
end

[m,n]=size(Q15);
for i=1:m
    for j=1:n    
       Q16(i,j)= Q15(i,j)*(max_1-min_1)/256+min_1;
    end
end

Q17=Q16(1:m/2,:);
Q18=Q16(m/2+1:m,:);



Q19=zeros(m,n);  %  恢复矩阵
Q20=zeros(m,n);  %  恢复矩阵
for i=1:n  %  列循环 
     x0=y_03+(0.00001)*i;y0=y_01; 
    [heon_x,heon_y]=henon(x0,y0,(m/2)*n-1);
     R3=reshape(heon_x,m/2,n);
     R4=reshape(heon_y,m/2,n);
     rec1=SL0(R3, Q17(:,i), 0.004);
     rec2=SL0(R4, Q18(:,i), 0.004);
%     rec1=omp(Q17(:,i),R3,m);
%     rec2=omp(Q18(:,i),R4,m);
    Q19(:,i)=rec1;
    Q20(:,i)=rec2;
end





Q21=ww'*sparse(Q19);  %  小波反变换
Q21=full(Q21); 
Q22=ww'*sparse(Q20);  %  小波反变换
Q22=full(Q22); 

%Q23=pixelUnEncryption(Q21,lorenze_x);
%Q24=pixelUnEncryption(Q22,lorenze_y);

% Q21=double(uint8(Q21));
% Q22=double(uint8(Q22));


figure(4);
subplot(2,2,1);imshow(uint8(X1));
subplot(2,2,2);imshow(uint8(X2));
subplot(2,2,3);imshow(uint8(Q21));
subplot(2,2,4);imshow(uint8(Q22));

figure(5);
imshow(uint8(Q21));
figure(6);
imshow(uint8(Q22));

% imwrite(uint8(Q21),'../Wavelet_OMP/解密0.00001DICOM1.bmp','bmp');   
% imwrite(uint8(Q22),'../Wavelet_OMP/解密0.00001DICOM2.bmp','bmp');   
 
 
imshow(uint8(Q22));

MSE1=(sum (sum(abs(uint8(X1)-uint8(Q21)).^2)))/(m*n)
psnr1=10*log10(255*255/(MSE1))   %  PSNR
MSE2=(sum (sum(abs(uint8(X2)-uint8(Q22)).^2)))/(m*n)
psnr2=10*log10(255*255/(MSE2))   %  PSNR

 SSIM1=Evaluation_indicators_SSIM(uint8(X1),uint8(Q21))
 SSIM2=Evaluation_indicators_SSIM(uint8(X2),uint8(Q22))

 
%%
%加密图像
T1_R=imhist(uint8(Q8_1));   %统计加密图像灰度值从0~255的分布情况，存至T1
S1_R=sum(T1_R);     %计算整幅图像R通道的灰度值
xxs1_R=0;           %加密图像相关性
%原始图像1
 T1_G=imhist(uint8(X1));
 S1_G=sum(T1_G);
 xxs1_G=0;
%原始图像2
T1_B=imhist(uint8(X2));
 S1_B=sum(T1_B);
 xxs1_B=0;

for i=1:256
    pp1_R=T1_R(i)/S1_R;   %每个灰度值占比，即每个灰度值的概率
    pp1_G=T1_G(i)/S1_G;
    pp1_B=T1_B(i)/S1_B;
    if pp1_R~=0
        xxs1_R=xxs1_R-pp1_R*log2(pp1_R);
    end
    if pp1_G~=0
        xxs1_G=xxs1_G-pp1_G*log2(pp1_G);
    end
    if pp1_B~=0
        xxs1_B=xxs1_B-pp1_B*log2(pp1_B);
    end
end


%% 原始图像相邻像素相关性分析
%{
先随机在0~M-1行和0~N-1列选中5000个像素点，
计算水平相关性时，选择每个点的相邻的右边的点；
计算垂直相关性时，选择每个点的相邻的下方的点；
计算对角线相关性时，选择每个点的相邻的右下方的点。
%}
NN=5000;    %随机取5000对像素点
x1=ceil(rand(1,NN)*(M-1));      %生成5000个1~M-1的随机整数作为行（ceil朝正无穷方向四舍五入）（rand产生在（0，1）之间均匀发布的随机数）
y1=ceil(rand(1,NN)*(N-1));      %生成5000个1~N-1的随机整数作为列
%预分配内存
XX_R_SP=zeros(1,NN);YY_R_SP=zeros(1,NN);        %水平
XX_G_SP=zeros(1,NN);YY_G_SP=zeros(1,NN);
XX_R_CZ=zeros(1,NN);YY_R_CZ=zeros(1,NN);        %垂直
XX_G_CZ=zeros(1,NN);YY_G_CZ=zeros(1,NN);
XX_R_DJX=zeros(1,NN);YY_R_DJX=zeros(1,NN);      %对角线
XX_G_DJX=zeros(1,NN);YY_G_DJX=zeros(1,NN);

for i=1:NN
    %水平
    XX_R_SP(i)=X1(x1(i),y1(i));
    YY_R_SP(i)=X1(x1(i)+1,y1(i));
    XX_G_SP(i)=X2(x1(i),y1(i));
    YY_G_SP(i)=X2(x1(i)+1,y1(i));
    %垂直
    XX_R_CZ(i)=X1(x1(i),y1(i));
    YY_R_CZ(i)=X1(x1(i),y1(i)+1);
    XX_G_CZ(i)=X2(x1(i),y1(i));
    YY_G_CZ(i)=X2(x1(i),y1(i)+1);
    %对角线
    XX_R_DJX(i)=X1(x1(i),y1(i));
    YY_R_DJX(i)=X1(x1(i)+1,y1(i)+1);
    XX_G_DJX(i)=X2(x1(i),y1(i));
    YY_G_DJX(i)=X2(x1(i)+1,y1(i)+1);
end
%水平
figure;scatter(XX_R_SP,YY_R_SP,10,'filled');xlabel('Pixel gray value on location(x,y)');ylabel('Pixel gray value on location(x+1,y)');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
figure;scatter(XX_G_SP,YY_G_SP,10,'filled');xlabel('Pixel gray value on location(x,y)');ylabel('Pixel gray value on location(x+1,y)');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
%垂直
figure;scatter(XX_R_CZ,YY_R_CZ,10,'filled');xlabel('Pixel gray value on location(x,y)');ylabel('Pixel gray value on location(x,y+1)');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
figure;scatter(XX_G_CZ,YY_G_CZ,10,'filled');xlabel('Pixel gray value on location(x,y)');ylabel('Pixel gray value on location(x,y+1)');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);

%对角线
figure;scatter(XX_R_DJX,YY_R_DJX,10,'filled');xlabel('Pixel gray value on location(x,y)');ylabel('Pixel gray value on location(x+1,y+1)');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
figure;scatter(XX_G_DJX,YY_G_DJX,10,'filled');xlabel('Pixel gray value on location(x,y)');ylabel('Pixel gray value on location(x+1,y+1)');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);

%原始图像1
EX1_R=0;EY1_SP_R=0;DX1_R=0;DY1_SP_R=0;COVXY1_SP_R=0;    %计算水平相关性时需要的变量
EY1_CZ_R=0;DY1_CZ_R=0;COVXY1_CZ_R=0;                %垂直
EY1_DJX_R=0;DY1_DJX_R=0;COVXY1_DJX_R=0;             %对角线
%原始图像2
EX1_G=0;EY1_SP_G=0;DX1_G=0;DY1_SP_G=0;COVXY1_SP_G=0;
EY1_CZ_G=0;DY1_CZ_G=0;COVXY1_CZ_G=0;
EY1_DJX_G=0;DY1_DJX_G=0;COVXY1_DJX_G=0;

for i=1:NN
    %第一个像素点的E，水平、垂直、对角线时计算得出的第一个像素点的E相同，统一用EX1表示
    EX1_R=EX1_R+X1(x1(i),y1(i)); 
    EX1_G=EX1_G+X2(x1(i),y1(i)); 
    %第二个像素点的E，水平、垂直、对角线的E分别对应EY1_SP、EY1_CZ、EY1_DJX
    %R通道
    EY1_SP_R=EY1_SP_R+X1(x1(i),y1(i)+1);
    EY1_CZ_R=EY1_CZ_R+X1(x1(i)+1,y1(i));
    EY1_DJX_R=EY1_DJX_R+X1(x1(i)+1,y1(i)+1);
    %G通道
    EY1_SP_G=EY1_SP_G+X2(x1(i),y1(i)+1);
    EY1_CZ_G=EY1_CZ_G+X2(x1(i)+1,y1(i));
    EY1_DJX_G=EY1_DJX_G+X2(x1(i)+1,y1(i)+1);
end
%统一在循环外除以像素点对数1000，可减少运算次数
% 原始图像1
EX1_R=EX1_R/NN;
EY1_SP_R=EY1_SP_R/NN;
EY1_CZ_R=EY1_CZ_R/NN;
EY1_DJX_R=EY1_DJX_R/NN;
% 原始图像2
EX1_G=EX1_G/NN;
EY1_SP_G=EY1_SP_G/NN;
EY1_CZ_G=EY1_CZ_G/NN;
EY1_DJX_G=EY1_DJX_G/NN;

for i=1:NN
    %第一个像素点的D，水平、垂直、对角线时计算得出第一个像素点的D相同，统一用DX表示
    DX1_R=DX1_R+(X1(x1(i),y1(i))-EX1_R)^2;
    DX1_G=DX1_G+(X2(x1(i),y1(i))-EX1_G)^2;

    %第二个像素点的D，水平、垂直、对角线的D分别对应DY1_SP、DY1_CZ、DY1_DJX
    %原始图像1
    DY1_SP_R=DY1_SP_R+(X1(x1(i),y1(i)+1)-EY1_SP_R)^2;
    DY1_CZ_R=DY1_CZ_R+(X1(x1(i)+1,y1(i))-EY1_CZ_R)^2;
    DY1_DJX_R=DY1_DJX_R+(X1(x1(i)+1,y1(i)+1)-EY1_DJX_R)^2;
    %原始图像2
    DY1_SP_G=DY1_SP_G+(X2(x1(i),y1(i)+1)-EY1_SP_G)^2;
    DY1_CZ_G=DY1_CZ_G+(X2(x1(i)+1,y1(i))-EY1_CZ_G)^2;
    DY1_DJX_G=DY1_DJX_G+(X2(x1(i)+1,y1(i)+1)-EY1_DJX_G)^2;

    %两个相邻像素点相关函数的计算，水平、垂直、对角线
    %原始图像1
    COVXY1_SP_R=COVXY1_SP_R+(X1(x1(i),y1(i))-EX1_R)*(X1(x1(i),y1(i)+1)-EY1_SP_R);
    COVXY1_CZ_R=COVXY1_CZ_R+(X1(x1(i),y1(i))-EX1_R)*(X1(x1(i)+1,y1(i))-EY1_CZ_R);
    COVXY1_DJX_R=COVXY1_DJX_R+(X1(x1(i),y1(i))-EX1_R)*(X1(x1(i)+1,y1(i)+1)-EY1_DJX_R);
    %原始图像2
    COVXY1_SP_G=COVXY1_SP_G+(X2(x1(i),y1(i))-EX1_G)*(X2(x1(i),y1(i)+1)-EY1_SP_G);
    COVXY1_CZ_G=COVXY1_CZ_G+(X2(x1(i),y1(i))-EX1_G)*(X2(x1(i)+1,y1(i))-EY1_CZ_G);
    COVXY1_DJX_G=COVXY1_DJX_G+(X2(x1(i),y1(i))-EX1_G)*(X2(x1(i)+1,y1(i)+1)-EY1_DJX_G);
end
%统一在循环外除以像素点对数1000，可减少运算次数
%原始图像1
DX1_R=DX1_R/NN;
DY1_SP_R=DY1_SP_R/NN;
DY1_CZ_R=DY1_CZ_R/NN;
DY1_DJX_R=DY1_DJX_R/NN;
COVXY1_SP_R=COVXY1_SP_R/NN;
COVXY1_CZ_R=COVXY1_CZ_R/NN;
COVXY1_DJX_R=COVXY1_DJX_R/NN;
%原始图像2
DX1_G=DX1_G/NN;
DY1_SP_G=DY1_SP_G/NN;
DY1_CZ_G=DY1_CZ_G/NN;
DY1_DJX_G=DY1_DJX_G/NN;
COVXY1_SP_G=COVXY1_SP_G/NN;
COVXY1_CZ_G=COVXY1_CZ_G/NN;
COVXY1_DJX_G=COVXY1_DJX_G/NN;
%水平、垂直、对角线的相关性
%原始图像1
RXY1_SP_R=COVXY1_SP_R/sqrt(DX1_R*DY1_SP_R);
RXY1_CZ_R=COVXY1_CZ_R/sqrt(DX1_R*DY1_CZ_R);
RXY1_DJX_R=COVXY1_DJX_R/sqrt(DX1_R*DY1_DJX_R);
%原始图像2
RXY1_SP_G=COVXY1_SP_G/sqrt(DX1_G*DY1_SP_G);
RXY1_CZ_G=COVXY1_CZ_G/sqrt(DX1_G*DY1_CZ_G);
RXY1_DJX_G=COVXY1_DJX_G/sqrt(DX1_G*DY1_DJX_G);


%% 加密图像相邻图像相关性分析
%{
先随机在0~M-1行和0~N-1列选中1000个像素点，
计算水平相关性时，选择每个点的相邻的右边的点；
计算垂直相关性时，选择每个点的相邻的下方的点；
计算对角线相关性时，选择每个点的相邻的右下方的点。
%}
%相关性曲线
%水平
XX_R_SP=zeros(1,NN);YY_R_SP=zeros(1,NN);  %预分配内存（NN随机5000对像素点）

%垂直
XX_R_CZ=zeros(1,NN);YY_R_CZ=zeros(1,NN);  %预分配内存

%对角线
XX_R_DJX=zeros(1,NN);YY_R_DJX=zeros(1,NN);  %预分配内存

for i=1:NN
    %水平
    XX_R_SP(i)=Q8_1(x1(i),y1(i));
    YY_R_SP(i)=Q8_1(x1(i)+1,y1(i));
    %垂直
    XX_R_CZ(i)=Q8_1(x1(i),y1(i));
    YY_R_CZ(i)=Q8_1(x1(i),y1(i)+1);
    %对角线
    XX_R_DJX(i)=Q8_1(x1(i),y1(i));
    YY_R_DJX(i)=Q8_1(x1(i)+1,y1(i)+1);
end
%水平
figure;scatter(XX_R_SP,YY_R_SP,10,'filled');xlabel('Pixel gray value on location(x,y)');ylabel('Pixel gray value on location(x+1,y)');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
%垂直
figure;scatter(XX_R_CZ,YY_R_CZ,10,'filled');xlabel('Pixel gray value on location(x,y)');ylabel('Pixel gray value on location(x,y+1)');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
%对角线
figure;scatter(XX_R_DJX,YY_R_DJX,10,'filled');xlabel('Pixel gray value on location(x,y)');ylabel('Pixel gray value on location(x+1,y+1)');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);
%加密图像
EX2_R=0;EY2_SP_R=0;DX2_R=0;DY2_SP_R=0;COVXY2_SP_R=0;    %水平
EY2_CZ_R=0;DY2_CZ_R=0;COVXY2_CZ_R=0;    %垂直
EY2_DJX_R=0;DY2_DJX_R=0;COVXY2_DJX_R=0;   %对角线

for i=1:NN
    %第一个像素点的E，水平、垂直、对角线时计算得出的第一个像素点的E相同，统一用EX2表示
    EX2_R=EX2_R+Q8_1(x1(i),y1(i));
    %第二个像素点的E，水平、垂直、对角线的E分别对应EY2_SP、EY2_CZ、EY2_DJX
    %R通道
    EY2_SP_R=EY2_SP_R+Q8_1(x1(i),y1(i)+1);
    EY2_CZ_R=EY2_CZ_R+Q8_1(x1(i)+1,y1(i));
    EY2_DJX_R=EY2_DJX_R+Q8_1(x1(i)+1,y1(i)+1);
end
%统一在循环外除以像素点对数1000，可减少运算次数
%加密图像
EX2_R=EX2_R/NN;
EY2_SP_R=EY2_SP_R/NN;
EY2_CZ_R=EY2_CZ_R/NN;
EY2_DJX_R=EY2_DJX_R/NN;

for i=1:NN
    %第一个像素点的D，水平、垂直、对角线时计算得出第一个像素点的D相同，统一用DX2表示
    DX2_R=DX2_R+(Q8_1(x1(i),y1(i))-EX2_R)^2;
    %第二个像素点的D，水平、垂直、对角线的D分别对应DY2_SP、DY2_CZ、DY2_DJX
    %加密图像
    DY2_SP_R=DY2_SP_R+(Q8_1(x1(i),y1(i)+1)-EY2_SP_R)^2;
    DY2_CZ_R=DY2_CZ_R+(Q8_1(x1(i)+1,y1(i))-EY2_CZ_R)^2;
    DY2_DJX_R=DY2_DJX_R+(Q8_1(x1(i)+1,y1(i)+1)-EY2_DJX_R)^2;
    %两个相邻像素点相关函数的计算，水平、垂直、对角线
    %加密图像
    COVXY2_SP_R=COVXY2_SP_R+(Q8_1(x1(i),y1(i))-EX2_R)*(Q8_1(x1(i),y1(i)+1)-EY2_SP_R);
    COVXY2_CZ_R=COVXY2_CZ_R+(Q8_1(x1(i),y1(i))-EX2_R)*(Q8_1(x1(i)+1,y1(i))-EY2_CZ_R);
    COVXY2_DJX_R=COVXY2_DJX_R+(Q8_1(x1(i),y1(i))-EX2_R)*(Q8_1(x1(i)+1,y1(i)+1)-EY2_DJX_R);
end
%统一在循环外除以像素点对数1000，可减少运算次数
%加密图像
DX2_R=DX2_R/NN;
DY2_SP_R=DY2_SP_R/NN;
DY2_CZ_R=DY2_CZ_R/NN;
DY2_DJX_R=DY2_DJX_R/NN;
COVXY2_SP_R=COVXY2_SP_R/NN;
COVXY2_CZ_R=COVXY2_CZ_R/NN;
COVXY2_DJX_R=COVXY2_DJX_R/NN;
%水平、垂直、对角线的相关性
%加密图像
RXY2_SP_R=COVXY2_SP_R/sqrt(DX2_R*DY2_SP_R);
RXY2_CZ_R=COVXY2_CZ_R/sqrt(DX2_R*DY2_CZ_R);
RXY2_DJX_R=COVXY2_DJX_R/sqrt(DX2_R*DY2_DJX_R);

%% NPCR UACI
 nu1=NPCRUACIBACI(uint8(X1),uint8(Q21))
 D=(X1~=Q21);
 nu(1)=sum(sum(D))/(M*N)*100;
 nu(2)=sum(sum(abs(X1-Q21)))/(255*M*N);
 NPCR1=nu1(1);
 UACI1=nu1(2);
 nu2=NPCRUACIBACI(X2,Q22)
 NPCR2=nu2(1);
 UACI2=nu2(2);

%%
disp('信息熵：'); 
disp(['原始图片1信息熵=',num2str(xxs1_G),'原始图片2信息熵=',num2str(xxs1_B)]); 
disp(['加密图片信息熵=',num2str(xxs1_R)]); 
disp(['图片1PSNR=',num2str(psnr1)  '  图片1PSNR=',num2str(psnr2) ]); 
disp(['图片1SSIM=',num2str(SSIM1)  '  图片1SSIM=',num2str(SSIM12) ]); 
disp(['加密图片相关性：','   水平相关性=',num2str(RXY2_SP_R),'    垂直相关性=',num2str(RXY2_CZ_R),'  对角线相关性=',num2str(RXY2_DJX_R)]); 
disp(['原始图片1相关性：','  水平相关性=',num2str(RXY1_SP_R),'    垂直相关性=',num2str(RXY1_CZ_R),'  对角线相关性=',num2str(RXY1_DJX_R)]);
disp(['原始图片2相关性：','  水平相关性=',num2str(RXY1_SP_G),'    垂直相关性=',num2str(RXY1_CZ_G),'  对角线相关性=',num2str(RXY1_DJX_G)]);


Q19=uint8(Q19);
Q20=uint8(Q20);
imwrite(Q19,'../Wavelet_OMP/Q19.bmp','bmp');    %将图片以bmp形式保存  
imwrite(Q20,'../Wavelet_OMP/Q20.bmp','bmp');    %将图片以bmp形式保存  

%%  OMP算法
X2=zeros(a,b);  %  恢复矩阵
for i=1:b  %  列循环       
    rec=omp(Y(:,i),measure,a);
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
X3=ww'*sparse(X2)*ww;  %  小波反变换
X3=full(X3);
imshow(uint8(X3));
title('恢复的图像');

%  误差(PSNR)
errorx=sum(sum(abs(X3-X).^2));        %  MSE误差
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
    
    if (norm(r_n)<9)                              %  残差足够小  2范数
        break;
    end
end
hat_y(pos_array)=aug_y;                           %  重构的向量
