%  ������ʵ��ͼ��LENA��ѹ������
%  �������ߣ�ɳ������۴�ѧ�������ӹ���ѧϵ��wsha@eee.hku.hk
%  �㷨��������ƥ�䷨���ο����� Joel A. Tropp and Anna C. Gilbert 
%  Signal Recovery From Random Measurements Via Orthogonal Matching
%  Pursuit��IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 53, NO. 12,
%  DECEMBER 2007.
%  �ó���û�о����κ��Ż�

function Wavelet_OMP

clc;clear

%  ���ļ�
X=imread('Lena512.bmp');
% X=X(:,:,1);        %R
%imwrite(X,'../Wavelet_OMP/walkbridge.bmp','bmp');    %��ͼƬ��bmp��ʽ����  

% I=ind2gray(X,map);
% figure(2);
% imshow(uint8(X));

%X=imread('4.2.06.tiff');
%X=imread('Rice.tif');
%X= rgb2gray(X);



X=double(X);
[a,b]=size(X);
%X=ones(a,b);

%  С���任��������
ww=DWT(a);
% wavename='db2';
% [cA2,cB2,cD2,cE2]=dwt2(X,wavename);


%  С���任��ͼ��ϡ�軯��ע��ò����ķ�ʱ�䣬���ǻ�����ϡ��ȣ�
X1=ww*sparse(X);   %   sparse������xת��Ϊϡ����󣬾�����ȥ����Ԫ��
X1=full(X1);           %   full ��ϡ�����ת����һ��ȫ����

%  �����������
%  M=256;
%  R=randn(M,a);    % randn����M*a�ı�׼��̬�ֲ���������� ��������
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

 
% u=3.9999;     %Logistic�����̣��Զ�Ϊ3.99
% logid=3;
% c=a;
% logix=zeros(1,logid*c+1000);        %Ԥ�����ڴ�
% 
% logix(1)=0.234;
% for i=1:logid*c+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
%     logix(i+1)=u*logix(i)*(1-logix(i));
% end
% logix=logix(1001:logid:length(logix));            %ȥ��ǰ1000�㣬��ø��õ�����ԣ�length���鳤�ȣ�
% 
% logiy(1)=0.345;
% for i=1:logid*c+999                 %����SUM+999��ѭ�������õ�SUM+1000�㣨������ֵ��
%     logiy(i+1)=u*logiy(i)*(1-logiy(i));
% end
% logiy=logiy(1001:logid:length(logiy));            %ȥ��ǰ1000�㣬��ø��õ�����ԣ�length���鳤�ȣ�
% 
% [~,Ux]=sort(logix,'descend');      %���Ծ�����н������У�UxΪ�����Ԫ����ԭ�������λ�ã�
% [~,Uy]=sort(logiy,'descend');  
% 
% for i=1:a   %���û�
%     temp = R(i,:);
%     R(i,:) = R(Ux(i),:);
%     R(Ux(i),:) = temp;    
% end
% for i=1:b   %���û�
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
% %  ����
% Y=R*X1;

%  OMP�㷨
X2=zeros(a,b);  %  �ָ�����
for i=1:b  %  ��ѭ��   
    x0=0.2+0.0001*i;y0=0.2; 
    [heon_x,heon_y]=henon(x0,y0,(512/2)*512-1);
    R=0.08838834756*reshape(heon_y,a/2,b);
   rec=SL0(R,Y(:,i), 0.004);
    % rec=omp(Y(:,i),R,a);
    X2(:,i)=rec;
end


%  ԭʼͼ��
figure(1);
imshow(uint8(X));
title('ԭʼͼ��');

%  �任ͼ��
figure(2);
imshow(uint8(X1));
title('С���任���ͼ��');

%  ѹ�����лָ���ͼ��
figure(3);
X3=ww'*sparse(X2);  %  С�����任
X3=full(X3);
subplot(2,2,2);imshow(uint8(X3));title('�ָ���ͼ��');

subplot(2,2,1);imshow(uint8(X));title('ԭʼ��ͼ��');

%  ���(PSNR)
errorx=sum (sum(abs(uint8(X3)-uint8(X)).^2))        %  MSE���
psnr=10*log10(255*255/(errorx/a/b))   %  PSNR
 

%  OMP�ĺ���
%  s-������T-�۲����N-������С
function hat_y=omp(s,T,N)

Size=size(T);                                     %  �۲�����С
M=Size(1);                                        %  ����
hat_y=zeros(1,N);                                 %  ���ع�������(�任��)����                     
Aug_t=[];                                         %  ��������(��ʼֵΪ�վ���)
r_n=s;                                            %  �в�ֵ

for times=1:M/4;                                  %  ��������(ϡ����ǲ�����1/4)
    for col=1:N;                                  %  �ָ����������������
        product(col)=abs(T(:,col)'*r_n);          %  �ָ�������������Ͳв��ͶӰϵ��(�ڻ�ֵ) 
    end
    [val,pos]=max(product);                       %  ���ͶӰϵ����Ӧ��λ��
    Aug_t=[Aug_t,T(:,pos)];                       %  ��������
    T(:,pos)=zeros(M,1);                          %  ѡ�е������㣨ʵ����Ӧ��ȥ����Ϊ�˼��Ұ������㣩
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  ��С����,ʹ�в���С
    r_n=s-Aug_t*aug_y;                            %  �в�
    pos_array(times)=pos;                         %  ��¼���ͶӰϵ����λ��
    
    if (norm(r_n)<2)                              %  �в��㹻С  2����
        break;
    end
end
hat_y(pos_array)=aug_y;                           %  �ع�������
