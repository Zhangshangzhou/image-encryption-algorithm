function y=zigzag(size)
% ����size*size�����zigzag��Χ˳��
% ��2*2���򷵻�1��3��2��4��ע��matlab��˳���������Ϊ��
seq = 1:1:size*size;
seq = reshape(seq,size,size);
y = zeros(1,size*size);
k = 1;
for i = 2:2*size
    for j = max(i-size,1):min(size,i-1)
        if mod(i,2)==1
            y(k) = seq(j,i-j);
        else
            y(k) = seq(i-j,j);
        end
        k = k+1;
    end
end
