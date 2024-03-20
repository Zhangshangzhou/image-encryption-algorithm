function y=zigzag(size)
% 返回size*size矩阵的zigzag范围顺序
% 如2*2，则返回1、3、2、4，注意matlab的顺序和是以列为主
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
