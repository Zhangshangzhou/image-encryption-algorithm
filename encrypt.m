function [ encrypted_symbolized ] = encrypt( original_dct,H )
[M,N] = size(original_dct);
%%Encryption

BLOCK = mat2cell(original_dct,ones(M/8,1)*8,ones(N/8,1)*8);   %one产生N*N的全1矩阵 mat2cell把矩阵变成元胞数组
BLOCKsize = size(BLOCK);
SCAMBLED = cell(BLOCKsize);

for i = 1:length(H)
    SCAMBLED(i) = BLOCK(H(i)); % scambling blocks
end

encrypted = cell2mat(SCAMBLED);
 encrypted_symbolized =encrypted;

% encrypted_symbolized = S.*encrypted;    

end

