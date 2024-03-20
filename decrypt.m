function [ dec_after_SCAMBLED ] = decrypt( dec_dct,H )
%%decryption
[M,N] = size(dec_dct);


%dec_dct = round(dec_dct./S);



dec_BLOCK = mat2cell(dec_dct,ones(M/8,1)*8,ones(N/8,1)*8);
BLOCKsize = size(dec_BLOCK);
dec_SCAMBLED = cell(BLOCKsize);
%%reverse scambling∑¥œÚº”»≈
for i = 1:length(H)
    dec_SCAMBLED(H(i)) = dec_BLOCK(i); % resume blocksª÷∏¥øÈ
end
dec_after_SCAMBLED = cell2mat(dec_SCAMBLED);


end

