function imgOutput=pixelEncryption(imgInput,base)
[rows,clos]=size(imgInput);
sequence=base;
imgOutput=imgInput;
%���ܣ���ÿ�����صĶ����ƽ���
for i=1:clos    %��
    seqClos=sequence((i-1)*8+1:i*8,1);
    [~,sortOrder]=sort(seqClos);
    [~,sortOrder]=sort(sortOrder);        %����
    for j=1:rows    %��
        pixelValue=imgInput(j,i);
        pixelBin=dec2bin(pixelValue,8); %8λ�Ķ�����
        tempBin=pixelBin;
        for k=1:8
            tempBin(1,k)=pixelBin(1,sortOrder(k));
        end
        imgOutput(j,i)=bin2dec(tempBin);
    end
end
