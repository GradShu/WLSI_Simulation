clear;close all;
handle = OCTFileOpen('D:\Desktop\zms1\714\Default_0026_Mode3D.oct');%�������ݰ�
nrRawData = OCTFileGetNrRawData( handle );
% [RawData, Spectrum] = OCTFileGetRawData(handle, 0);%��һ��
%SData(1:2048,1:512)=RawData;
SData(1:2048,1:1000,1:1000)=0;
for i=1:nrRawData-1
[RawData, Spectrum] = OCTFileGetRawData(handle, i);%����һ�ж�ȡ�������Ǵ��
SData(:,:,i)=RawData;
% figure(4);clf;
% plot(Spectrum);
%  figure(5);clf;
%  imagesc(RawData);
% figure(6);clf;
% plot(RawData(:,512));
end
 save('D:\Desktop\zms1\Surfacemian1000_03.mat','SData', '-v7.3') ;
 