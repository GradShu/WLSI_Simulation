clear;close all;
handle = OCTFileOpen('D:\Desktop\zms3\jiaodai\Default_0094_Mode2D.oct');%调用数据包
 %nrRawData = OCTFileGetNrRawData( handle );
[RawData, Spectrum] = OCTFileGetRawData(handle, 0);%第一列
% SData(1:2048,1:512)=RawData;
% SData(1:2048,1:512,1:512)=0;
% for i=0:nrRawData-1
% [RawData, Spectrum] = OCTFileGetRawData(handle, i);% 0~511 B-scan
% SData(:,:,i+1)=RawData;
% % figure(4);clf;
% % plot(Spectrum);
% % figure(5);clf;
% % imagesc(RawData);
% % figure(6);clf;
% % plot(RawData(:,512));
% end
save('D:\Desktop\zms3\jiaodai\Surfacejiaodai_9.mat','RawData', '-v7.3') ;