% close all;clc;clear;
tic
if exist('RawData','var')
elseif exist('Default_0188_Mode3D.mat','file')
   load('Default_0188_Mode3D.mat','RawData','Spectrum','RangeZ','RangeX','NrRawData');
else
    handle  = OCTFileOpen('Default_0188_Mode3D.oct');
    disp( OCTFileGetProperty(handle, 'AcquisitionMode') );
    disp( OCTFileGetProperty(handle, 'RefractiveIndex') );
    disp( OCTFileGetProperty(handle, 'Comment') );
    disp( OCTFileGetProperty(handle, 'Study') );
    disp( OCTFileGetProperty(handle, 'ExperimentNumber') );
    %获取X方向和z方向的数据长度
    thisList = handle.head.DataFiles;
    node = thisList.DataFile{3};
    SizeZ = str2double(node.Attributes.SizeZ);
    SizeX = str2double(node.Attributes.SizeX);
    RangeZ = str2double(node.Attributes.RangeZ);
    RangeX = str2double(node.Attributes.RangeX);
    NrRawData = OCTFileGetNrRawData(handle);
    %reading spectral raw data
    RawData = cell(NrRawData-1,1);
    for i = 1:NrRawData-1
         [Data, Spectrum] = OCTFileGetRawData(handle, i);
         RawData{i} = Data; 
    end
    save('Default_0188_Mode3D.mat','RawData','Spectrum','RangeZ','RangeX','NrRawData','RangeX','RangeZ');
end
toc
%--------------------------------------------------------------------------
load('Wavelength.mat','Wavelength');
%reading spectral raw data
% NrRawData = OCTFileGetNrRawData(handle);
% [RawData, Spectrum] = OCTFileGetRawData(handle, 0);
RawData1 = RawData{1};
[row, column] = size(RawData1);
Nifft = 2^15;
%二维层析图像初始化
gray_value = zeros(Nifft/2,column);
%波数空间采样及线性矫正
[sigma,sigma_inter] = linear_correction(Wavelength);
for i = 1:column
   [z,intensity_ifft] = A_scan(sigma,sigma_inter,Spectrum,RawData1(:,i),Nifft);
   grayscale = inten2gray(intensity_ifft);
   gray_value(:,i) = grayscale';
end
%计算峰值位置和FWHM（采用的最后一组数据）
threshold = max(intensity_ifft)/6;
Measuredpeaks = calculate_FWHM(z, intensity_ifft, 0.00005, threshold, 7,20,3);
figure(1);plot(z,intensity_ifft);
text(Measuredpeaks(:,2),Measuredpeaks(:,3),num2str(Measuredpeaks(:,1)));
%获取X方向和z方向的范围
% thisList = handle.head.DataFiles;
% node = thisList.DataFile{3};
% RangeZ = str2double(node.Attributes.RangeZ);
% RangeX = str2double(node.Attributes.RangeX);
%显示图像
figure(2);imshow(gray_value,[],'XData', [0, RangeX], 'YData', [0, RangeZ],'InitialMagnification','fit');
set(gca, 'YDir', 'reverse');
%将小Tick打开
set(gca,'XMinorTick','on'),xlabel('X(mm)');
set(gca,'YMinorTick','on'),ylabel('Z(mm)');
axis on;
axis image

%--------------------------------------------------------------------------
function [sigma, sigma_inter] = linear_correction(Wavelength)
    %波数空间采样及线性矫正    
    Wavelength = Wavelength';
    wavelength(length(Wavelength):-1:1) = Wavelength(1:length(Wavelength));
    %波长映射到波数
    sigma = 1./wavelength;
    %波数均匀采样
    sigma_inter = linspace(sigma(1),sigma(end),length(sigma));
end

function [z,intensity_ifft] = A_scan(sigma,sigma_inter,Spectrum,RawData,Nifft)
    spectrum = RawData(:,1);
    %消除直流项对ifft的影响
    spectrum = spectrum - Spectrum;
    %波数空间采样及线性矫正
    %将光谱强度逆序排列
    intensity(length(spectrum):-1:1) = spectrum(1:length(spectrum));
    intensity_inter = spline(sigma,intensity,sigma_inter);
    %干涉光谱整形（采用高斯窗）
    intensity_shaping = intensity_inter.* (gausswin(length(intensity_inter)))';
    %对光谱进行傅里叶逆变换
    interference_signal = ifft(intensity_shaping,Nifft);
%     interference_signal = ifftshift(interference_signal);
    deltaSigma = (sigma_inter(end)-sigma_inter(1))/(length(sigma_inter)-1);
    deltaZ = 1/(2*Nifft*deltaSigma);
    z = (1:Nifft)*deltaZ;
    intensity_ifft = abs(interference_signal(1:Nifft/2));
    %坐标变换
    z = z-z(end)/2;z = z(Nifft/2+1:end);
%     plot(z,intensity_ifft);
end

function grayscale = inten2gray(intensity_ifft)
    %幅度值转化为灰度值
    intensitymax_ifft = max(intensity_ifft);
    intensitymin_ifft = min(intensity_ifft);
    delta_grayscale = (intensitymax_ifft-intensitymin_ifft)/256;
    grayscale = round(intensity_ifft/delta_grayscale);
end

function Measuredpeaks = calculate_FWHM(x, y, SlopeThreshold, AmpThreshold, SmoothWidth,FitWidth,smoothtype)
%计算Full width at half maximum(FWHM)
figure(1);plot(x,y)  % Graph the signal in red
title('Detected peaks are numbered. Peak table is printed in Command Window')

% Initial values of variable parameters

% SlopeThreshold=.0001;
% AmpThreshold=1;
% SmoothWidth=20;
% FitWidth=18;
% Label the x-axis with the parameter values
xlabel(['SlopeThresh. = ' num2str(SlopeThreshold) '    AmpThresh. = ' num2str(AmpThreshold) '    SmoothWidth = ' num2str(SmoothWidth) '    FitWidth = ' num2str(FitWidth) ])

% Find the peaks
% tic;
Measuredpeaks=findpeaksG(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,smoothtype);
% ElapsedTime=toc;
% PeaksPerSecond=length(Measuredpeaks)/ElapsedTime;

% Display results
disp('---------------------------------------------------------')
disp(['SlopeThreshold = ' num2str(SlopeThreshold) ] )
disp(['AmpThreshold = ' num2str(AmpThreshold) ] )
disp(['SmoothWidth = ' num2str(SmoothWidth) ] )
disp(['FitWidth = ' num2str(FitWidth) ] )
% disp(['Speed = ' num2str(round(PeaksPerSecond)) ' Peaks Per Second' ] )
disp('         Peak #     Position      Height      Width       Area')
disp(Measuredpeaks)  % Display table of peaks
figure(1);text(Measuredpeaks(:,2),Measuredpeaks(:,3),num2str(Measuredpeaks(:,1)))  % Number the peaks found on the graph
end
