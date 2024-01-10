close all;clc;clear;
handle  = OCTFileOpen('Default_0186_Mode2D.oct');
disp( OCTFileGetProperty(handle, 'AcquisitionMode') );
disp( OCTFileGetProperty(handle, 'RefractiveIndex') );
disp( OCTFileGetProperty(handle, 'Comment') );
disp( OCTFileGetProperty(handle, 'Study') );
disp( OCTFileGetProperty(handle, 'ExperimentNumber') );
%-----------------------------------------------------------
%reading spectral raw data
NrRawData = OCTFileGetNrRawData(handle);
[RawData, Spectrum] = OCTFileGetRawData(handle, 0);
% plot(1:2048,RawData(:,1));

load('Wavelength.mat','Wavelength');
[row, column] = size(RawData);
Nifft = 2^11;
gray_value = zeros(Nifft/2,column);
for i = 1:column
   [z,intensity_ifft] = A_scan(Wavelength,Spectrum,RawData(:,i),Nifft);
   grayscale = inten2gray(intensity_ifft);
   gray_value(:,i) = grayscale';
end
% image(gray_value);
% I = mat2gray(gray_value);
imshow(gray_value,[]);
imwrite(uint8(gray_value), 'matrix.tif' );

%二维层析图像绘制

function [z,intensity_ifft] = A_scan(Wavelength,Spectrum,RawData,Nifft)
    Wavelength = Wavelength';
%     plot(1:length(Wavelength),Wavelength);
%     plot(Wavelength,Spectrum);
%     plot(Wavelength,RawData(:,1));
    spectrum = RawData(:,1);
    %消除直流项对ifft的影响
    spectrum = spectrum - Spectrum;
    %波数空间采样及线性矫正
    %将波长和光谱强度逆序排列
    intensity = zeros(1,length(spectrum));
    wavelength = zeros(1,length(Wavelength));
    for i = 1:length(Wavelength)
        wavelength(length(wavelength)+1-i) = Wavelength(i);
        intensity(length(intensity)+1-i) = spectrum(i);
    end
    %波长映射到波数
    sigma = 1./wavelength;
    %波数均匀采样
    sigma_inter = linspace(sigma(1),sigma(end),length(sigma));
    %光谱强度三次样条插值
    intensity_inter = spline(sigma,intensity,sigma_inter);
%     plot(sigma_inter,intensity_inter);
    %干涉光谱整形（采用汉宁窗）
    intensity_shaping = intensity_inter.* (hann(length(intensity_inter)))';
%     plot(sigma_inter,intensity_shaping);
    %对光谱进行傅里叶逆变换
%     Nifft = 2^15;
    interference_signal = ifft(intensity_shaping,Nifft);
%     interference_signal = ifftshift(interference_signal);
    deltaSigma = (sigma_inter(end)-sigma_inter(1))/(length(sigma_inter)-1);
    deltaZ = 1/(2*Nifft*deltaSigma);
    z = (1:Nifft)*deltaZ;
    intensity_ifft = abs(interference_signal(1:Nifft/2));
    %坐标变换
    z = z-z(end)/2;
%     plot(z,intensity_ifft);
end

function grayscale = inten2gray(intensity_ifft)
    %幅度值转化为灰度值
    intensitymax_ifft = max(intensity_ifft);
    intensitymin_ifft = min(intensity_ifft);
    delta_grayscale = (intensitymax_ifft-intensitymin_ifft)/256;
    grayscale = round(intensity_ifft/delta_grayscale);
end
