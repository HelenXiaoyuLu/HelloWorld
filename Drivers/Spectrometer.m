clear
close all
D = readcell('D:\SpectraMeter\Spectrum_manual','FileType','spreadsheet');
Mat = cell2mat(D(2:end,:));
Head = D(1,3:end);
waveLen = zeros(1, length(Head)/2);
for i = [1:2:length(Head)-1]
    waveLen(ceil(i/2))= Head{i};
end

Matbgcrt = zeros(size(Mat,1),size(Mat,2)-2);
figure()
hold on
waveLenMeas = zeros(1,length(waveLen));
for i = 4:2:size(Mat,2)
    Matbgcrt(:,(i-2)/2) = Mat(:,i) - Mat(:,2);
    [~,idx] = max(Matbgcrt(:,(i-2)/2));
    waveLenMeas((i-2)/2) = Mat(idx,1);
    plot(Mat(:,1),Matbgcrt(:,(i-2)/2))
end
figure()
hold on
plot(waveLen,waveLen,'-')
scatter(waveLen,waveLenMeas,100,'+')
legend('Expected','SpectroMeter')
xlabel('wavelength')
ylabel('wavelength')

%%
clear
close all
D = readcell('D:\SpectraMeter\Spectrum_manual','FileType','spreadsheet','Sheet','JOBs');
Mat = cell2mat(D(2:end,:));
Head = D(1,3:end);
waveLen = zeros(1, length(Head)/2);
for i = [1:2:length(Head)-1]
    waveLen(ceil(i/2))= Head{i};
end

Matbgcrt = zeros(size(Mat,1),size(Mat,2)-2);
figure()
hold on
waveLenMeas = zeros(1,length(waveLen));
for i = 4:2:size(Mat,2)
    Matbgcrt(:,(i-2)/2) = Mat(:,i) - Mat(:,2);
    [~,idx] = max(Matbgcrt(:,(i-2)/2));
    waveLenMeas((i-2)/2) = Mat(idx,1);
    plot(Mat(:,1),Matbgcrt(:,(i-2)/2))
end
figure()
hold on
plot(waveLen,waveLen,'-')
scatter(waveLen,waveLenMeas,100,'+')
legend('Expected','SpectroMeter')
xlabel('wavelength')
ylabel('wavelength')
