%% Init
close all
clear
clc

WaveLengths = 700:20:1080;
PowerPcts = 0:10:100;
PowerReads = zeros(length(WaveLengths),11);
%%
%close all
clc
[~,data] = xlsread("20200121 sample960nm 512x12 PFSout no interval.csv");
wavelength = strsplit(data{4,:},';');
wavelength = str2num(wavelength{2});
datavalue = string(data(8:end-20));
value = split(datavalue,';');
value_num = double(value(:,3))*10^3;
figure()
hold on
plot(1:length(value_num),value_num)
xlabel("sample#")
ylabel("Power/mW")

[pks,locs] = findpeaks(-value_num,0.5,'MinPeakProminence',20,'MinPeakDistance',12200);
loc2 = floor(locs/2);
 
%plot(loc2,value_num(loc2),"o");
title(strcat("WaveLength = ",num2str(wavelength),"nm"))

peakvalues = zeros(1,length(loc2));
powerTrace = zeros(1,length(value_num));
frame1 = find(PowerSpecDiff>4,1);
for i = 1:length(loc2)
    peakvalues(i) = max(value_num(loc2(i)-5000:loc2(i)));
    if i>1
        powerTrace(loc2(i-1):loc2(i))=peakvalues(i);
    else
        powerTrace(frame1:loc2(i))=peakvalues(i);
    end
end
plot(1:length(value_num),powerTrace);

Waven = (wavelength-680)/20;
PowerReads(Waven,2:end) = peakvalues;

%%
figure(2)
hold on
subplot(5,4,Waven)
hold on
plot(1:length(value_num),value_num,'Color',[0, 0.4470, 0.7410])
xlabel("sample#")
ylabel("Power/mW")
%plot(loc2,value_num(loc2),"ro");
plot(loc2,peakvalues,"-",'Color',[0.8500, 0.3250, 0.0980]);
plot(loc2,peakvalues,"*",'Color',[0.9290, 0.6940, 0.1250]);
title(strcat("WaveLength = ",num2str(wavelength),"nm"))


%%
figure(10)
clf
hold on
plot(0:10:100,PowerReads,"*-")
legend('700nm','720nm','740nm','760nm','780nm','800nm','820nm','840nm','860nm','880nm','900nm','920nm','940nm','960nm','980nm','1000nm','1020nm','1040nm','1060nm','1080nm')

%%
figure(11)
hold on
Pfit = zeros(20,2);
for i = 1:20;
	Pfit(i,:) = polyfit(PowerPcts(2:11),PowerReads(i,2:11),1);
    plot(PowerPcts(1:11),PowerPcts(1:11).*Pfit(i,1)+Pfit(i,2),'o-')
end

%%
figure(13)
clf
hold on
plot(0:10:100,PowerReads,"*-")
legend('700nm','720nm','740nm','760nm','780nm','800nm','820nm','840nm','860nm','880nm','900nm','920nm','940nm','960nm','980nm','1000nm','1020nm','1040nm','1060nm','1080nm')
Pfit_2 = zeros(20,2);
for i = 1:20;
	Pfit_2(i,:) = polyfit(PowerPcts(1:11),PowerReads(i,1:11),1);
    plot(PowerPcts(1:11),PowerPcts(1:11).*Pfit_2(i,1)+Pfit_2(i,2),'o-')
end
%Pct_for_ten = (10*ones(20,1)-Pfit_2(:,2))./Pfit_2(:,1)
Pct_for_ten = 10*ones(20,1)./PowerReads(:,2)*0.1
%% 
figure(14)
clf
hold on
for i = 1:8
plot(700:20:1080,PowerReads(:,i+1)','LineWidth',2)
end
legend('10%','20%','30%','40%','50%','60%','70%','80%')
xlabel("Wavelength/nm")
ylabel("Power at sample plane/mW ")