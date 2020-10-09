close all
clear
f1 = openfig('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200918 Power ramping during warm up\Power ramping during start up');
f2 = openfig('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200928 Power adjust 960\Power ramping\Power ramping during start up');
h1 = findobj(f1, '-class', 'matlab.graphics.chart.primitive.Line');
h2 = findobj(f2, '-class', 'matlab.graphics.chart.primitive.Line');
tr1(1,:) = [minutes(h1.XData)- minutes(h1.XData(1))]';
tr2(1,:) = [minutes(h2.XData)- minutes(h2.XData(1))]';

tr1(2,:) = h1.YData./max(h1.YData);
tr2(2,:) = h2.YData./max(h2.YData);
figure()
hold on
plot(tr1(1,:), tr1(2,:))
plot(tr2(1,:), tr2(2,:))
xlabel("Minutes")
ylabel("Power relative")
legend('Trial 1','Trail 2')
