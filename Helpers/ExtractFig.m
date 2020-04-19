fig = openfig('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\manual vs jobs test4.fig');
fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;
x = dataObjs(1).XData;
y = dataObjs(1).YData;
figure()
plot(x,y)
peakvaluest = [700	10.567
720	10.391
740	10.818
760	10.257
780	9.9644
800	9.913
820	10.619
840	10.263
860	10.699
880	10.439
900	10.082
920	9.7973
940	9.4976
960	9.6885
980	9.6965
1000	9.8622
1020	12.109
1040	14.421
1060	13.584
1080	13.492]';
save('D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\WL-Power-Jobs.mat','peakvaluest')