%% Problem 1: Three-neuron winner-take-all
close all
clear
options = odeset('RelTol',1e-5,'NonNegative',[1,2,3]);
tau = 10;
tLen = 1000; % simulation length
[t,y]=ode45(@(t,E) wta(t,E,tau),[0 tLen],[10; 20; 20],options);
figure(1) %ode method
%subplot(2,2,1)
plot(t,y(:,1),'-',t,y(:,2),'-',t,y(:,3),'-')
xlabel('t')
ylabel('E(t)')
legend('E_{1}(t)','E_{2}(t)','E_{3}(t)')
%title('Deterministic linear recurrent network')
%title(strcat('a = ',num2str(a)))

%%
clear
cmap = colormap(winter);
K1 = 80;
K2 = 80;
K3 = 80;
tau = 10;
[r1,r2,r3] = meshgrid(0:10:100);
u = 1/tau*(-r1+ASfcn(K1-5*r2-5*r3));
v = 1/tau*(-r2+ASfcn(K2-5*r1-5*r3));
w = 1/tau*(-r3+ASfcn(K3-5*r1-5*r2));

figure()
quiver3(r1,r2,r3,u,v,w,'MarkerSize',10);
axis equal
axis([0 100 0 100 0 100])
hold on
sx1 = 25; sy1 = 10; sz1 = 10;
sx2 = 30; sy2 = 60; sz2 = 60;
sx3 = 60; sy3 = 10; sz3 = 100;
hlines1 = streamline(stream3(r1,r2,r3,u,v,w,sx1,sy1,sz1));
set(hlines1,'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
hlines2 = streamline(stream3(r1,r2,r3,u,v,w,sx2,sy2,sz2));
set(hlines2,'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250]);
hlines3 = streamline(stream3(r1,r2,r3,u,v,w,sx3,sy3,sz3));
set(hlines3,'LineWidth',2,'Color',[0.4660, 0.6740, 0.1880]);
FixedPoints = [80,0,0; 0,80,0; 0,0,80; 6.91, 6.91,6.91];
scatter3([sx1,sx2,sx3],[sy1,sy2,sy3],[sz1,sz2,sz3],80,'+r','LineWidth',1.25);
scatter3(FixedPoints(:,1),FixedPoints(:,2),FixedPoints(:,3),'r','LineWidth',1.25);
xlabel('E1')
ylabel('E2')
zlabel('E3')
legend('Vector Field','Trajectory_1','Trajectory_2','Trajectory_3','Initial Conditions','End Points')
set(gca,'linewidth', 2,'fontsize',15) % Sets the width of the axis lines, font size, font

figure()

hold on
%quiver3(r1,r2,r3,u,v,w,3,'MarkerSize',10);
% when E3 = 0
r1 = 0:1:100;
r2 = ASfcn(K1-5*r1);
r3 = zeros(size(r1));
plot3(r1,r2,r3,'r','LineWidth',2)
plot3(r2,r1,r3,'b','LineWidth',2)
plot3(r2,r3,r1,'g','LineWidth',2)
plot3(r1,r3,r2,'r','LineWidth',2)
plot3(r3,r1,r2,'b','LineWidth',2)
plot3(r3,r2,r1,'g','LineWidth',2)
[r1,r2,r3] = meshgrid(0:10:100);

[sx, sy, sz] = meshgrid(0:20:100, 0:20:100, 0:20:100);
hlines = streamline(stream3(r1,r2,r3,u,v,w,sx,sy,sz));
set(hlines,'LineWidth',1,'Color',[0.3010, 0.7450, 0.9330]);

X = [0;100;100];
Y = [0;100;100];
Z = [0;100;0];
fill3(X,Y,Z,'g','FaceAlpha',0.2)
X = [0;100;100];
Y = [0;100;0];
Z = [0;100;100];
fill3(X,Y,Z,'r','FaceAlpha',0.2)
X = [0;100;0];
Y = [0;100;100];
Z = [0;100;100];
fill3(X,Y,Z,'b','FaceAlpha',0.2)
legend('E_{1} nullclines','E_{2} nullclines','E_{3} nullclines')

axis equal
axis([0 100 0 100 0 100])
xlabel('E1')
ylabel('E2')
zlabel('E3')
set(gca,'linewidth', 2,'fontsize',15) % Sets the width of the axis lines, font size, font

%%
    syms E1 E2 E3
    tau = 10;
    Efix1 = 6.91;
    Efix2 = 6.91;
    Efix3 = 6.91;
    K1 = 80;
    K2 = 80;
    K3 = 80;
    defE1 = 1/tau*(-E1+(100*(K1-5*E2-5*E3)^2)/(40^2+(K1-5*E2-5*E3)^2));
    defE2 = 1/tau*(-E2+(100*(K2-5*E1-5*E3)^2)/(40^2+(K2-5*E1-5*E3)^2));
    defE3 = 1/tau*(-E3+(100*(K3-5*E1-5*E2)^2)/(40^2+(K3-5*E1-5*E2)^2));
    J = jacobian([defE1, defE2, defE3],[E1,E2,E3]);
    E1 = Efix1;
    E2 = Efix2;
    E3 = Efix3;
    evalJ = [subs(J(1,1)),subs(J(1,2)),subs(J(1,3));subs(J(2,1)),subs(J(2,2)),subs(J(2,3));subs(J(3,1)),subs(J(3,2)),subs(J(3,3))];
    V = eig(evalJ);
    double(evalJ)
    double(V)
M = 0.1*[-1 -2 -2; 0 -1 0; 0 0 -1]
eig(M)

function Resp = wta(t,E,tau)
    K1 = 80;
    K2 = 80;
    K3 = 80;
    Resp = 1/tau*[-E(1)+Sfcn(K1-5*E(2)-5*E(3)); -E(2)+Sfcn(K2-5*E(1)-5*E(3)); -E(3)+Sfcn(K3-5*E(1)-5*E(2))];
end

function gSfcn = Sfcn(x)
    if x < 0
        gSfcn = 0;
    else
        gSfcn = 100*x.^2./(40^2 + x.^2);
    end
end

function gASfcn = ASfcn(x)
    gASfcn = zeros(size(x));
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            for k = 1:size(x,3)
                if x(i,j,k) < 0;
                    gASfcn(i,j,k) = 0;
                else
                    gASfcn(i,j,k) = 100*x(i,j,k)^2/(40^2 +x(i,j,k)^2);
                end
            end
        end
    end
end

