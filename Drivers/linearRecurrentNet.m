%% Problem 1: Deterministic Recurrent Network
close all
options = odeset('RelTol',1e-5,'NonNegative',[1,2,3]);
a = 3
tLen = 20; % simulation length
[t,y]=ode45(@(t,r) lrn(t,r,a),[0 tLen],[1; 1]);
figure(1) %ode method
%subplot(2,2,1)
plot(t,y(:,1),'-',t,y(:,2),'-')
xlabel('t')
ylabel('r(t)')
legend('r_{1}(t)','r_{2}(t)')
title('Deterministic linear recurrent network')
%title(strcat('a = ',num2str(a)))

figure(2)
%subplot(2,2,1)
plot(y(:,1),y(:,2))
xlabel('r_1')
ylabel('r_2')
title('Deterministic linear recurrent network')
%title(strcat('a = ',num2str(a)))

figure(3)
hold on
r1 = exp(-t).*(cos(a*t)-sin(a*t));
r2 = exp(-t).*(cos(a*t)+sin(a*t));
plot(t,r1)
plot(t,r2)
xlabel('t')
ylabel('r(t)')
legend('r_{1}(t)','r_{2}(t)')
title('Deterministic linear recurrent network')

%% Stochastic recurrent Network
%close all
sig = 1;
tau = 0.01;  
tRange = 0:tau:tLen;
W = [0, -a; a, 0];
A = W - eye(2);
K = 1000;
% Initiate the neurons
r = ones(K,2,numel(tRange));
for i = 1:K
    for tt = 2:length(tRange)
        r(i,:,tt)= A*reshape(r(i,:,tt-1),2,1)*tau+reshape(r(i,:,tt-1),2,1)+ normrnd(0,sig,[2,1]);
    end  
end
mu = reshape(mean(r,1),2,numel(tRange));
sigma = cell(1,numel(tRange));
for i = 1:numel(tRange)
    sigmaT = zeros(2,2);
    for j = 1:K
        sigmaT = sigmaT +(r(j,:,i)' - mu(:,i))*(r(j,:,i)' - mu(:,i))';
    end
    sigma{i} = sigmaT/K;
end

figure(4)
subplot(1,3,1)
hold on
plot(tRange,mu(1,:),'-',tRange,mu(2,:),'-')
xlabel('t')
ylabel('r')
legend('r1','r2')
title({'r vs t';strcat('sigma = ',num2str(sig),' tau = ',num2str(tau),' a = ',num2str(a))})

figure(4)
subplot(1,3,2)
hold on
plot(mu(1,:),mu(2,:),'lineWidth',1.5)
xlabel('r1')
ylabel('r2')
title({'r_{1} vs r_{2}';strcat('sigma = ',num2str(sig),' tau = ',num2str(tau),' a = ',num2str(a))})

for i = 1:numel(tRange)
    plotCovarianceEllipse(mu(:,i),sigma{i})
end
%axis equal

figure(4)
subplot(1,3,3)
lbd = zeros(2,numel(tRange));
for i = 1:numel(tRange)
    D = eig(sigma{i});
    lbd(1,i) = max(D);
    lbd(2,i) = min(D);
end
hold on
plot(tRange,lbd(1,:))
plot(tRange,lbd(2,:))
xlabel('t')
ylabel('eigenvalues of \Sigma_{t}')
legend('major axis length','minor axis length')
title({'Axis length of ellipse';strcat('sigma = ',num2str(sig),' tau = ',num2str(tau),' a = ',num2str(a))})

%plot(tRange,lbd(1,:).*lbd(2,:))
    
% question 1 ode
function LV = lrn(t,r,a)
    LV = [-r(1)-a*r(2); a*r(1)-r(2)];
end

% question 2 plot Ellipse
function plotCovarianceEllipse(mu,Sigma) 
    [V,D]=eig(Sigma); % eigenvects and vals 
    theta_grid = linspace(0,2*pi ,40); % 40 points around ellipse 
    a=sqrt(D(2,2)); % major axis length 
    b=sqrt(D(1,1)); % minor axis length 
    ellipse_x = a*cos( theta_grid ); % axis -aligned ellipse 
    ellipse_y = b*sin( theta_grid ); 
    phi = atan2(V(2,2),V(2,1)); % angle of big eigenvector 
    R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ]; % rotation matrix 
    ellipse = [ellipse_x;ellipse_y]' * R; % rotated ellipse
    plot(ellipse(:,1) + mu(1), ellipse(:,2) + mu(2), '-'); % draw the ellipse
end
