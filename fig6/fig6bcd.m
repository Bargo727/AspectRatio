clear;

global lambda kappa eps N gamma D a td b n

lambda = 1;
kappa  = 0.1;
eps    = 0.1;
gamma  = .1;
N      = 30;
D      = 1;
a      = 10000;
b      = 10000;
td     = 55;
n      = 50;

x_init  = N/3;
h1_init = 0;
h2_init = 0;
z_init  = N/4;

T = 2000;

[t,y] = ode23s(@EffectiveOscillation,[0,T],[x_init;h1_init;h2_init;z_init]);

% z = eps + (1-eps)*(N - (y(:,1)))./(1 + N - y(:,1));
% l = eps + (1-eps)*(y(:,1) )./(1 + (y(:,1)));
% 
% zz = 1./(1 + N - y(:,1));
% ll = 1./(1 + y(:,1));
% 
% figure(4)
% plot(t,zz,'r',t,ll,'b','LineWidth',4)
% set(gca,'fontsize',20)
% xlabel('time')
% ylabel('prot')
% title('h to v')
% h = legend('prot1','prot2');
% set(h,'box','off')
% 
% 
% figure(3)
% plot(t,z,'r',t,l,'b','LineWidth',4)
% set(gca,'fontsize',20)
% xlabel('time')
% ylabel('prot')
% title('v to h')
% h = legend('prot1','prot2');
% set(h,'box','off')

figure(3)
plot(t,y(:,1),'r','LineWidth',4)
set(gca,'fontsize',20)
xlabel('time')
ylabel('position')
h = legend('x(t)');
set(h,'box','off')

%axis([0 T 0 N])

% figure(2)
% plot(t,y(:,2),'r',t,y(:,3),'b','LineWidth',4)
% set(gca,'fontsize',20)
% xlabel('time')
% ylabel('fraction of horizontal cells')
% h = legend('strain 1','strain 2');
% set(h,'box','off')