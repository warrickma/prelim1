x1=Balances(A,S,1,660); % steady state for first 600 and pahse one for 600-660
t1=0:660;
figure (1);
plot(t1,x1(4,:),'k','LineWidth',2,'DisplayName','P1');
hold on;
plot(t1,x1(5,:),'b','LineWidth',2,'DisplayName','P2');
hold on;
plot(t1,x1(6,:),'r','LineWidth',2,'DisplayName','P3');
hold off;
title('Phase 1')
xlabel('Time (min)');
ylabel('Protein concentration (nmol/gDW)');
xlim([0 660]);

x0=x1(:,end);
x2=Balances2(A,S,1,300,10,x0);
t2=661:961;
figure (2);
plot(t2,x2(4,:),'k','LineWidth',2,'DisplayName','P1');
hold on;
plot(t2,x2(5,:),'b','LineWidth',2,'DisplayName','P2');
hold on;
plot(t2,x2(6,:),'r','LineWidth',2,'DisplayName','P3');
hold off;
xlabel('Time (min)');
ylabel('Protein concentration (nmol/gDW)');
xlim([661 961]);
title('phase 2');
x=[x1 x2]; % put all concentration data into one matrix x

s1=perturb_phase1(A,S,x);
s2=perturb_phase2(A,S,x);

M1=trapz(s1,2)/600;
M2=trapz(s2,2)/750;

N1=zeros(3,6);
N1(:,1)=M1(1:3,1);
N1(:,2)=M1(4:6,1);
N1(:,3)=M1(7:9,1);
N1(:,4)=M1(10:12,1);
N1(:,5)=M1(13:15,1);
N1(:,6)=M1(16:18,1);

N2=zeros(3,6);
N2(:,1)=M2(1:3,1);
N2(:,2)=M2(4:6,1);
N2(:,3)=M2(7:9,1);
N2(:,4)=M2(10:12,1);
N2(:,5)=M2(13:15,1);
N2(:,6)=M2(16:18,1);

[U1,S1,V1]=svd(N1);
[U2,S2,V2]=svd(N2);