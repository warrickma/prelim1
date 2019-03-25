function x=Balances(A,S,tStart,tEnd)
pd=-log(0.5)/(60*24); % assume 24 hr for proein half life
L(1)=1200; L(2)=2400; L(3)=600;
xd=-log(0.5)/2.1;% mRNA degradation rate constant based on 2.1 min half life unit second ^(-1)
N_A=6.022*10^23;
p_conc=200/N_A*10^9/(2.8*10^(-13)*0.3); % gene p concentration in nmol/gDW
rnap_conc=1150/N_A*10^9/(2.8*10^(-13)*0.3); % rnap concentration in nmol/gDW
ribo_conc=45000/N_A*10^9/(2.8*10^(-13)*0.3); % ribosome concentration in nmol/gDW
e_x=60*60;ke1=e_x/L(1); % 60 nt/sec; compute for enlongation rate constant ke 
tau1=ke1/20;% tau is unitless;k_i (initiation rate constant) is from bionumber
ke2=e_x/L(2);
tau2=ke2/20;
ke3=e_x/L(3);
tau3=ke3/20;
K_x=0.24;%dissociation constant of RNAP in nmol/gDW
K_L=454.64;%dissociation constant of ribosome in nmol/gDW
n=1.5;Kd=0.3; % Kd (disocciation constant of inducer) is in mM
x=zeros(6,tEnd);
for t=tStart:tEnd
f1=0^n/(Kd^n+0^n);
f2=x(4,t)^n/(1000^n+x(4,t)^n);
f13=x(4,t)^n/(1000^n+x(4,t)^n);
f23=x(5,t)^10/(100000^10+x(5,t)^10);
u1=(0.0000001+100*f1)/(1+0.0000001+100*f1);
u2=(0.0000001+10*f2)/(1+0.0000001+10*f2);
u3=(0.0000001+5*f13+25*f23)/(1+0.0000001+5*f13+25*f23);
tx1=(rnap_conc*ke1/tau1*(p_conc/(p_conc+K_x)))*u1;
tx2=(rnap_conc*ke2/tau2*(p_conc/(p_conc+K_x)))*u2;
tx3=(rnap_conc*ke3/tau3*(p_conc/(p_conc+K_x)))*u3;
tl1=x(1,t)*ribo_conc/K_L*(16.5*60/L(1)*3);
tl2=x(2,t)*ribo_conc/K_L*(16.5*60/L(2)*3);
tl3=x(3,t)*ribo_conc/K_L*(16.5*60/L(3)*3);
r=[tx1;tx2;tx3;tl1;tl2;tl3];
x(:,t+1)=expm(A*1)*x(:,t)+inv(A)*[expm(A*1)-eye(6)]*S*r;
end
