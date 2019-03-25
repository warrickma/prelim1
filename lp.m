xd=8.35;% mRNA degradation rate constant based on 2.1 min half life unit h ^(-1)
N_A=6.022*10^23;
pd=9.9*10^(-3); %h^(-1)
p_conc=5*10^(-3); % gene p concentration in microM
rnap_conc=0.15; % rnap concentration in microM
ribo_conc=1.6; % ribosome concentration in microM
e_x=60*3600;ke1=e_x/924; % 60 nt/sec; compute for enlongation rate constant ke 
tau1=2.7;% time constant for transcription
tau2=0.8;%time constant for translation
K_x=0.3;%dissociation constant of RNAP in microM
K_L=57;%dissociation constant of ribosome in microM
r_x=rnap_conc*ke1/tau1*(p_conc/(p_conc+K_x));

syms v1 v2 v3 v4 v5 v6 b1 b2 b3 b4 b5 b6 b7 b8 b9
r={v1;v2;v3;v4;v5;v6;b1;b2;b3;b4;b5;b6;b7;b8; b9};
S1=sym(S);
sv=S1*r;
lb=[0,0,0,0,0,0,-100000,-100000,-100000,-100000,-100000,-100000,-100000,-100000,-100000];
f=[0 0.52 0 0 0 0 0 0 0 0 0 0 0 0 0];
v=zeros(1000,1);
I=zeros(1000,1);
for i=1:1000
    f1=(i/100)^1.5/(0.3^1.5+(i/100)^1.5);
    u1=(0.26+300*f1)/(1+0.26+300*f1);
    mRNA=r_x*u1/xd;
    
up=16.5*3600/300*1.6*0.8*(mRNA/(57+mRNA)); % set upper bound
ub=[1000,50,10,1000,up,10,100000,100000,100000,100000,100000,100000,100000,100000,100000];
Aeq=[-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 -924 0 0 0 0 0 -1 0 0 0 0 0 0 0;
    0 0 924 0 0 0 0 0 0 1 0 0 0 0 0;
    -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 1 -1 -1 1 0 0 0 0 0 0 0 0 0 0;
    0 1 -1 -1 1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 308 -1 0 0 0 0 0 0 0 0 0;
    0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 -616 0 0 0 0 0 0 0 -1 0 0;
    0 0 0 0 -308 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 1 0 0 0 0 0 0;
    0 0 0 0 616 0 0 0 0 0 0 0 0 -1 0;
    0 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0;
    0 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 1 0 0 0 0 0 1 0 0 0]; % sv=Aeq
x=linprog(-f,[],[],Aeq,zeros(1,17),lb,ub);
v(i)=x(2);
I(i)=i/100;
end
semilogx(I,v,'LineWidth',2,'Color',[0 0 0]);
xlabel('I');
ylabel('Protein Steady State Concentration (nM)');