function s1=perturb_phase1(A,S,x)
% Phase 1

[p1,p2,p3]=perturb1(A,S,1,650,10,1,1,1,1,1); % rnap conc * 10
%calculate partial derivative
N_A=6.022*10^23;

deltap_a1=(p1-x(4,550:600))./(1150/N_A*10^9/(2.8*10^(-13)*0.3)*9);%calculate partial derivative
sc_a1=1150/N_A*10^9/(2.8*10^(-13)*0.3)./x(4,550:600); % calculate scale factor
s_a1=sc_a1.*deltap_a1;% compute sensitivity

deltap_a2=(p2-x(5,550:600))./(1150/N_A*10^9/(2.8*10^(-13)*0.3)*9);%calculate partial derivative
sc_a2=1150/N_A*10^9/(2.8*10^(-13)*0.3)./x(5,550:600); % calculate scale factor
s_a2=sc_a2.*deltap_a2;% compute sensitivity

deltap_a3=(p3-x(6,550:600))./(1150/N_A*10^9/(2.8*10^(-13)*0.3)*9);%calculate partial derivative
sc_a3=1150/N_A*10^9/(2.8*10^(-13)*0.3)./x(6,550:600); % calculate scale factor
s_a3=sc_a3.*deltap_a3;% compute sensitivity

% perturb b ribosome conc 10 folds increase
[p1,p2,p3]=perturb1(A,S,1,650,1,10,1,1,1,1);

deltap_b1=(p1-x(4,550:600))./(45000/N_A*10^9/(2.8*10^(-13)*0.3)*9);%calculate partial derivative
sc_b1=45000/N_A*10^9/(2.8*10^(-13)*0.3)./x(4,550:600); % calculate scale factor
s_b1=sc_b1.*deltap_b1;% compute sensitivity

deltap_b2=(p2-x(5,550:600))./(45000/N_A*10^9/(2.8*10^(-13)*0.3)*9);%calculate partial derivative
sc_b2=45000/N_A*10^9/(2.8*10^(-13)*0.3)./x(5,550:600); % calculate scale factor
s_b2=sc_b2.*deltap_b2;% compute sensitivity

deltap_b3=(p3-x(6,550:600))./(45000/N_A*10^9/(2.8*10^(-13)*0.3)*9);%calculate partial derivative
sc_b3=45000/N_A*10^9/(2.8*10^(-13)*0.3)./x(6,550:600); % calculate scale factor

s_b3=sc_b3.*deltap_b3;% compute sensitivity

% perturb c transcription elongation rate increases 10 folds
[p1,p2,p3]=perturb1(A,S,1,650,1,1,10,1,1,1);

deltap_c1=(p1-x(4,550:600))./(3600*9);%calculate partial derivative
sc_c1=3600./x(4,550:600); % calculate scale factor
s_c1=sc_c1.*deltap_c1;% compute sensitivity

deltap_c2=(p2-x(5,550:600))./(3600*9);%calculate partial derivative
sc_c2=3600./x(5,550:600); % calculate scale factor
s_c2=sc_c2.*deltap_c2;% compute sensitivity

deltap_c3=(p3-x(6,550:600))./(3600*9);%calculate partial derivative
sc_c3=3600./x(6,550:600); % calculate scale factor

s_c3=sc_c3.*deltap_c3;% compute sensitivity

% perturb d translation elongation rate increases 10 folds
[p1,p2,p3]=perturb1(A,S,1,650,1,1,10,1,1,1);

deltap_d1=(p1-x(4,550:600))./(16.5*60*9);%calculate partial derivative
sc_d1=16.5*60./x(4,550:600); % calculate scale factor
s_d1=sc_d1.*deltap_d1;% compute sensitivity

deltap_d2=(p2-x(5,550:600))./(16.5*60*9);%calculate partial derivative
sc_d2=16.5*60./x(5,550:600); % calculate scale factor
s_d2=sc_d2.*deltap_d2;% compute sensitivity

deltap_d3=(p3-x(6,550:600))./(16.5*60*9);%calculate partial derivative
sc_d3=16.5*60./x(6,550:600); % calculate scale factor

s_d3=sc_d3.*deltap_d3;% compute sensitivity

% perturb d translation elongation rate increases 10 folds
[p1,p2,p3]=perturb1(A,S,1,650,1,1,10,1,1,1);

deltap_d1=(p1-x(4,550:600))./(16.5*60*9);%calculate partial derivative
sc_d1=16.5*60./x(4,550:600); % calculate scale factor
s_d1=sc_d1.*deltap_d1;% compute sensitivity

deltap_d2=(p2-x(5,550:600))./(16.5*60*9);%calculate partial derivative
sc_d2=16.5*60./x(5,550:600); % calculate scale factor
s_d2=sc_d2.*deltap_d2;% compute sensitivity

deltap_d3=(p3-x(6,550:600))./(16.5*60*9);%calculate partial derivative
sc_d3=16.5*60./x(6,550:600); % calculate scale factor

s_d3=sc_d3.*deltap_d3;% compute sensitivity% perturb d translation elongation rate increases 10 folds
[p1,p2,p3]=perturb1(A,S,1,650,1,1,10,1,1,1);

% perturb e Kx increases 10 folds
deltap_e1=(p1-x(4,550:600))./(0.24*9);%calculate partial derivative
sc_e1=0.24./x(4,550:600); % calculate scale factor
s_e1=sc_e1.*deltap_e1;% compute sensitivity

deltap_e2=(p2-x(5,550:600))./(0.24*9);%calculate partial derivative
sc_e2=0.24./x(5,550:600); % calculate scale factor
s_e2=sc_e2.*deltap_e2;% compute sensitivity

deltap_e3=(p3-x(6,550:600))./(0.24*9);%calculate partial derivative
sc_e3=0.24./x(6,550:600); % calculate scale factor

s_e3=sc_e3.*deltap_e3;% compute sensitivity

% perturb f KL increases 10 folds
deltap_f1=(p1-x(4,550:600))./(454.64*9);%calculate partial derivative
sc_f1=454.64./x(4,550:600); % calculate scale factor
s_f1=sc_f1.*deltap_f1;% compute sensitivity

deltap_f2=(p2-x(5,550:600))./(454.64*9);%calculate partial derivative
sc_f2=454.64./x(5,550:600); % calculate scale factor
s_f2=sc_f2.*deltap_f2;% compute sensitivity

deltap_f3=(p3-x(6,550:600))./(454.64*9);%calculate partial derivative
sc_f3=454.64./x(6,550:600); % calculate scale factor

s_f3=sc_f3.*deltap_f3;% compute sensitivity

s1=vertcat(s_a1,s_a2,s_a3,s_b1,s_b2,s_b3,s_c1,s_c2,s_c3,s_d1,s_d2,s_d3,s_e1,s_e2,s_e3,s_f1,s_f2,s_f3);
 
end