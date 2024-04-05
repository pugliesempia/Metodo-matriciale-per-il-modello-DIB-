%a%    Studio i tempi computazionali ed anche la convergenza
       
       
%b%    Stampando le soluzioni per i casi: 
       % Nx=50 Lx=10 ottengo un pattern intrinseco 
       % Nx=200 Lx=100 coglierà la peculiarità del modello
       % Nx=50 Lx=100 pattern spurio
 
%Script per la risoluzione del modello DIB 
par=struct();
par.d1=1; par.d2=20;
par.ro=1;
par.A1=10; par.A2=30;
par.k2=2.5; par.k3=1.5;  par.a=0.5; par.gamma=0.2; par.D=2.4545;
par.B=85; par.C=3;

F1= @(u, v, par) par.ro*(par.A1*(1-v).*u-par.A2*u.^3- par.B*(v-par.a));
F2= @(u, v, par) par.ro*(par.C*(1+par.k2*u).*(1-v).*(1-par.gamma*(1-v))-par.D*v.*(1+par.k3*u).*(1+par.gamma*v));

theta=1/2; %Theta varia tra [0.5,1](cambiano i tempi di convergenza)
ue=0;ve=0.5;

N=200;Lx=100;
Tf=100; ht=1e-2;
hs =Lx/(N+1);
Nt=fix(Tf/ht);
d=20;
x=linspace(0,Lx,N+2); y=x; 
xi=x(2:end-1); yi=xi; 
%[X,Y]=meshgrid(x,y); 
[Xi,Yi]=meshgrid(xi,yi);  

uno=ones(N,1);
A2=spdiags([uno,-2*uno,uno],-1:1,N,N);
B_nc=spalloc(N,N,4);
B_nc(1,1)=2;     B_nc(end,end)=2;
B_nc(1,2)=-1/2;  B_nc(end,end-1)=-1/2; 
B_nc=2/3*B_nc;
T=A2+B_nc;
I=speye(N);
T=(1/hs^2)*T;
tic

%ho due eq. di sylvester: A1 Un+1 +B Un+1= C1
%                         A2 Vn+1 +B Vn+1= C2             

A1=I-ht*theta*T; 
A2=I-d*ht*theta*T;

%calcolo gli autovettori, gli autovettori e l'inversa
[Qa1,Ra1]=eig(full(A1));   I_Qa1=inv(Qa1);       %coeff A1 della eq di syvester in U
[Qa2,Ra2]=eig(full(A2));   I_Qa2=inv(Qa2);       %coeff della eq di syvester in V
B=-ht*theta*T'; 
[Qb,Rb]=eig(full(B));      I_Qb=inv(Qb);

%estraggo gli autovalori
dA1=diag(Ra1); dA2=diag(Ra2); dB=diag(Rb);

%Calcolo Lij cappuccio
Inveig_u=1./(dA1*ones(1,N)+ones(N,1)*dB');
Inveig_v=1./(dA2*ones(1,N)+ones(N,1)*d*dB');

rng('default');
Rs=rand(N); 
U=ue+Rs*1e-05; 
V=ve+Rs*1e-05;
Incremento=zeros(Nt,1);
Media=zeros(Nt,1);
tempi_array=zeros(Nt,1);
%Ciclo temporale
for n=1:Nt
    C_u=(I+ht*(1-theta)*T)*U+ U*(ht*(1-theta)*T')+ht* F1(U,V,par);
    C_v=(I+d*ht*(1-theta)*T)*V+ V*(d*ht*(1-theta)*T')+ ht* F2(U,V,par);

    Cc_u= I_Qa1*C_u*Qb;
    Cc_v= I_Qa2*C_v*Qb;
    
    Xx= Cc_u.*Inveig_u;
    Yy= Cc_v.*Inveig_v;

    U_new= Qa1*Xx*I_Qb;
    V_new= Qa2*Yy*I_Qb;
    Incremento(n)=norm(U_new-U,'fro');
    Media(n)=hs^2*sum(sum(U))/(Lx^2);
     
    U= U_new; 
    V= V_new;
    tempi_array(n)=ht*n;
   
end
tempo=toc;

% Visualizzazione dei risultati
  figure, surf(Xi,Yi,U,'edgecolor','none','facecolor','interp'), view(2),colorbar
  figure, surf(Xi,Yi,V,'edgecolor','none','facecolor','interp'), view(2),colorbar
  
  figure
  subplot(121),semilogy(tempi_array,Incremento)
  subplot(122),plot(tempi_array,Media)
  
  
  