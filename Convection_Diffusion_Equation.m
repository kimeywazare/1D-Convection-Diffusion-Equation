% Convection Diffusion Equation
% a*du/dx + nu*d^2u/dx^2 = 0 in [0,L]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc ; close all; clear all;
%-------------------------------------------------------------------------%
                        %PRE-PROCESS%
%-------------------------------------------------------------------------%
% Define Input Variables
%Change a and nu for Peclet No.[PE - 0,0.5,1,,10,50,500,1000,1000] 
a = 2;                           % Convection Coefficient 
nu = 0.0001;                       % Diffusion Coefficient
dom = [0,1];                     % Computational domain
example.dom = dom; 
nElem = 10;                      % Number of Elements
nPt = nElem + 1;                 % Number of Nodes
%Domain Discretization
h = (dom(2) - dom(1))/nElem;     %h = 0.1
X = (dom(1):h:dom(2))'; 
T = [1:nPt-1; 2:nPt]';  
x=X;L=1;
% Reference element: numerical quadrature and shape functions
p = 1;                                  %gauss point
referenceElement = SetRefereceElement(p); 
example.a = a; 
example.nu = nu; 
Pe = a*h/(2*nu);
%disp(' ')
disp(['Peclet number: ',num2str(Pe)]);
%-------------------------------------------------------------------------%
                    %MAIN (SIMULATION)%
%-------------------------------------------------------------------------%
%---------------------%Exact Solution%------------------------------------%
u_exact = (1-(exp((a*(x-L)/nu))))/(1-exp(-(a*L)/nu));

%---------------------%GALERKIN (G) METHOD%-------------------------------%
disp('Galerkin formulation')
example.methodName = 'Galerkin'; 
example.tau = 0; 
[K,F] = Galerkin_system(X,T,referenceElement,example); 

%Lagrangian Multiplier
A = zeros(2,nPt);
A(1,1) = 1; A(2,nPt) = 1;
b = zeros(nPt,1);          
bd = [1;0];

Ad = [K A';A zeros(2,2)];
B = [F;bd];
sol = Ad\B;
u_g= sol(1:nPt); %Solution of Galerkin Method

%---------------------%Full Upwind (UP) METHOD%---------------------------%
disp('Full-Upwind formulation')
example.methodName = 'Full Upwind'; 
alpha = 1;
tau= alpha*h/(2*a);
example.tau = cinput ('Stabilization parameter to be used',tau);
[K,F] = Upwind_system(X,T,referenceElement,example); 
      
%Lagrangian Multiplier
A = zeros(2,nPt);
A(1,1) = 1; A(2,nPt) = 1;
b = zeros(nPt,1);          
bd = [1;0];

Ad = [K A';A zeros(2,2)];
B = [F;bd];
sol = Ad\B;
u_up= sol(1:nPt); %Solution of Full-Upwind Method

%-------------------%Petro-Galerkin (PG) METHOD%--------------------------%
disp('Petro-Galerkin formulation')
example.methodName = 'PG'; 
alpha = coth(Pe)-1/Pe;                %alpha=1 recover Full-Upwind
tau = alpha*h/(2*a);
example.tau = cinput ('Stabilization parameter to be used',tau);
[K,F] = PG_system(X,T,referenceElement,example);
    
%Lagrangian Multiplier
A = zeros(2,nPt);
A(1,1) = 1; A(2,nPt) = 1;
b = zeros(nPt,1);          
bd = [1;0];

Ad = [K A';A zeros(2,2)];
B = [F;bd];
sol = Ad\B;
u_pg= sol(1:nPt); %Solution of Petro-Galerkin Method
%-------------------------------------------------------------------------%
                       %POST-PROCESS%
%-------------------------------------------------------------------------%
%Plot graphs of methods for differnet peclet number. 
plot(X,u_g,'b-*',X,u_up,'g-d',X,u_pg,'r-h', X,u_exact,'k--','LineWidth',1,'MarkerSize',5)
%plot(X,u_up,'g-d',X,u_pg,'r-h', X,u_exact,'k--','LineWidth',2,'MarkerSize',5)
l = legend('Galarkin','Upwind', 'PG','Exact'); 
%l = legend('Upwind', 'PG','Exact'); 
set(l,'FontSize',18);set(gca, 'FontSize',18);xlim([0 1]); grid on;
title(['Pe = ',num2str(Pe)],'fontsize',20);
%title(['Pe = 0'],'fontsize',20);
xlabel('Domain (X)','fontsize',18); 
ylabel('Scalar Variable (u)','fontsize',18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
