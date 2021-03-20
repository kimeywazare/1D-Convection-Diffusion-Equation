function [K,f] = Galerkin_system(X,T,referenceElement,example)
%
% [K,f] = Galerkin_system(X,referenceElement,example)
% System obtained by discretizing the weak form associated to
% the convection-diffusion equation
%               a ux - nu uxx = f
% using the Galerkin formulation and linear interpolation.
% Boundary conditions are not considered


% reference element information
nen = referenceElement.nen; 
ngaus = referenceElement.ngaus; 
wgp = referenceElement.GaussWeigths; 
N = referenceElement.N; 
Nxi = referenceElement.Nxi; 

% example properties
a = example.a; 
nu = example.nu; 

% Number of nodes and elements
nPt = length(X); 
nElem = size(T,1); 

K = zeros(nPt,nPt);
f = zeros(nPt,1);

% Loop on elements
for ielem=1:nElem
    Te = T(ielem,:); 
    Xe = X(Te); 
    h = Xe(end) - Xe(1);
    Ke = zeros(nen); 
    fe = zeros(nen,1); 
    % Loop on Gauss points
    for ig = 1:ngaus
        N_ig = N(ig,:);
        Nx_ig = Nxi(ig,:)*2/h;
        w_ig = wgp(ig)*h/2;        
        Ke = Ke + w_ig*(N_ig'*(a*Nx_ig) + Nx_ig'*(nu*Nx_ig));
        x = N_ig*Xe; % x-coordinate of the gauss point
        s = SourceTerm(x,example); 
        fe = fe + w_ig*(N_ig')*s;
    end
    % Assmebly
    K(Te,Te) = K(Te,Te) + Ke; 
    f(Te) = f(Te) + fe; 
end