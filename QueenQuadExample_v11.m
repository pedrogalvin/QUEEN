% This script provides a procedure that allows to compute the boundary integral equations by numerical quadratures.
% The methodology is based on the obtain of the weights of the quadrature rules by the solution of an underdetermined
% system of equations in the minimum norm sense.

%% Inputs
% ElementType=options;
% options:
% 'cheby1'                       % Chebyshev first kind (interior nodes)
% 'cheby2'                       % Chebyshev second kind
% 'lgl'                          % LGL integration points
% 'equi'                         % equidistant nodes
% ElementOrder=Integer;          % Element order
% iNode=Integer;                 % Collocation point
% Integrand
% options:
% f0= @(Xi) 1;                   % integralcase=0
% Singular functions
% f1= @(Xi) log(abs(Xi_p-Xi));   % integralcase=1
% f2 =@(Xi) 1./(Xi_p-Xi);        % integralcase=2
% f3 =@(Xi) 1./(Xi_p-Xi).^2;     % integralcase=3
%% Outputs
% ElementOrder                   % Element order
% ElementType                    % Element type
% ElementNodes                   % Nodal locations at natural coordinates
% TotalElementNodes              % Number of nodes
% ControlPoints                  % Control points used to define the Lagrange polynomial
% IntegrationPoints              % Integration points of the quadrature rules [8xTotalElementNodes]
% IntegrationWeights             % Integration weights of the quadrature rules [8xTotalElementNodes,TotalElementNodes]
% IntegrationMoments             % Generalised moments [3xTotalElementNodes,TotalElementNodes]
%% Syntax
% Element=QueenQuad(ElementOrder,ElementType);
%
%% Author Info
% Rocio Velazquez, Antonio Romero, Jose Dominguez, Antonio Tadeu, Pedro Galvin - Universidad de Sevilla (Spain), ITECONS (Portugal) - last modified 11/04/2022
%
%% References
% Rocio Velazquez, Antonio Romero, Jose Dominguez, Antonio Tadeu, Pedro Galvin,
% A novel high-performance quadrature rule for BEM formulations,
% Engineering Analysis with Boundary Elements, 2022, ISSN 0955-7997,
% http://dx.doi.org/j.enganabound.2022.04.036
%
% Rocio Velazquez, Antonio Romero, Antonio Tadeu, Pedro Galvin,
% Quadrature rule for BEM formulations involving singularities up to order 1/r2,

%
% You can download the QUEEN 1.1 source file (November 11th 2022) here:
% https://personal.us.es/pedrogalvin/queen.en.html

%% -------------------------------------

% Example
close all
ElementOrder=4;
ElementType='cheby1';
obj=QueenQuad_v11(ElementOrder,ElementType);
integralcase=3;

% Collocation point
iNode=3;

%%
% Natural coordinates of the collocation point
Xi_p=obj.ElementNodes(iNode);

% Integrand
% Reference
% http://dx.doi.org/j.enganabound.2022.04.036
f0= @(Xi) 1;                   %integralcase=0
% Singular functions
f1= @(Xi) log(abs(Xi_p-Xi));   %integralcase=1
f2 =@(Xi) 1./(Xi_p-Xi);        %integralcase=2
% Reference
% http://dx.doi.org/j.enganabound.
f3 =@(Xi) 1./(Xi_p-Xi).^2;     %integralcase=3

if integralcase==0
    Fun=@(Xi) bsxfun(@times,obj.ShapeFunction(obj,Xi),f0(Xi));
elseif integralcase==1
    Fun=@(Xi) bsxfun(@times,obj.ShapeFunction(obj,Xi),f1(Xi));
elseif integralcase==2
    Fun=@(Xi) bsxfun(@times,obj.ShapeFunction(obj,Xi),f2(Xi));
elseif integralcase==3
    Fun=@(Xi) bsxfun(@times,obj.ShapeFunction(obj,Xi),f3(Xi));
else
    error('integralcase variable should be 1, 2 or 3')
end

% Plot Shape Functions and Integrands
Xi=-1:1e-3:1;
figure(1)
box on
hold on
plot(Xi,obj.ShapeFunction(obj,Xi))
ylim([-0.5 1.5])
xlabel('$\xi$','Interpreter','latex')
ylabel('$\phi$','Interpreter','latex')
set(gca,'YScale','lin','YMinorTick','off','TickLabelInterpreter','latex')

figure(2)
box on
hold on
plot(Xi,Fun(Xi))
xlabel('$\xi$','Interpreter','latex')
ylabel('$\psi$','Interpreter','latex')
set(gca,'YScale','lin','YMinorTick','off','TickLabelInterpreter','latex')
if integralcase==0
    ylim([-0.5 1.5])
elseif integralcase==1
    ylim([-2 1])
elseif integralcase==2
    ylim([-10 10])
elseif integralcase==3
    ylim([-1 100])
end

% Integral computed by QueenQuad
tic;
Int_Q=Fun(obj.IntegrationPoints')*obj.IntegrationWeights(:,iNode);
elapsedTime = toc;
fprintf('Elapsed time (QueenQuad) %d seconds\n',elapsedTime)

% Exact value by moments
Int_M=obj.ControlPoints*obj.IntegrationMoments((1:obj.ElementOrder+1)+integralcase*(obj.ElementOrder+1),iNode);

% Integral computed by built-in integral function
tic;
Int_I=integral(Fun,-1,1,'ArrayValued',true);
elapsedTime = toc;
fprintf('Elapsed time (integral) %d seconds\n',elapsedTime)

% Results
for iorder=1:obj.ElementOrder+1
    fprintf('Shape function %d x Singular function. Singularity at node %d\n',iorder,iNode)
    fprintf('%s *** %s *** %s \n','Exact value','Queen','Integral')
    fprintf('%f *** %f *** %f \r\n',[Int_M(iorder);Int_Q(iorder);Int_I(iorder)])
end