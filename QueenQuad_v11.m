% QUEEN QUADRATURE (velazQez domingUez romEro tadEu galviN)
% Author Info
% Rocio Velazquez, Antonio Romero, Jose Dominguez, Antonio Tadeu, Pedro Galvin - Universidad de Sevilla (Spain), ITECONS (Portugal) - last modified 25/08/2022
%
% References

% Rocio Velazquez, Antonio Romero, Jose Dominguez, Antonio Tadeu, Pedro Galvin,
% A novel high-performance quadrature rule for BEM formulations,
% Engineering Analysis with Boundary Elements, 2022, ISSN 0955-7997,
% http://dx.doi.org/j.enganabound.2022.04.036

% Rocio Velazquez, Antonio Romero, Antonio Tadeu, Pedro Galvin,
% Quadrature rule for BEM formulations involving singularities up to order 1/r2,

%
% You can download the QUEEN 1.1 source file (November 11th 2022) here:
% https://personal.us.es/pedrogalvin/queen.en.html

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
% OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

classdef QueenQuad_v11
    %
    properties
        ElementOrder
        ElementType
        ElementNodes
        TotalElementNodes
        ControlPoints
        IntegrationPoints
        IntegrationWeights
        IntegrationMoments
    end
    %
    properties (Access=private)
        mQ=1    % parameter for singular integral
        nQ=4    % parameter for singular integral

    end
    %
    methods
        %--------------
        function [obj]=QueenQuad_v11(ElementOrder,ElementType)
            % constructor
            obj.ElementOrder=ElementOrder;
            obj.ElementType=ElementType;

            %
            obj=ComputeControlPoints(obj);
            obj=QuadratureRule(obj);

        end
        %--------------
        function TotalElementNodes=get.TotalElementNodes(obj)
            TotalElementNodes=obj.ElementOrder+1;
        end
        %--------------
        function ElementNodes=get.ElementNodes(obj)
            % This function define the element nodes

            N=obj.ElementOrder;

            if N==0
                ElementNodes=0;
            else

                % number of nodes
                j=0:N;

                switch obj.ElementType

                    case 'cheby1' % Chebyshev first kind (interior nodes)
                        ElementNodes=cos((2*j'+1)*pi/(2*N+2));
                    case 'cheby2' % Chebyshev second kind
                        ElementNodes=cos(pi*j'/N);
                    case 'lgl' % LGL integration points
                        ElementNodes=lglwt(N);
                    case 'equi' %  equidistant nodes
                        ElementNodes=linspace(-1,1,N+1)';
                end

                ElementNodes=sort(ElementNodes);

            end
        end
        %--------------
        function obj=ComputeControlPoints(obj)
            % This function return the control points of the Lagrange
            % polinomial in Bernstein basis

            N=obj.ElementOrder;

            % interpolation data
            Y=eye(N+1);

            % nodes [0,1]
            X=obj.Coordinate_t(obj.ElementNodes);

            % control point computation
            C=zeros(N+1);

            for l=1:N+1

                % interpolation data
                f=Y(:,l);

                % divided differences
                F=zeros(N+1);

                % control points
                w=zeros(N+1);
                c=zeros(N+1);

                for k=0:N
                    % divided differences
                    for j=k:N
                        if k==0
                            F(:,k+1)=f;
                        else
                            F(j+1,k+1)=(F(j+1,k)-F(j,k))/(X(j+1)-X(j+1-k));
                        end
                    end
                    % control points
                    for j=0:k
                        if k==0
                            w(k+1,j+1)=1;
                            c(k+1,j+1)=f(1);
                        else
                            if j==0
                                w(k+1,j+1)=-(1-j/k)*w(k,j+1)*X(k);
                                c(k+1,j+1)=(1-j/k)*c(k,j+1)+w(k+1,j+1)*F(k+1,k+1);
                            else
                                w(k+1,j+1)=j/k*w(k,j)*(1-X(k))-(1-j/k)*w(k,j+1)*X(k);
                                c(k+1,j+1)=j/k*c(k,j)+(1-j/k)*c(k,j+1)+w(k+1,j+1)*F(k+1,k+1);
                            end
                        end
                    end
                end
                C(l,:)=c(k+1,:);
            end

            obj.ControlPoints=C;
        end
        %--------------
        function obj=QuadratureRule(obj)

            % This subroutine computes N integration weights for
            % integration accuracy of M-1, when the collocation point is given by the
            % natural coordinate x.

            N=obj.ElementOrder;

            % Legendre-Gauss nodes
            [Xi,~]=lgwt(obj.nQ*(N+1),-1,1);
            t=obj.Coordinate_t(Xi);

            % weights and moments
            obj.IntegrationWeights=zeros(obj.nQ*(N+1),obj.TotalElementNodes);
            obj.IntegrationMoments=zeros(4*(N+1),obj.TotalElementNodes); 

            for iNode=1:obj.TotalElementNodes

                % Singular point
                Xi_p=obj.ElementNodes(iNode);
                t_p=obj.Coordinate_t(Xi_p);

                % Singular functions
                f1=@(t) log(abs(Xi_p-(2*t-1)));
                f2=@(t) 1./(Xi_p-(2*t-1));
                f3=@(t) 1./(Xi_p-(2*t-1)).^2;

                % Compute function phi at each t
                P1=zeros(N+1,obj.nQ*(N+1));
                P2=zeros(N+1,obj.nQ*(N+1));
                P3=zeros(N+1,obj.nQ*(N+1));
                P4=zeros(N+1,obj.nQ*(N+1));
                %
                for k=0:N
                    P1(k+1,:)=obj.Bern(N,k,t)*2;
                    P2(k+1,:)=obj.Bern(N,k,t).*f1(t)*2;
                    P3(k+1,:)=obj.Bern(N,k,t).*f2(t)*2;
                    P4(k+1,:)=obj.Bern(N,k,t).*f3(t)*2;
                end
                Phi=[P1; P2; P3; P4];

                % Moments computation
                Int1=zeros(N+1,1);
                Int2=zeros(N+1,1);
                Int3=zeros(N+1,1);
                Int4=zeros(N+1,1);
                %
                for k=0:N
                    Int1(k+1)=obj.BernIntegral(N);

                    Int2(k+1)=obj.LogBernIntegral(N,k,Xi_p);

                    Int3(k+1)=obj.Bern(N,k,t_p)*obj.FinitePart2a(Xi_p) +...
                        integral(@(t) ((obj.Bern(N,k,t)-obj.Bern(N,k,t_p)).*f2(t))*2,0,1);

                    Int4(k+1,:)=obj.Bern(N,k,t_p)*obj.FinitePart2b(Xi_p) - ...
                        0.5*obj.BernDerivate(N,k,t_p)*obj.FinitePart2a(Xi_p)+...
                        quadgk(@(t) ((obj.Bern(N,k,t)-obj.Bern(N,k,t_p)+0.5*obj.BernDerivate(N,k,t_p).*(Xi_p-(2*t-1))).*f3(t))*2,0,1);   %integral

                end
                obj.IntegrationMoments(:,iNode)=[Int1; Int2; Int3; Int4];

                % Solve
                obj.IntegrationWeights(:,iNode)=...
                    2*lsqminnorm(Phi,obj.IntegrationMoments(:,iNode));
            end
            obj.IntegrationPoints=Xi;
        end
        %--------------
    end


    methods (Static)
        %--------------
        function [t]=Coordinate_t(xi)
            t=0.5*(xi+1);
        end
        %--------------
        function P=ShapeFunction(obj,Xi,k)
            % This function computes the element shape function at Xi
            % points. Only return the k-th function if the input argument k
            % is given.

            N=obj.ElementOrder;
            %
            nXi=length(Xi);
            t=obj.Coordinate_t(Xi);
            %
            A=zeros(N+1,nXi);
            for i=0:N
                A(i+1,:)=obj.Bern(N,i,t);
            end
            P=obj.ControlPoints*A;
            %
            if nargin==3
                P=P(k,:);
            end
        end
        %--------------
        function B=Bern(n,k,t)
            a=0;
            b=1;
            bin=factorial(n)/(factorial(n-k)*factorial(k));
            B=(b-a)^-n*bin*(t-a).^k.*(b-t).^(n-k);
        end
        %--------------
        function Dif=BernDerivate(n,k,t)
            % This function computes the Bernstein polinomial derivate of
            % order n
            a=0;
            b=1;
            if n==0
                Dif=0;
            else
                if k==0
                    Dif=n/(b-a)*(-Bern(n-1,k,t));
                elseif n==k
                    Dif=n/(b-a)*(Bern(n-1,k-1,t));
                else
                    Dif=n/(b-a)*(Bern(n-1,k-1,t)-Bern(n-1,k,t));
                end
            end
        end
        %--------------
        function Int=BernIntegral(n)
            % This function computes the Bernstein polinomial integral of
            % order n from -1 to +1
            a=0;
            b=1;
            Int=(b-a)*(n+1)^-1*2;
        end
        %--------------
        function Int=LogBernIntegral(n,k,Xi_p)
            % This function computes the integral B_k^n*log(abs(x-t)
            fun=@(t) Bern(n,k,t).*log(abs(Xi_p-(2*t-1)))*2;
            Int=quadgk(fun,0,1);
        end
        %--------------
        function Int=FinitePart2a(x)
            % This function computes the finite part of integral 1/(x-t)
            if abs(x)==1
                Int=x*log(2);
            else
                Int=log(abs((x+1)./(1-x)));
            end
        end
        %--------------
        function Int2b=FinitePart2b(x)
            % This function computes the finite part of integral 1/(x-t)^2
            if abs(x)==1
                Int2b=-0.5;
            else
                Int2b=2./(x.^2-1);
            end
        end
        %--------------
    end
    %
end

%--------------
function B=Bern(n,k,t)
bin=factorial(n)/(factorial(n-k)*factorial(k));
B=bin*t.^k.*(1-t).^(n-k);
end

%--------------

function [x,w,P]=lglwt(N)

% Copyright (c) 2009, Greg von Winckel
% All rights reserved.
% Downloaded from https://es.mathworks.com/matlabcentral/fileexchange/4775-legende-gauss-lobatto-nodes-and-weights

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% The downloaded function has been modified

% Truncation + 1
N1=N+1;

% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x=-cos(pi*(0:N)/N)';

% The Legendre Vandermonde Matrix
P=zeros(N1,N1);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and
% update x using the Newton-Raphson method.

xold=2;

while max(abs(x-xold))>eps

    xold=x;

    P(:,1)=1;    P(:,2)=x;

    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end

    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );

end

w=2./(N*N1*P(:,N1).^2);
end
%--------------

function [x,w]=lgwt(N,a,b)

% Copyright (c) 2009, Greg von Winckel
% All rights reserved.
% Downloaded from https://es.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% The downloaded function has been modified

N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps


    L(:,1)=1;
    Lp(:,1)=0;

    L(:,2)=y;
    Lp(:,2)=1;

    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end

    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);

    y0=y;
    y=y0-L(:,N2)./Lp;

end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;

% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end
%--------------