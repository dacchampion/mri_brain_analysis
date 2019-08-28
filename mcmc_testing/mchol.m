function [L,D,E,pneg]=mchol(G)
%
%  n gives the size of the matrix.
%
n=size(G,1);
%
%  gamma, zi, nu, and beta2 are quantities used by the algorithm.  
%
gamma=max(diag(G));
zi=max(max(G-diag(diag(G))));
nu=max([1,sqrt(n^2-1)]);
beta2=max([gamma, zi/nu, 1.0E-15]);
%
%  Initialize diag(C) to diag(G).
%
C=diag(diag(G));
%
%  Loop through, calculating column j of L for j=1:n
%

L=zeros(n);
D=zeros(n);
E=zeros(n);

for j=1:n,
    bb=[1:j-1];
    ee=[j+1:n];

    %
    %  Calculate the jth row of L.  
    %
    if (j > 1),
        L(j,bb)=C(j,bb)./diag(D(bb,bb))';
    end;
    %
    %  Update the jth column of C.
    %
    if (j >= 2),
        if (j < n), 
            C(ee,j)=G(ee,j)-(L(j,bb)*C(ee,bb)')';
        end;
    else
        C(ee,j)=G(ee,j);
    end;
    %
    % Update theta. 
    %
    if (j == n)
        theta(j)=0;
    else
        theta(j)=max(abs(C(ee,j)));
    end;
    %
    %  Update D
    %
    D(j,j)=max([eps,abs(C(j,j)),theta(j)^2/beta2]');
    %
    % Update E.
    %
    E(j,j)=D(j,j)-C(j,j);

    
    %
    %  Update C again...
    % 
    %%%%%%%% M.Zibulevsky: begin of changes, old version is commented %%%%%%%%%%%%%
    
    %for i=j+1:n,
    %    C(i,i)=C(i,i)-C(i,j)^2/D(j,j);
    %end;
    
    ind=[j*(n+1)+1 : n+1 : n*n]';
    C(ind)=C(ind)-(1/D(j,j))*C(ee,j).^2;


end;

%
% Put 1's on the diagonal of L
%
%for j=1:n,
%    L(j,j)=1;
%end;

ind=[1 : n+1 : n*n]';
L(ind)=1;

%%%%%%%%%%%%%%%%%%%%%%%% M.Zibulevsky: end of changes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  if needed, find a descent direction.  
%
if ((nargout == 4) & (min(diag(C)) < 0.0))
    [m,col]=min(diag(C));
    rhs=zeros(n,1);
    rhs(col)=1;
    pneg=L'\rhs;
else
  pneg=[];
end;


return