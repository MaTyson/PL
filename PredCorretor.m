%Método de pontos interiores Preditor-Corretor%
function [x,y,z]=PredCorretor(A,b,c)
%Inicialização%
[m, n]=size(A);
xtil=A'*((A*A')\b);
e2=100;
e1=min([-min(xtil),e2,norm(b,1)/(e2*norm(A,1))]);
e3=1+norm(c,1);

%Ponto inicial%
x=max(xtil,e1);
y=zeros(m,1);
z=zeros(n,1);
z(c>=0)=c(c>=0)+e3;
z(c<=-e3)=-c(c<=-e3);
aux=min(c,0);
z((aux>-e3))=e3;
clear aux;
%fim inicialização%

%parâmetros%
tau=0.99995;
eps=1e-8;
e=ones(n,1);
k=0;
itmax=1e5;

%loop até convergência ou limite de iterações itmax%
while k<=itmax
    
    %rp, rd e ra%
    rp=b-A*x;
    rd=c-A'*y-z;
    ra=-diag(sparse(x))*diag(sparse(z))*e;
    
    %definimos D%
    D=diag(sparse(z.^-1))*diag(sparse(x));
    
    %fazemos a fatoração de cholesky de ADA'%
    B=sparse(A*D*A');
    [R, ~]=chol(B);
    
    %primeiro cálculo do d%
    dy=R\(R'\(rp+A*D*(rd-diag(sparse(x.^-1))*ra)));
    dx=D*(A'*dy-rd+diag(sparse(x.^-1))*ra);
    dz=diag(sparse(x.^-1))*(ra-diag(sparse(z))*dx);
    
    %rhop e rhod%
    rhop=inf;
    rhod=inf;
    for i=1:n
        if(dx(i)<0)
            rhop=min(rhop,-x(i)/dx(i));
        end
        if(dz(i)<0)
            rhod=min(rhod,-z(i)/dz(i));
        end
    end
    
    %alphap e alphad%
    alphap=min(1,(tau*rhop));
    alphad=min(1,(tau*rhod));
    
    %gap gamma e gammatil$
    gama=z'*x;
    gamatil=(x+alphap*dx)'*(z+alphad*dz);
    
    %diferentes valores de mi%
    %mi=((gamatil/gama)^2)*(gamatil/n);
    %mi=gama/n^2;
    if(gama<1)
        mi=(gama/n)^2;
    else
        mi=((gamatil/gama)^2)*(gamatil/n); 
    end
    %devemos escolher manualmente qual mi usar no calculo de rs
    rs= mi*e + ra - diag(sparse(dx))*diag(sparse(dz))*e;
    
    %cálculo do d%
    dy=R\(R'\(rp+A*D*(rd-diag(sparse(x.^-1))*rs)));
    dx=D*(A'*dy-rd+diag(sparse(x.^-1))*rs);
    dz=diag(sparse(x.^-1))*(rs-diag(sparse(z))*dx);
    
    %rhop e rhod%
    rhop=inf;
    rhod=inf;
    for i=1:n
        if(dx(i)<0)
            rhop=min(rhop,-x(i)/dx(i));
        end
        if(dz(i)<0)
            rhod=min(rhod,-z(i)/dz(i));
        end
    end
    
    %alphap e alphad%
    alphap=min(1,(tau*rhop));
    alphad=min(1,(tau*rhod));
    
    %atualização%
    x=x+alphap*dx;
    y=y+alphad*dy;
    z=z+alphad*dz;
    
    %iterações%
    k=k+1;
    
    %otimalidade%
    Fp=norm((b-A*x),1)/(1+norm(b,1));
    Fd=norm((c-A'*y-z),1)/(1+norm(c,1));
    Fa=gama/(1+norm(c'*x,1)+norm(b'*y,1));
    
    %critério de parada%
    if((Fp<eps)&&(Fd<eps)&&(Fa<eps))
        break
    end    
     
    %gap por iteração%
    fprintf(1,' \n');
    fprintf(1,'valor do gap: %f \n',gama);
end

%impressão dos resultados%
fprintf(1,' \n');
fprintf(1,'valor do gap: %f \n',gama);
fprintf(1,' \n');
fprintf(1,'Mínimo da função objetivo Primal (c^t*x): %f \n',c'*x);
fprintf(1,' \n');
fprintf(1,'Iterações Realizadas: %d \n',k);
fprintf(1,' \n');
end