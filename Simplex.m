%min f(x)=c'*x
%s.a: Ax=b, x>=0
clear
clc
fprintf('Simplex\nEntrada de dados\n')
%inicializamos a contagem k e o Big M
k=0;
M=1e5;
%pedimos ao usuario as dimensoes de A, o vetor dos custos, 
%o vetor dos recursos e a propria A
m=input('Digite o numero de linhas de A:')
n=input('Digite o numero de colunas de A:')
c=zeros(n,1);
c=input('Digite o vetor dos custos:' )
b=input('Digite o vetor dos recursos:')
A=input('Digite a matriz A:')
%inicializamos a matriz B como Identidade de dimensão m
B=eye(m);
%fazemos N ser A
N=A;
%fazemos A=[B N]
A=[B N];
%inicializamos o vetor dos indices basicos e nao basicos
indb=linspace(1,m,m);
indn=linspace(m+1,m+n,n);
%criamos m variaveis artificiais de custo M
cb=M*ones(m,1)
%fazemos o custo nao basico o proprio vetor dos custos inicial
cn=c
%procuramos solucao basica factivel Bxb=b
xn=zeros(n,1);
xb=B\b
%se alguma componente de xb for negativa, retornamos como 
%problema infactivel
for i=1:m
  if(xb(i)<0)
    fprintf('Problema infactivel\n')
    return
  endif
endfor
%calculamos f de xb e imprimimos seu valor
f=cb'*xb;
fprintf('o valor de f e %f\n',f)

%calculamos o vetor multiplicador simplex B'lambda=cb
lambda=B'\cb

%calculamos os custos relativos e os colocamos num vetor r
r=zeros(n,1);
for i=1:n
  r(i)=cn(i)-lambda'*A(:,indn(i));
endfor
r
%achamos o indice da variavel que entra na base e colocamos em entra
%pegamos o primeiro custo relativo que for o menor em caso de 
%mais de um menor valor
entra=find(r==min(r),1,'first')
%caso todos os custos forem nao negativos, imprimimos a solucao otima
%juntamos xb e xn em um vetor x e o reordenamos à ordem original
%das variaveis iniciais já excluindo as variaveis artificiais, i.e.,
% x=[x1 x2 x3 ...xm]
if(r>=0)
  x=[xb;xn];
  indices=[indb indn];
  z=zeros(m+n,1);
  for i=1:m+n
    for j=1:m+n
      if((indices(i))==j)
        z(j)=x(i);
        break
      endif
    endfor
  endfor
  x=z(m+1:m+n,1);
  fprintf('estamos na solucao otima\nx*=\n')
  disp(x)
  fprintf('f=%f',f)
  return
endif
%direcao simplex By=an
y=B\A(:,indn(entra))
%se todos os y forem nao positivos retornamos que nao ha solucao otima finita
if(y<=0)
fprintf('nao existe solucao otima finita\n')
return
endif
%calculamos o passo pra ver quem sai da base
eps=M;
for i=1:m
  if(y(i)>0)
    eps=min(eps,xb(i)/y(i));
  endif
endfor
eps
%encontramos o indicie da variavel que sai da base e colocamos em sai
%pegamos o primeiro para evitar ter que escolher entre mais de uma 
%variavel candidata a sair
sai=find((xb./y)==eps,1,'first')
%por fim, atualizamos A, B, N e os indices
att=A(:,indn(entra));
A(:,indn(entra))=A(:,indb(sai));
A(:,indb(sai))=att
troca=cb(sai);
cb(sai)=cn(entra)
cn(entra)=troca
troca=indb(sai);
indb(sai)=indn(entra);
indn(entra)=troca;
B=A(:,1:m)
N=A(:,m+1:m+n)
%incremento na contagem
k=k+1
%aqui repetimos o processo acima para k>0
while(k>0)
xn=zeros(n,1);
xb=B\b
f=cb'*xb;
fprintf('o valor de f e %f\n',f)

%vetor multiplicador simplex B'lambda=cb
lambda=B'\cb

%custos relativos
r=zeros(n,1);
for i=1:n
  r(i)=cn(i)-lambda'*A(:,indn(i));
endfor
r
%variavel que entra na base
entra=find(r==min(r),1,'first')
%solucao otima
if(r>=0)
  x=[xb;xn];
  indices=[indb indn];
  z=zeros(m+n,1);
  for i=1:m+n
    for j=1:m+n
      if((indices(i))==j)
        z(j)=x(i);
        break
      endif
    endfor
  endfor
  x=z(m+1:m+n,1);
  fprintf('estamos na solucao otima\nx*=\n')
  disp(x)
  fprintf('f=%f\n',f)
  return
endif
%direcao simplex By=an
y=B\A(:,indn(entra))
%teste solucao ilimitada
if(y<=0)
fprintf('nao existe solucao otima finita\n')
return
endif
%passo pra ver quem sai da base
eps=M;
for i=1:m
  if(y(i)>0)
    eps=min(eps,xb(i)/y(i));
  endif
endfor
eps
%variavel que sai da base
sai=find((xb./y)==eps,1,'first')
%atualizacao
att=A(:,indn(entra));
A(:,indn(entra))=A(:,indb(sai));
A(:,indb(sai))=att
troca=cb(sai);
cb(sai)=cn(entra)
cn(entra)=troca
troca=indb(sai);
indb(sai)=indn(entra);
indn(entra)=troca;
B=A(:,1:m)
N=A(:,m+1:m+n)
%incremento
k=k+1
%acrescentamos uma condicao de parada para limite de iteracoes
%a principio estabelecemos 2000 iteracoes visto
%os problemas testados foram resolvidos em poucas iteracoes
%entao muitas iteracoes significa algum problema nao tratado nesse
%programa, como ciclagem da base
if(k>=2e3)
  fprintf('iteracoes demais :(')
  return
endif
endwhile
