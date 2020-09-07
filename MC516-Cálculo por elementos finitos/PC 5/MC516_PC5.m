format long
K = [];
NodesCondition = [];
ForcesCondition = [];
Elements = [];
E = 3e5;
M_A=2e5;
I = 315e4;

Nodes = 101;
NumberOfElement = 100;


PosNodes = [];
L = [];
Elements = zeros(NumberOfElement,2);

for i=1:Nodes
    PosNodes = [PosNodes,(i-1)*1000/(Nodes-1)];
end
Nodes = 2*Nodes;
for i=1:NumberOfElement
    Elements(i,1) = i;
    Elements(i,2) = i+1;
    L = [L, PosNodes(i+1) - PosNodes(i)];
end

for i=1:NumberOfElement
    aux = ElementStiffness(E,I,L(i));
    K = [K;aux];
end


StiffnessMatrix = zeros(Nodes,Nodes);
U = zeros(Nodes,1);
F = zeros(Nodes,1);

StiffnessMatrix = Initialize(StiffnessMatrix,U,F,K,Elements);

[U,NodesCondition] = UBoundaryCondition(U,0,2*0+1,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,2*NumberOfElement+1,NodesCondition);


[F,ForcesCondition] = FBoundaryCondition(F,M_A,2*700*NumberOfElement/1000+2,ForcesCondition);

[U, F] = Solve(StiffnessMatrix, U, F, NodesCondition);

disp("Stiffness Matrix")
disp(StiffnessMatrix)
disp("Displacements")
disp(U)
disp("Forces")
disp(F)


function [nU,NodesCondition] = UBoundaryCondition(nU,u,i,NodesCondition)
  NodesCondition = [NodesCondition i];
  nU(i,1) = u;
end

function [nF,ForcesCondition] = FBoundaryCondition(nF,f,i,ForcesCondition)
  nF(i,1) = nF(i,1)+f;
  ForcesCondition = [ForcesCondition i];
end

function [nF,ForcesCondition] = DFBoundaryCondition(nF,w,e,ForcesCondition)
  f = Elements(i,1);
  s = Elements(i,2);
  l = L(e);
  ForcesCondition = FBoundaryCondition(nF,w*l/2,2*f+0);
  ForcesCondition = FBoundaryCondition(nF,w*l*l/12,2*f+1);
  ForcesCondition = FBoundaryCondition(nF,w*l/2,2*s+0);
  ForcesCondition = FBoundaryCondition(nF,-w*l*l/12,2*s+1);
end

function idk = ElementStiffness(E, I, L)
    aux = zeros(4,4);
    s1 = [12,6*L,-12,6*L];
    s2 = [6*L,4*L*L,-6*L,2*L*L];
    for i=1:2
        for j=1:4
            aux(2*i-1,j) = s1(j);
            aux(2*i,j) = s2(j);
            if i==2
                aux(2*i-1,j) = -1*aux(2*i-1,j);
            end
        end
    end
    aux(4,2) = 2*L*L;
    aux(4,4) = 4*L*L;
    aux = aux*E*I/(L*L*L);
    idk = aux;
end

function nS = AssemblyStiffness(nS, z, k, i, j)
  for p = 1:2
    for m = 1:2
      nS(2*(i-1)+p,2*(i-1)+m) = nS(2*(i-1)+p,2*(i-1)+m)+k(4*(z-1)+p,m);
      nS(2*(i-1)+p,2*(j-1)+m) = nS(2*(i-1)+p,2*(j-1)+m)+k(4*(z-1)+p,2+m);
      nS(2*(j-1)+p,2*(i-1)+m) = nS(2*(j-1)+p,2*(i-1)+m)+k(4*(z-1)+p+2,m);
      nS(2*(j-1)+p,2*(j-1)+m) = nS(2*(j-1)+p,2*(j-1)+m)+k(4*(z-1)+p+2,2+m);
    end
  end
end

function nS = Initialize(nS,nU,nF, K, Elm)
  Nodes = size(nU,1);
  val1 = size(K,1);
  val2 = size(K,2);
  NumberOfElement = val1/val2;
  Elements = Elm;
  for i = 1:Nodes
    nU(i,1) = 0;
    nF(i,1) = 0;
  end
  for i = 1:NumberOfElement
    a = AssemblyStiffness(nS,i,K,Elements(i,1),Elements(i,2));
    nS = a;
  end
end

function newS = PreSolvingS(S,NodesCondition)
  Nodes = size(S,1);
  nsize = Nodes-size(NodesCondition,2);
  newS = zeros(nsize,nsize);
  contr = 0;
  for i=1:Nodes
    contc = 0;
    flagr = 0;
    for k=1:size(NodesCondition,2)
      if i==NodesCondition(k)
        flagr = 1;
        break
      end
    end
    if flagr == 1
      continue
    end
    contr = contr+1;
    for j=1:Nodes
      flagc = 0;
      for k=1:size(NodesCondition,2)
        if j==NodesCondition(k)
          flagc = 1;
          break
        end
      end
      if flagc == 1
        continue
      end
      contc = contc+1;
      newS(contr,contc) = S(i,j);
    end
  end
end

function newF = PreSolvingF(nF,nS,nU,NodesCondition)
  Nodes = size(nF,1);
  nsize = Nodes-size(NodesCondition,2);
  newF = zeros(nsize,1);
  contr = 0;
  for i=1:Nodes
    flagr = 0;
    for k=1:size(NodesCondition,2)
      if(i == NodesCondition(k))
        flagr = 1;
        break
      end
    end
    if flagr == 1
      for k=1:Nodes
        nF(k,1) = nF(k,1) - nS(k,1)*nU(i,1);
      end
      continue
    end
  end
    for i=1:Nodes
      flagr = 0;
      for k=1:size(NodesCondition,2)
        if i==NodesCondition(k)
          flagr = 1;
          break
        end
      end
      if flagr == 1
        continue
      end
      contr = contr+1;
      newF(contr,1) = nF(i,1);
    end
end

function [nU, nF] = Solve(nS, nU, nF, NodesCondition)
  Nodes = size(nS,1);
  newS = PreSolvingS(nS,NodesCondition);
  newF = PreSolvingF(nF,nS,nU,NodesCondition);
  u = newS\newF;
  contr = 0;
  for i=1:Nodes
    flagr = 0;
    for k=1:size(NodesCondition,2)
      if i == NodesCondition(k)
        flagr = 1;
        break
      end
    end
    if flagr == 1
      continue
    end
    contr = contr+1;
    nU(i,1) = u(contr,1);
  end
  nF = nS*nU;
end
