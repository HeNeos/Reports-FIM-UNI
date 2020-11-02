format long
A = [];
K = [];
NodesCondition = [];
ForcesCondition = [];
Elements = [];

function aux2 = DistNodes(f1,f2,s1,s2)
  f = [f1,f2];
  s = [s1,s2];
  if abs(s(1)- f(1)) < 1e-6
    aux = pi/2;
    if s(2) < f(2)
      aux = aux*-1;
    end
    aux2 = sqrt((s(1)-f(1))^2 + (s(2)-f(2))^2);
  else
    aux = atan((s(2)-f(2))/(s(1)-f(1)));
    if aux < 0 && s(2) > f(2)
      aux = aux + pi;
    end
    if s(2) < f(2)
      aux = aux + pi;
      if s(1) > f(1)
        aux = aux + pi;
      end
    end
    aux2 = sqrt((s(1)-f(1))^2 + (s(2)-f(2))^2);
  end
end

function aux = DistNodes2(f1,f2,s1,s2)
  f = [f1,f2];
  s = [s1,s2];
  if s(1) == f(1)
    aux = pi/2;
    if s(2) < f(2)
      aux = aux*-1;
    end
    aux2 = sqrt((s(1)-f(1))^2 + (s(2)-f(2))^2);
  else
    aux = atan((s(2)-f(2))/(s(1)-f(1)));
    if aux < 0 && s(2) > f(2)
      aux = aux + pi;
    end
    if s(2) < f(2)
      aux = aux + pi;
      if s(1) > f(1)
        aux = aux + pi;
      end
    end
    aux2 = sqrt((s(1)-f(1))^2 + (s(2)-f(2))^2);
  end
end

function [nU,NodesCondition] = UBoundaryCondition(nU,u,i,NodesCondition)
  NodesCondition = [NodesCondition i];
  nU(i,1) = u;
end

function [nF,ForcesCondition] = FBoundaryCondition(nF,f,i,ForcesCondition)
  nF(i,1) = nF(i,1)+f;
  ForcesCondition = [ForcesCondition i];
end

function nS = AssemblyStiffness(nS, z, k, i, j)
  for p = 1:2
    for m = 1:2
      nS(2*(i-1)+p,2*(i-1)+m) += k(4*(z-1)+p,m);
      nS(2*(i-1)+p,2*(j-1)+m) += k(4*(z-1)+p,2+m);
      nS(2*(j-1)+p,2*(i-1)+m) += k(4*(z-1)+p+2,m);
      nS(2*(j-1)+p,2*(j-1)+m) += k(4*(z-1)+p+2,2+m);
    end
  end
end

function nS = Initialize(nS,nU,nF, K, Elm)
  Nodes = size(nU)(1);
  NumberOfElement = size(K)(1)/size(K)(2);
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
  Nodes = size(S)(1);
  nsize = Nodes-size(NodesCondition)(2);
  newS = zeros(nsize,nsize);
  contr = 0;
  for i=1:Nodes
    contc = 0;
    flagr = 0;
    for k=1:size(NodesCondition)(2)
      if i==NodesCondition(k)
        flagr = 1;
        break
      end
    end
    if flagr == 1
      continue
    end
    contr += 1;
    for j=1:Nodes
      flagc = 0;
      for k=1:size(NodesCondition)(2)
        if j==NodesCondition(k)
          flagc = 1;
          break
        end
      end
      if flagc == 1
        continue
      end
      contc += 1;
      newS(contr,contc) = S(i,j);
    end
  end
end

function newF = PreSolvingF(nF,nS,nU,NodesCondition)
  Nodes = size(nF)(1);
  nsize = Nodes-size(NodesCondition)(2);
  newF = zeros(nsize,1);
  contr = 0;
  for i=1:Nodes
    flagr = 0;
    for k=1:size(NodesCondition)(2)
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
      for k=1:size(NodesCondition)(2)
        if i==NodesCondition(k)
          flagr = 1;
          break
        end
      end
      if flagr == 1
        continue
      end
      contr += 1;
      newF(contr,1) = nF(i,1);
    end
end

function [nU, nF] = Solve(nS, nU, nF, NodesCondition)
  Nodes = size(nS)(1);
  newS = PreSolvingS(nS,NodesCondition);
  newF = PreSolvingF(nF,nS,nU,NodesCondition);
  u = newS\newF;
  contr = 0;
  for i=1:Nodes
    flagr = 0;
    for k=1:size(NodesCondition)(2)
      if i == NodesCondition(k)
        flagr = 1;
        break
      end
    end
    if flagr == 1
      continue
    end
    contr += 1;
    nU(i,1) = u(contr,1);
  end
  nF = nS*nU;
end



Nodes = 5;
Nodes = 2*Nodes;
NumberOfElement = 6;

h = 1500e-3;
E = 3.2e8;
A = (0.25*pi*(50e-3)**2);
L = [];
L2 = [];
P_A = 5000e-3;
P_B = 4200e-3;
P_C = 2500e-3;
P_E = 3000e-3;

PosNodes = [[0,0];[h,0];[0,h];[h,h];[h,2*h]];
Elements = [[1,3];[2,3];[2,4];[3,4];[3,5];[4,5]];

for i=1:NumberOfElement
  p1 = Elements(i,1);
  p2 = Elements(i,2);
  L = [L DistNodes(PosNodes(p1,1),PosNodes(p1,2),PosNodes(p2,1),PosNodes(p2,2))];
  L2 = [L2 DistNodes2(PosNodes(p1,1),PosNodes(p1,2),PosNodes(p2,1),PosNodes(p2,2))];
end

for i=1:NumberOfElement
  aux = zeros(4);
  angle = L2(i);
  rows = [cos(angle),sin(angle),-cos(angle),-sin(angle)];
  for j=1:4
    for k=1:4
      aux(j,k) = rows(j)*rows(k);
    end
  end
  aux = aux*E*A/L(i);
  K = [K;aux];
end

StiffnessMatrix = zeros(Nodes,Nodes);
U = zeros(Nodes,1);
F = zeros(Nodes,1);

StiffnessMatrix = Initialize(StiffnessMatrix,U,F,K,Elements);

[U,NodesCondition] = UBoundaryCondition(U,0,2*0+1,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,2*0+2,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,2*1+1,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,2*1+2,NodesCondition);

[F,ForcesCondition] = FBoundaryCondition(F,-P_C,2*2+1,ForcesCondition);
[F,ForcesCondition] = FBoundaryCondition(F,-P_E,2*3+1,ForcesCondition);
[F,ForcesCondition] = FBoundaryCondition(F,-P_B,2*4+1,ForcesCondition);
[F,ForcesCondition] = FBoundaryCondition(F,P_A,2*4+2,ForcesCondition);

[U, F] = Solve(StiffnessMatrix, U, F, NodesCondition);

disp("Stiffness Matrix")
disp(StiffnessMatrix)
disp("Displacements")
disp(U)
disp("Forces")
