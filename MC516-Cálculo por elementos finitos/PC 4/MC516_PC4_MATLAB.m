format long
A = [];
K = [];
NodesCondition = [];
ForcesCondition = [];
Elements = [];

function aux = AllAngleC(f,s)
    L = DistanceNodes(f(1),f(2),f(3),s(1),s(2),s(3));
    C = [];
    for i=1:3
        C = [C (s(i)-f(i))/L];
    end
    C = acos(C);
    aux = C;
end

function aux = AllAngleL(f,s)
    L = DistanceNodes(f(1),f(2),f(3),s(1),s(2),s(3));
    C = [];
    for i=1:3
        C = [C (s(i)-f(i))/L];
    end
    C = acos(C);
    aux = L;
end

function aux = DistanceNodes(f1,f2,f3,s1,s2,s3)
    aux = sqrt((s1-f1)^2+(s2-f2)^2+(s3-f3)^2);
end

function aux = SingleAngle(f1,f2,s1,s2)
  f = [f1,f2];
  s = [s1,s2];
  if s(1) == f(1)
    aux = pi/2;
    if s(2) < f(2)
      aux = aux*-1;
    end
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
  for p = 1:3
    for m = 1:3
      nS(3*(i-1)+p,3*(i-1)+m) = nS(3*(i-1)+p,3*(i-1)+m)+k(6*(z-1)+p,m);
      nS(3*(i-1)+p,3*(j-1)+m) = nS(3*(i-1)+p,3*(j-1)+m)+k(6*(z-1)+p,3+m);
      nS(3*(j-1)+p,3*(i-1)+m) = nS(3*(j-1)+p,3*(i-1)+m)+k(6*(z-1)+p+3,m);
      nS(3*(j-1)+p,3*(j-1)+m) = nS(3*(j-1)+p,3*(j-1)+m)+k(6*(z-1)+p+3,3+m);
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


Nodes = 12;
Nodes = 3*Nodes;
NumberOfElement = 33;

E = 2.1e5;
K = [];
L = [];
C = [];
P_A = 10000;
P_B = 8000;
l1 = 600;
l2 = 500;
alpha = 30*pi/180;
beta = 70*pi/180;
A = pi*25*25;


PosNodes = [[0,0,0];[l1,0,0];[2*l1,0,0];[3*l1,0,0];[l1,-l1*tan(alpha),0];[2*l1,-l1*tan(alpha),0];[0,0,-l2];[l1,0,-l2];[2*l1,0,-l2];[3*l1,0,-l2];[l1,-l1*tan(alpha),-l2];[2*l1,-l1*tan(alpha),-l2]];

Elements = [[0,1];[0,4];[0,6];[0,10];
                     [1,2];[1,4];[1,6];[1,7];[1,8];[1,10];
                     [2,3];[2,4];[2,5];[2,8];[2,10];[2,11];
                     [3,5];[3,8];[3,9];[3,11];
                     [4,5];[4,10];
                     [5,10];[5,11];
                     [6,7];[6,10];
                     [7,8];[7,10];
                     [8,9];[8,10];[8,11];
                     [9,11];
                     [10,11]];

    
Elements = Elements+1;
for i=1:NumberOfElement
  p1 = PosNodes(Elements(i,1),1:3);
  p2 = PosNodes(Elements(i,2),1:3);
  L = [L;AllAngleL(p1,p2)];
  C = [C;AllAngleC(p1,p2)];
end

for i=1:NumberOfElement
  aux = zeros(6);
  angle = cos(C(i,1:3));
  rows = zeros(3);
  l = L(i);
  for j=1:3
      for k=1:3
          rows(j,k) = angle(j)*angle(k);
      end
  end
  for j=1:6
    for k=1:6
        s=1;
        if j > 3
            s = s*-1;
        end
        if k > 3
            s = s*-1;
        end
        aux(j,k) = rows(1+mod(j-1,3),1+mod(k-1,3))*s;
    end
  end
  aux = aux*E*A/l;
  K = [K;aux];
end

StiffnessMatrix = zeros(Nodes,Nodes);
U = zeros(Nodes,1);
F = zeros(Nodes,1);

StiffnessMatrix = Initialize(StiffnessMatrix,U,F,K,Elements);

[U,NodesCondition] = UBoundaryCondition(U,0,3*0+1,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,3*0+2,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,3*0+3,NodesCondition);

[U,NodesCondition] = UBoundaryCondition(U,0,3*6+1,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,3*6+2,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,3*6+3,NodesCondition);

[U,NodesCondition] = UBoundaryCondition(U,0,3*3+2,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,3*3+3,NodesCondition);

[U,NodesCondition] = UBoundaryCondition(U,0,3*9+2,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,3*9+3,NodesCondition);


[F,ForcesCondition] = FBoundaryCondition(F,-P_A/2,3*4+2,ForcesCondition);
[F,ForcesCondition] = FBoundaryCondition(F,-P_A/2,3*10+2,ForcesCondition);

[F,ForcesCondition] = FBoundaryCondition(F,P_B*sin(beta)/2,3*5+1,ForcesCondition);
[F,ForcesCondition] = FBoundaryCondition(F,P_B*sin(beta)/2,3*11+1,ForcesCondition);

[F,ForcesCondition] = FBoundaryCondition(F,-P_B*cos(beta)/2,3*5+2,ForcesCondition);
[F,ForcesCondition] = FBoundaryCondition(F,-P_B*cos(beta)/2,3*11+2,ForcesCondition);

[U, F] = Solve(StiffnessMatrix, U, F, NodesCondition);

disp("Stiffness Matrix")
disp(StiffnessMatrix)
disp("Displacements")
disp(U)
disp("Forces")
disp(F)


