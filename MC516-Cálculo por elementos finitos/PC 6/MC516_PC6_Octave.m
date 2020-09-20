format long
K = [];
L = [];
NodesCondition = [];
ForcesCondition = [];
Elements = [];
E = 3.2e5;
P_A = 5000;
P_B = 4200;
P_C = 2500;
P_E = 3000;
A = (0.25*pi*(50)^2);
I = (pi*50^4)/64;
h = 1500;

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

function idk = ElementStiffness(l, angle, A, I, E)
    w = [A*cos(angle)^2+12*I*sin(angle)^2/(l*l),
         A*sin(angle)^2+12*I*cos(angle)^2/(l*l),
         (A-12*I/(l*l))*cos(angle)*sin(angle),
         6*I*sin(angle)/l,
         6*I*cos(angle)/l];
    aux = [[w(1),w(3),-w(4),-w(1),-w(3),-w(4)];
           [w(3),w(2),w(5),-w(3),-w(2),w(5)];
           [-w(4),w(5),4*I,w(4),-w(5),2*I];
           [-w(1),-w(3),w(4),w(1),w(3),w(4)];
           [-w(3),-w(2),-w(5),w(3),w(2),-w(5)];
           [-w(4),w(5),2*I,w(4),-w(5),4*I]];
    aux = aux*E/l;
    idk = aux;
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

Nodes = 5;
Nodes = 3*Nodes;
NumberOfElement = 6;
n = 10;
Nodes = Nodes+3*(n-2)*NumberOfElement;
NumberOfElement = NumberOfElement*(n-1);
L = [];
L2 = [];

PosNodes = [];
Elements = [];

for i=0:n-1
  PosNodes = [PosNodes; [i*h/(n-1),0]];
end

for i=1:n-1
  PosNodes = [PosNodes; [h,i*h/(n-1)]];
end

for i=1:n-1
  PosNodes = [PosNodes; [(h-i*h/(n-1)),h]];
end

for i=1:n-2
    PosNodes = [PosNodes; [i*h/(n-1),h-i*h/(n-1)]];
end

for i=1:n-1
    PosNodes = [PosNodes; [h+i*h/(n-1),h]];
end

for i=1:n-2
    PosNodes = [PosNodes; [2*h-i*h/(n-1),h-i*h/(n-1)]];
end

for i=0:4*n-6
  Elements = [Elements; [i+1,i+2]];
end

Elements = [Elements; [4*n-4,n]];
Elements = [Elements; [2*n-1,4*n-3]];

for i=4*n-3:6*n-8
   Elements = [Elements; [i,i+1]];
end
Elements = [Elements; [6*n-7,n]];

for i=1:NumberOfElement
  p1 = Elements(i,1);
  p2 = Elements(i,2);
  L = [L DistNodes(PosNodes(p1,1),PosNodes(p1,2),PosNodes(p2,1),PosNodes(p2,2))];
  L2 = [L2 DistNodes2(PosNodes(p1,1),PosNodes(p1,2),PosNodes(p2,1),PosNodes(p2,2))];
end

for i=1:NumberOfElement
  K = [K;ElementStiffness(L(i),L2(i),A,I,E)];
end

StiffnessMatrix = zeros(Nodes,Nodes);
U = zeros(Nodes,1);
F = zeros(Nodes,1);


StiffnessMatrix = Initialize(StiffnessMatrix,U,F,K,Elements);


[U,NodesCondition] = UBoundaryCondition(U,0,3*0+1,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,3*0+2,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,3*0+3,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,3*(3*n-3)+1,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,3*(3*n-3)+2,NodesCondition);
[U,NodesCondition] = UBoundaryCondition(U,0,3*(3*n-3)+3,NodesCondition);


[F,ForcesCondition] = FBoundaryCondition(F,-P_C,3*(n-1)+2,ForcesCondition);
[F,ForcesCondition] = FBoundaryCondition(F,-P_E,3*(2*n-2)+2,ForcesCondition);
[F,ForcesCondition] = FBoundaryCondition(F,-P_B,3*(5*n-6)+2,ForcesCondition);
[F,ForcesCondition] = FBoundaryCondition(F,P_A,3*(5*n-6)+1,ForcesCondition);

[U, F] = Solve(StiffnessMatrix, U, F, NodesCondition);

disp("Stiffness Matrix")
%disp(StiffnessMatrix)
disp("Displacements")
disp(U)
disp("Forces")
disp(F)

%colormap('jet')
%imagesc(StiffnessMatrix)

for i=1:NumberOfElement
  f = 0;
  inp = Elements(i,1);
  out = Elements(i,2);  
  x1 = -1*PosNodes(inp,1)-f*U(3*inp-2,1);
  y1 = PosNodes(inp,2)+f*U(3*inp-1,1);
  x2 = -1*PosNodes(out,1)-f*U(3*out-2,1);
  y2 = PosNodes(out,2)+f*U(3*out-1,1);
  plot(linspace(x1,x2,2),linspace(y1,y2,2),'--r','Linewidth',2)
  hold on
end

for i=1:NumberOfElement
  f = 4000;
  inp = Elements(i,1);
  out = Elements(i,2);  
  x1 = -1*PosNodes(inp,1)-f*U(3*inp-2,1);
  y1 = PosNodes(inp,2)+f*U(3*inp-1,1);
  x2 = -1*PosNodes(out,1)-f*U(3*out-2,1);
  y2 = PosNodes(out,2)+f*U(3*out-1,1);
  plot(linspace(x1,x2,2),linspace(y1,y2,2),'-bo','Linewidth',2.5)
  hold on
end
