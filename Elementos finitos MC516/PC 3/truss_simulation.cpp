#pragma GCC optimize("Ofast")
#pragma GCC target("sse,sse2,sse3,ssse3")

#include <bits/stdc++.h>
#include </home/heneos/Sparse-Matrix/src/SparseMatrix/SparseMatrix.cpp>
using namespace std;
#define mp make_pair
#define pb push_back
#define FIFO ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0)
#define MOD 1000000007
#define PI 3.141592653589793
#define eps 1e-6

typedef long long ll;
typedef pair<int, int> pii;
typedef pair<double,double> pdd;
typedef pair<ll,ll> pll;
typedef vector<int> vi;
typedef vector<ll> vl;
typedef vector<pii> vii;
typedef vector<pll> vll;
typedef vector <vector <double> > vvd;

vector <vvd> K;
vector <double> A;
vector <pii> Elements;

vector <int> NodesCondition;
vector <int> ForcesCondition;
int NumberOfElement;
int Nodes;

vvd multiply(vvd &a, vvd &b){
	vvd x;
	int r = a.size();
	int c = b[0].size();
	for(int i=0; i<r; i++){
		vector <double> aux;
		for(int j=0; j<c; j++) aux.pb(0);
		x.pb(aux);
	}
	for(int i=0; i<r; i++){
		for(int j=0; j<c; j++){
			double ac = 0;
			for(int k=0; k<b.size(); k++) ac += a[i][k]*b[k][j];
			x[i][j] = ac;
		}
	}
	return x;
}

vvd conjugate_grad(vvd A, vvd b){
	int n = b.size();
	int what = A.size();
	int gg = A[0].size();
	vvd x;
	for(int i=0; i<n; i++){
		vector <double> aux;
		aux.pb(0);
		x.pb(aux);
	}
	vvd r = multiply(A, x);
	SparseMatrix <double> nA(what,gg);
	for(int i=0; i<what; i++){
		for(int j=0; j<gg; j++){
			nA.set(A[i][j],i+1,j+1);
		}
	}
	vector <double> p;
	double r_k_norm = 0;
	for(int i=0; i<n; i++){
		r[i][0] -= b[i][0];
		r_k_norm += (r[i][0]*r[i][0]);
		p.pb(-r[i][0]);
	}
	for(int i=0; i<2*n; i++){
		vector <double> Ap = nA*p;
		double aux = 0;
		for(int i=0; i<n; i++) aux += p[i]*Ap[i];
		double alpha = r_k_norm/aux;
		for(int i=0; i<n; i++){
			x[i][0] += alpha*p[i];
			r[i][0] += alpha*Ap[i];
		}
		double r_kplus_norm = 0; 
		for(int i=0; i<n; i++){
			r_kplus_norm += (r[i][0]*r[i][0]);
		}
		if(sqrt(r_kplus_norm) < 1e-6){
			break;
		}
		double beta = r_kplus_norm/r_k_norm;
		r_k_norm = r_kplus_norm;
		for(int i=0; i<n; i++){
			p[i] = beta*p[i]-r[i][0];
		}
	}
	return x;
}

pdd DistNodes(pdd f, pdd s){
	double aux;
	double x1 = f.first;
	double y1 = f.second;
	double x2 = s.first;
	double y2 = s.second;
	if(abs(x1-x2) < eps){
		aux = PI/2;
		if(y2 < y1) aux *= -1;
		return {sqrt((double)pow((x2-x1),2)+pow((y2-y1),2)),aux};
	}
	else{
		aux = atan((double)((y2-y1)/(x2-x1)));
		if(aux < 0 and y2 > y1) aux += PI;
		if(y2 < y1){
			aux += PI;
			if(x2 > x1) aux += PI;
		}
		return {sqrt((double)(pow((x2-x1),2)+pow((y2-y1),2))),aux};
	}
}

void UBoundaryCondition(vvd &nU, double u, int i){
	nU[i][0] = u;
	NodesCondition.pb(i);
}

void FBoundaryCondition(vvd &nF, double f, int i){
	nF[i][0] += f;
	ForcesCondition.pb(i);
}

void AssemblyStiffness(vvd &nS, vvd &k, int i, int j){
	for(int p=0; p<2; p++){
		for(int m=0; m<2; m++){
			nS[2*i+p][2*i+m] += k[p][m];
			nS[2*i+p][2*j+m] += k[p][2+m];
			nS[2*j+p][2*i+m] += k[p+2][m];
			nS[2*j+p][2*j+m] += k[p+2][m+2];
		}
	}
}

void Initialize(vvd  &nS, vvd &nU, vvd &nF){
	for(int i=0; i<Nodes; i++){
		nU[i][0] = 0;
		nF[i][0] = 0;
	}
	for(int i=0; i<NumberOfElement; i++){
		AssemblyStiffness(nS, K[i], Elements[i].first, Elements[i].second);	
	}
}

vvd PreSolvingS(vvd &S){
	int nsize = Nodes-NodesCondition.size();
	vvd newS;
	for(int i=0; i<nsize; i++){
		vector <double> aux;
		for(int j=0; j<nsize; j++) aux.pb(0);
		newS.pb(aux);
	}
	int contr = -1;
	for(int i=0; i<Nodes; i++){
		int contc = -1;
		bool flagr = false;
		for(int k=0; k < NodesCondition.size(); k++){
			if(i == NodesCondition[k]){
				flagr = true;
				break;
			}
		}
		if(flagr) continue;
		contr += 1;
		for(int j=0; j<Nodes; j++){
			bool flagc = false;
			for(int k=0; k<NodesCondition.size(); k++){
				if(j == NodesCondition[k]){
					flagc = true;
					break;
				}
			}
			if(flagc) continue;
			contc += 1;
			newS[contr][contc] = S[i][j];
		}
	}
	return newS;
}

vvd PreSolvingF(vvd nF, vvd &nS, vvd &nU){
	int nsize = Nodes-NodesCondition.size();
	vvd newF;
	for(int i=0; i<nsize; i++){
		vector <double> aux;
		aux.pb(0);
		newF.pb(aux);
	}
	int contr = -1;
	for(int i=0; i<Nodes; i++){
		bool flagr = false;
		for(int k=0; k<NodesCondition.size(); k++){
			if(i == NodesCondition[k]){
				flagr = true;
				break;
			}
		}
		if(flagr){
			for(int k=0; k<Nodes; k++) nF[k][0] = nF[k][0] - nS[k][i]*nU[i][0];
			continue;
		}
	}
	for(int i=0; i<Nodes; i++){
		bool flagr = false;
		for(int k=0; k<NodesCondition.size(); k++){
			if(i == NodesCondition[k]){
				flagr = true;
				break;
			}
		}
		if(flagr) continue;
		contr++;
		newF[contr][0] = nF[i][0];
	}
	return newF;
}

void Solve(vvd &nS, vvd &nU, vvd &nF){
	vvd newS = PreSolvingS(nS);
	vvd newF = PreSolvingF(nF,nS,nU);
	
	vvd u = conjugate_grad(newS, newF);
	
	int contr = -1;

	for(int i=0; i<Nodes; i++){
		bool flagr = false;
		for(int k=0; k<NodesCondition.size(); k++){
			if(i == NodesCondition[k]){
				flagr = true;
				break;
			}
		}
		if(flagr) continue;
		contr++;
		nU[i][0] = u[contr][0];
	}
	nF.clear();
	nF = multiply(nS,nU);
}

int main(){
	Nodes = 5;
	Nodes *= 2;
	NumberOfElement = 6;

	double E = 3.2e8;
	double h = 1500e-3;
	double d = 50e-3;
	double P_A = 5000e-3;
	double P_B = 4200e-3;
	double P_C = 2500e-3;
	double P_E = 3000e-3;

	for(int i=0; i<Nodes; i++) A.pb(0.25*PI*pow(d,2));
	
	vector <pdd> PosNodes = {mp(0,0),mp(h,0),mp(0,h),mp(h,h),mp(h,2*h)};
	Elements = {mp(0,2),mp(1,2),mp(1,3),mp(2,3),mp(2,4),mp(3,4)};
		

	pdd L[NumberOfElement];
	for(int i=0; i<NumberOfElement; i++)
		L[i] = DistNodes(PosNodes[Elements[i].first],PosNodes[Elements[i].second]);

	for(int i=0; i<NumberOfElement; i++){
		double aux[4][4];
		double angle = L[i].second;
		double rows[4] = {cos(angle),sin(angle),-cos(angle),-sin(angle)};
		for(int j=0; j<4; j++){
			for(int k=0; k<4; k++) aux[j][k] = rows[j]*rows[k];
		}
		vvd ans;
		for(int j=0; j<4; j++){
			vector <double> aux2;
			for(int k=0; k<4; k++) aux2.pb(E*aux[j][k]*A[i]/L[i].first);
			ans.push_back(aux2);
		}
		K.push_back(ans);
	}

	vvd StiffnessMatrix;

	for(int i=0; i<Nodes; i++){
		vector <double> aux;
		for(int j=0; j<Nodes; j++) aux.pb(0);
		StiffnessMatrix.pb(aux);
	}
		
	vvd U;
	vvd F;
	for(int i=0; i<Nodes;i++){
		vector <double> aux;
		aux.pb(0);
		U.pb(aux);
		F.pb(aux);
	}
	
	Initialize(StiffnessMatrix,U,F);
	
	//Node in UBoundary = Node*2+(x=0,y=1)
	UBoundaryCondition(U,0,2*0+0); //Nodo 0 en X
	UBoundaryCondition(U,0,2*0+1); //Nodo 0 en Y
	UBoundaryCondition(U,0,2*1+0); //Nodo 3 en X
	UBoundaryCondition(U,0,2*1+1); //Nodo 3 en Y

	FBoundaryCondition(F,-P_C,2*2+0); //Nodo 2 en X
	FBoundaryCondition(F,-P_E,2*3+0); //Nodo 3 en X
	FBoundaryCondition(F,-P_B,2*4+0); //Nodo 4 en X
	FBoundaryCondition(F,P_A,2*4+1); //Nodo 4 en Y
	
	
	Solve(StiffnessMatrix, U, F);
	
	for(int i=0; i<(int)U.size(); i++){
		cout << "[";
		for(int j=0; j<(int)U.size(); j++){
			printf("%+.5e",StiffnessMatrix[i][j]);
			if(j != (int)U.size()-1) cout << ", ";
		}
		cout << "]";
		cout << '\n';
	}	
	cout << "Displacements:\n";
	for(int i=0; i<(int)U.size(); i++) printf("%+.12e\n",U[i][0]);
	cout << "Forces:\n";
	for(int i=0; i<(int)F.size(); i++) printf("%+.12e\n",F[i][0]);
	
	return 0;
}
