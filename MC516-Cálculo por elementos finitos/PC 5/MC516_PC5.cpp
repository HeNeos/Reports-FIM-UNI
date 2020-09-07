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
vector <pii> Elements;
vector <double> PosNodes;
vector <double> L;

vector <int> NodesCondition;
vector <int> ForcesCondition;
unordered_map <int,int> SearchCondition;
int NumberOfElement;
int Nodes;

vvd multiply(vvd &a, vvd &b){
	int r = a.size();
	int c = b[0].size();
	vvd x(r, vector <double>(c,0));
	for(int i=0; i<r; i++){
		for(int j=0; j<c; j++){
			double ac = 0;
			for(int k=0; k<(int)b.size(); k++) ac += a[i][k]*b[k][j];
			x[i][j] = ac;
		}
	}
	return x;
}

vvd conjugate_grad(vvd A, vvd b){
	int n = b.size();
	int rows = A.size();
	int columns = A[0].size();
	vvd x(n, vector <double> (1,0));	
	vvd r = multiply(A, x);
	SparseMatrix::SparseMatrix <double> nA(rows,columns);
	
	for(int i=0; i<rows; i++){
		for(int j=0; j<columns; j++)
			nA.set(A[i][j],i+1,j+1);
	}
	
	vector <double> p(n);
	double r_k_norm = 0;
	
	for(int i=0; i<n; i++){
		r[i][0] -= b[i][0];
		r_k_norm += (r[i][0]*r[i][0]);
		p[i] = -r[i][0];
	}

	for(int i=0; i<n*n; i++){
		vector <double> Ap = nA*p;
		double aux = 0;
		for(int i=0; i<n; i++) aux += p[i]*Ap[i];
		double alpha = r_k_norm/aux;
		for(int i=0; i<n; i++){
			x[i][0] += alpha*p[i];
			r[i][0] += alpha*Ap[i];
		}

		double r_kplus_norm = 0; 
		
		for(int i=0; i<n; i++)
			r_kplus_norm += (r[i][0]*r[i][0]);

		if(sqrt(r_kplus_norm) < 1e-6)
			break;
		
		double beta = r_kplus_norm/r_k_norm;
		r_k_norm = r_kplus_norm;
		
		for(int i=0; i<n; i++)
			p[i] = beta*p[i]-r[i][0];
	}
	return x;
}

void UBoundaryCondition(vvd &nU, double u, int i){
	nU[i][0] = u;
	NodesCondition.pb(i);
	SearchCondition[i] = 1;
}

void FBoundaryCondition(vvd &nF, double f, int i){
	nF[i][0] += f;
	ForcesCondition.pb(i);
}

void DFBoundaryCondition(vvd &nF, double w, int e){
	int f = Elements[e].first;
	int s = Elements[e].second;
	double l = L[e];
	FBoundaryCondition(nF,w*l/2,2*f+0);
	FBoundaryCondition(nF,w*l*l/12,2*f+1);
	FBoundaryCondition(nF,w*l/2,2*s+0);
	FBoundaryCondition(nF,-w*l*l/12,2*s+1);
}

vvd ElementStiffness(double E, double I, double l){
	vvd ES(4,vector <double> (4,0));
	double s1[4] = {12,6*l,-12,6*l};
	double s2[4] = {6*l,4*l*l,-6*l,2*l*l};
	for(int i=0; i<2; i++){
		for(int j=0; j<4; j++){
			ES[2*i][j] = s1[j];
			ES[2*i+1][j] = s2[j];
			if(i == 1) ES[2*i][j] *= -1;
		}
	}
	ES[3][1] = 2*l*l;
	ES[3][3] = 4*l*l;
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++) ES[i][j] *= (E*I/(l*l*l));
	}
	return ES;
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
	for(int i=0; i<NumberOfElement; i++)
		AssemblyStiffness(nS, K[i], Elements[i].first, Elements[i].second);	
}

vvd PreSolvingS(vvd &S){
	int nsize = Nodes-NodesCondition.size();
	vvd newS(nsize, vector <double> (nsize,0));	
	int contr = -1;
	for(int i=0; i<Nodes; i++){
		int contc = -1;
		if(SearchCondition.find(i) != SearchCondition.end()) continue;
		contr += 1;
		for(int j=0; j<Nodes; j++){
			if(SearchCondition.find(j) != SearchCondition.end()) continue;
			contc += 1;
			newS[contr][contc] = S[i][j];
		}
	}
	return newS;
}

vvd PreSolvingF(vvd nF, vvd &nS, vvd &nU){
	int nsize = Nodes-NodesCondition.size();
	vvd newF(nsize, vector <double>(nsize,1));
	int contr = -1;
	for(int i=0; i<Nodes; i++){
		if(SearchCondition.find(i) != SearchCondition.end()){
			for(int k=0; k<Nodes; k++) nF[k][0] = nF[k][0] - nS[k][i]*nU[i][0];
			continue;
		}
	}
	for(int i=0; i<Nodes; i++){
		if(SearchCondition.find(i) != SearchCondition.end()) continue;
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
		if(SearchCondition.find(i) != SearchCondition.end()) continue;
		contr++;
		nU[i][0] = u[contr][0];
	}
	nF.clear();
	nF = multiply(nS,nU);
}

int main(){
	Nodes = 11;
	Nodes *= 2;
	NumberOfElement = 10;

	double E = 3e5; //MPa
	double M_A = 2e5; //N-mm
	double I = 314.2222e4; //mm^4

	for(int i=0; i<Nodes/2; i++)
		PosNodes.pb(i*1000.0/(Nodes/2-1));	

	for(int i=0; i<NumberOfElement; i++)
		Elements.pb(mp(i,i+1));

	for(int i=0; i<NumberOfElement; i++)
		L.pb(PosNodes[Elements[i].second]-PosNodes[Elements[i].first]);


	for(int i=0; i<NumberOfElement; i++)
		K.pb(ElementStiffness(E,I,L[i]));

	vvd StiffnessMatrix(Nodes, vector <double> (Nodes,0));
	vvd U(Nodes, vector <double> (1,0));
	vvd F(Nodes, vector <double> (1,0));
	
	Initialize(StiffnessMatrix,U,F);
	
	UBoundaryCondition(U,0,2*0+0);
	UBoundaryCondition(U,0,2*NumberOfElement+0);
	FBoundaryCondition(F,M_A,2*NumberOfElement*7/10+1);

	Solve(StiffnessMatrix, U, F);

	cout << "Stiffness Matrix:\n";
	for(int i=0; i<(int)U.size(); i++){
		cout << "[";
		for(int j=0; j<(int)U.size(); j++){
			printf("%+.3e",StiffnessMatrix[i][j]);
			if(j != (int)U.size()-1) cout << ", ";
		}
		cout << "]";
		cout << '\n';
	}
	cout << "Displacements:\n";
	for(int i=0; i<(int)(U.size()/2); i++) printf("%+.12e %+.12e\n",U[2*i][0],U[2*i+1][0]);
	cout << "Forces:\n";
	for(int i=0; i<(int)(F.size()/2); i++) printf("%+.12e %+.12e\n",F[2*i][0],F[2*i+1][0]);
	return 0;
}
