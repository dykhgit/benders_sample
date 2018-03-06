#pragma once
#include <vector>
#include "type.h"
using namespace std;
class Maxflow
{
private:
	int s;
	int t;
	double eps;
	bool *side;
	int *buff;
	vector<vector<double> > data;
	void copyData(const vector<vector<double> > &g);
	void copyData(const IloNumMatrix &g);
	double maxflow(int s,int t,double eps,double up,bool *side);
	void dfs(int i,bool *side);

	void dfs(int i,bool *side,int *a,int val);
public:
	Maxflow(const vector<vector<double> > &g);
	Maxflow(){side=NULL;buff=NULL;};
	Maxflow(const IloNumMatrix &g);
	void resetGraph(const vector<vector<double> > &g);
	void resetGraph(const IloNumMatrix &g);
	double flow(int s,int t);
	void getSep(vector<vector<int> > &sep);
	void getSubGraph(vector<vector<int> > &sep);
	~Maxflow();
};

