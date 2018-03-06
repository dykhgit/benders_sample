#include "Maxflow.h"


void Maxflow::copyData(const vector<vector<double> > &g)
{
	int i,j;
	if(data.size()<g.size()){
		data.resize(g.size());
		for(i=0;i<g.size();i++){
			data[i].resize(g[i].size());
		}
	}
	
	for(i=0;i<g.size();i++){
		for(j=0;j<g[i].size();j++){
			data[i][j]=g[i][j];
		}
	}
}

Maxflow::Maxflow(const vector<vector<double> > &g)
{
	data.resize(g.size());
	for(int i=0;i<g.size();i++){
		data[i].resize(g[i].size());
	}
	copyData(g);
	side=new bool[g.size()];
	buff=new int[g.size()];
}
void Maxflow::copyData(const IloNumMatrix &g){
	int i,j;
	if(data.size()<g.getSize()){
		data.resize(g.getSize());
		for(i=0;i<g.getSize();i++){
			data[i].resize(g[i].getSize());
		}
	}
	for(i=0;i<g.getSize();i++){
		for(j=0;j<g[i].getSize();j++){
			data[i][j]=g[i][j];
		}
	}
}

Maxflow::Maxflow(const IloNumMatrix &g)
{
	data.resize(g.getSize());
	for(int i=0;i<g.getSize();i++){
		data[i].resize(g[i].getSize());
	}
	copyData(g);
	side=new bool[g.getSize()];
	buff=new int[g.getSize()];
}

void Maxflow::resetGraph(const IloNumMatrix &g)
{
	if(side==NULL){
		side=new bool[g.getSize()];
		buff=new int[g.getSize()];
	}
	copyData(g);
}

void Maxflow::resetGraph(const vector<vector<double> > &g)
{
	if(side==NULL){
		side=new bool[g.size()];
		buff=new int[g.size()];
	}
	copyData(g);
}

double Maxflow::maxflow(int s,int t,double eps,double up,bool *side)  //计算最大流
{
	if(s==t) return up;
	double left=up,curflow=-1;
	side[s]=true;
	for(int i=0;i<data[s].size();i++)
	{
		if(data[s][i]>eps&&!side[i])
		{
			curflow=min(left,data[s][i]);
			curflow=maxflow(i,t,eps,curflow,side);
			left-=curflow;
			data[s][i]-=curflow;
			data[i][s]+=curflow;
		}
	}
	side[s]=false;
	return up-left;
}

double Maxflow::flow(int s,int t)
{
	this->s=s;
	this->t=t;
	this->eps=1e-3;
	fill(side,side+data.size(),false);
	return maxflow(s,t,eps,1e5,side);
}

void Maxflow::dfs(int i,bool *side)
{
	int j;
	side[i]=true;
	for(j=0;j<data.size();j++){
		if(data[i][j]>eps&&!side[j]){
			dfs(j,side);
		}
	}
}

void Maxflow::dfs(int i,bool *side,int *a,int val)
{
	int j;
	side[i]=true;
	a[i]=val;
	for(j=0;j<data.size();j++){
		if(data[i][j]>eps&&!side[j]){
			dfs(j,side,a,val);
		}
	}
}

void Maxflow::getSep(vector<vector<int> > &sep)
{
	fill(side,side+data.size(),false);
	dfs(s,side);
	if(sep.size()<2){
		sep.resize(2);
	}
	sep[0].clear();
	sep[1].clear();
	for(int i=0;i<data.size();i++){
		if(side[i]){
			sep[0].push_back(i);
		}else{
			sep[1].push_back(i);
		}
	}
}

void Maxflow::getSubGraph(vector<vector<int> > &sep)
{
	int cur=0,i;
	fill(side,side+data.size(),false);
	this->eps=1e-1;
	for(i=0;i<data.size();i++){
		if(!side[i]){
			dfs(i,side,buff,cur);
			cur+=1;
		}
	}
	this->eps=1e-3;
	sep.clear();
	sep.resize(cur);
	/*cout<<cur<<"--"<<sep.size()<<endl;
	 for(i=0;i<data.size();i++){
	 cout<<buff[i]<<" ";
	 }
	 cout<<endl;*/
	for(i=0;i<data.size();i++){
		sep[buff[i]].push_back(i);
	}
}

Maxflow::~Maxflow()
{
	delete [] side;
	delete [] buff;
}
