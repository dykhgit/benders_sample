#include <iostream>
#include "Model.h"
#include <ilcplex/ilocplex.h>
using namespace std;

void test()
{
	IloEnv env;
	IloCplex cplex(env);
	IloModel model(env);

	IloNumVar x1(env),x2(env);
	IloExtractable ext1=model.add(IloRange(env,0,x1,2,"c1"));
	IloRange rng1;
	if(ext1.isConstraint()) {
		IloRangeI *impl = dynamic_cast<IloRangeI *>(ext1.asConstraint().getImpl());
		if (impl) {
			rng1=IloRange(impl);
		}
		else {
		}
	}

	IloExtractable ext2=model.add(IloRange(env,0,x2,2,"c2"));
	IloRange rng2;
	if(ext1.isConstraint()) {
		IloRangeI *impl = dynamic_cast<IloRangeI *>(ext2.asConstraint().getImpl());
		if (impl) {
			rng2=IloRange(impl);
		}
		else {
		}
	}
	IloExtractable ext3=model.add(IloRange(env,-IloInfinity,x1+2*x2,-1,"c3"));
	IloRange rng3;
	if(ext1.isConstraint()) {
		IloRangeI *impl = dynamic_cast<IloRangeI *>(ext3.asConstraint().getImpl());
		if (impl) {
			rng3=IloRange(impl);
		}
		else {
		}
	}
	model.add(IloMaximize(env,x1+x2));
	cplex.setParam(IloCplex::Param::Preprocessing::Presolve,IloFalse);
	cplex.extract(model);
	cplex.exportModel("tmp.lp");
	cplex.solve();
	//cplex.dualFarkas()
	IloNum num1=cplex.getDual(rng1);
	IloNum num2=cplex.getDual(rng2);
	IloNum num3=cplex.getDual(rng3);
	cout<<num1<<endl;
	cout<<num2<<endl;
	cout<<num3<<endl;
	model.end();
	cplex.end();
	env.end();
}
void test2()
{
	IloEnv env;
	IloCplex cplex(env);
	IloModel model(env);

	IloNumVar x1(env),x2(env),x3(env);
	model.add(IloMinimize(env,2*x1+2*x2-x3));
	model.add(x1+x3>=1);
	model.add(x2+2*x3>=1);
	IloNumVarArray varArray(env);
	IloNumArray numArray(env);
	cplex.setParam(IloCplex::Param::Preprocessing::Presolve,IloFalse);
	cplex.extract(model);
	cplex.solve();
	cout<<cplex.getStatus()<<endl;
	cplex.getRay(numArray,varArray);
	
	cout<<"size:"<<numArray.getSize()<<endl;
	for(int i=0;i<numArray.getSize();i++){
		if(varArray[i].getId()==x1.getId()){
			cout<<"x1=";
		}
		else if((varArray[i].getId()==x2.getId())){
			cout<<"x2=";
		}else{
			cout<<"x3=";
		}
		cout<<numArray[i]<<endl;
	}
	cout<<"print the value:"<<endl;
	cout<<cplex.getValue(x1)<<endl;
	cout<<cplex.getValue(x2)<<endl;
	cout<<cplex.getValue(x3)<<endl;
	model.end();
	cplex.end();
	env.end();
}

void test3()
{
	IloEnv env;
	IloCplex cplex(env);
	IloModel model(env);
	int *p;
	*p=100;
	IloNumVar x1(env),x2(env);
	IloRange rng1(env,0,x1,2,"c1");
	IloRange rng2(env,0,x2,2,"c2");
	x1.setObject(p);
	model.add(rng1);
	model.add(rng2);

	model.add(IloMaximize(env,x1+x2));
	cplex.setParam(IloCplex::Param::Preprocessing::Presolve,IloFalse);
	cplex.extract(model);
	
	cplex.solve();
	cout<<cplex.getObjValue()<<endl;

	rng1.setBounds(0,1);
	rng2.setBounds(0,1);
	cplex.solve();
	cout<<cplex.getObjValue()<<endl;

	int *p2=(int *)(x1.getObject());
	cout<<"the value is:"<<*p2<<endl;
	IloAny p3=x2.getObject();
	if(p3==NULL){
		cout<<"111"<<endl;
	}else{
		cout<<"222"<<endl;
	}
	model.end();
	cplex.end();
	env.end();
}
void test4()
{
	IloEnv env;
	IloCplex cplex(env);
	IloModel model(env);

	IloNumVarArray x(env,2,0,IloInfinity);
	IloRange c1(env,2,x[0]+x[1]);
	IloRange c2(env,x[1],1);
	model.add(c1);
	model.add(c2);

	model.add(IloMaximize(env,-2*x[0]-x[1]));
	cplex.setParam(IloCplex::Param::RootAlgorithm,IloCplex::Algorithm::Dual);
	cplex.extract(model);
	cplex.solve();
	double result=cplex.getObjValue();
	cout<<result<<endl;
	double u=cplex.getDual(c1);
	double v=cplex.getDual(c2);
	cout<<"u="<<u<<",v="<<v<<endl;
	model.end();
	cplex.end();
	env.end();
	int i;
	cin>>i;
}

//不检查长度是否合适，可能溢出
void getLogFileName(char *o,char *log){
	int i=-1;
	for(i=0;o[i]!='.'&&o[i]!='\0';i++){
		log[i]=o[i];
	}
	log[i]='.';
	log[i+1]='l';
	log[i+2]='o';
	log[i+3]='g';
	log[i+4]='\0';
}
int main(int argc,char*argv[])
{
	char *filename=NULL,*outfile=new char[100];
	int i;
	if(argc==1){
		filename=new char[100];
		cout<<"please input the filename:";
		cin>>filename;
		bool haveSuffix=false;
		for(i=0;filename[i]!='\0';i++){
			if(filename[i]=='.'){
				haveSuffix=true;
				break;
			}
		}
		if(!haveSuffix){
			filename[i]='.';
			filename[i+1]='t';
			filename[i+2]='x';
			filename[i+3]='t';
			filename[i+4]='\0';
		}
		getLogFileName(filename,outfile);
		Model m(filename,outfile);
		delete [] filename;
		cin>>i;
	}else{
		for(i=1;i<argc;i++){
			filename=argv[i];
			getLogFileName(filename,outfile);
			Model m(filename,outfile);
			
		}
	}
	cout<<"---------------------------------  end  ---------------------------"<<endl;
	delete [] outfile;
	return 0;
}