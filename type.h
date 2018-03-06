#pragma once
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

	//  类型
typedef IloArray<IloNumArray>            IloNumMatrix;
typedef IloArray<IloNumMatrix>           IloNumMatrix3;
typedef IloArray<IloNumMatrix3>          IloNumMatrix4;

typedef IloArray<IloIntArray>            IloIntMatrix;
typedef IloArray<IloIntMatrix>           IloIntMatrix3;
typedef IloArray<IloNumVarArray>         IloNumVarMatrix;
typedef IloArray<IloNumVarMatrix>        IloNumVarMatrix3;
typedef IloArray<IloNumVarMatrix3>       IloNumVarMatrix4;
typedef IloArray<IloNumVarMatrix4>       IloNumVarMatrix5;
typedef IloArray<IloIntVarArray>         IloIntVarMatrix;
typedef IloArray<IloIntVarMatrix>        IloIntVarMatrix3;
typedef IloArray<IloIntVarMatrix3>       IloIntVarMatrix4;
typedef IloArray<IloIntVarMatrix4>       IloIntVarMatrix5;

typedef IloArray<IloRangeArray>          RangeMatrix;
typedef IloArray<RangeMatrix>            RangeMatrix3;
typedef IloArray<RangeMatrix3>           RangeMatrix4;

// 常量定义

