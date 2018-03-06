#include "Wait.h"

bool ConstantWaiting::isAdd()
{
	curTime+=1;
	if(curTime>=nextTime){
		curTime=0;
		return true;
	}
	return false;
}

bool LinearWaiting::isAdd()
{
	curTime+=1;
	if(curTime>=nextTime){
		curTime=0;
		nextTime+=1;
		return true;
	}
	return false;
}

bool QuadraticWaiting::isAdd()
{
	curTime+=1;
	if(curTime>=nextTime*nextTime){
		curTime=0;
		nextTime+=1;
		return true;
	}
	return false;
}
