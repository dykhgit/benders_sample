#pragma once

//控制 在分数点添加 subtour cut的频率

class Wait
{
public:
	Wait(int next){curTime=0;nextTime=next;}
	virtual bool isAdd() =0;
	virtual ~Wait(){}
protected:
	int curTime;
	int nextTime;
};

class ConstantWaiting:public Wait
{
public:
	virtual bool isAdd();
	ConstantWaiting(int c):Wait(c){}
	virtual ~ConstantWaiting(){}
};

class LinearWaiting:public Wait
{
public:
	virtual bool isAdd();
	LinearWaiting():Wait(1){}
	virtual ~LinearWaiting(){}
};
class  QuadraticWaiting:public Wait
{
public:
	virtual bool isAdd();
	QuadraticWaiting():Wait(1){}
	virtual ~QuadraticWaiting(){}
};

