//
// Iteration.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>  
#include <iomanip>  
#include <string>
#include <vector>
#include <stdio.h>
using namespace std;

#define	MAX_n	 100
#define PRECISION	0.000001
#define MAX_Number	1000

//求数组中的最大值  
double MaxOfList(vector<double>x)
{
	double max = x[0];
	int n = x.size();
	for (int i = 0; i < n; i++){
		if (x[i] > max){
			max = x[i];
		}
	}
	return max;
}

//雅可比迭代公式算法实现  
void Jacobi(vector<vector<double> > A, vector<double> B, int n)
{
	vector<double> X(n, 0);
	vector<double> Y(n, 0);
	vector<double> D(n, 0);
	int k = 0; //记录循环次数  
	do{
		X = Y;
		for (int i = 0; i<n; i++){
			double tem = 0;
			for (int j = 0; j<n; j++){
				if (i != j) tem += A[i][j] * X[j];
			}
			Y[i] = (B[i] - tem) / A[i][i];
			cout << left << setw(8) << Y[i] << " ";
		}
		cout << endl;
		k++;
		if (k>100){
			cout << "迭代失败！（可能是函数不收敛）" << endl;
			return;
		}

		for (int a = 0; a<n; a++){
			D[a] = X[a] - Y[a];
		}
	} while (MaxOfList(D)>0.000001 || MaxOfList(D)<-0.000001);

	return;
}

//使用雅克比迭代算法解方程组
void Call_Jacobi()
{
	cout << "您正在执行的迭代算法是: Jacobi" << endl;
	
	int n;
	cout << "请输入方程组未知数的个数n：";
	cin >> n;
	cout << endl;

	vector<vector<double> >A(n, vector<double>(n, 0));
	vector<double>B(n, 0);

	cout << "请输入方程组的系数矩阵：" << endl;
	for (int i = 0; i<n; i++){
		for (int j = 0; j<n; j++){
			cin >> A[i][j];
		}
	}
	cout << endl;

	cout << "请输入方程组的值向量：" << endl;
	for (int k = 0; k<n; k++){
		cin >> B[k];
	}
	cout << endl;

	cout << "您输入的方程组为：" << endl;
	for (int a = 0; a<n; a++){
		for (int b = 0; b<n; b++){
			cout << A[a][b] << " ";
		}
		cout << "    " << B[a] << endl;
	}
	cout << endl;
	cout << "由雅可比迭代公式求的方程组的解为：" << endl;
	Jacobi(A, B, n);

	return;
}

//输入初始向量
void VectorInput(float x[], int n)
{
	int i;

	for (i = 1; i <= n; ++i)
	{
		printf("x[%d]=", i);
		scanf_s("%f", &x[i]);
	}
}

//输入增广矩阵
void MatrixInput(float A[][MAX_n], int m, int n)
{
	int  i, j;
	printf("\n输入系数矩阵：\n");
	for (i = 1; i <= m; ++i)
	{
		printf("增广矩阵行数%d : ", i);
		for (j = 1; j <= n; ++j){
			scanf_s("%f", &A[i][j]);
		}
	}
}

//输出向量
void VectorOutput(float x[], int n)
{
	int i;
	for (i = 1; i <= n; ++i){
		printf("\nx[%d]=%f", i, x[i]);
	}
}

//判断是否在规定精度内
int IsSatisfyPricision(float x1[], float x2[], int n)
{
	int i;

	for (i = 1; i <= n; ++i){
		if (fabs(x1[i] - x2[i]) > PRECISION) {
			return 1;
		}
	}
	return 0;
}

//Gauss-Seidel迭代算法实现
int GaussSeidel(float A[][MAX_n], float x[], int n)
{
	float x_former[MAX_n];
	int i, j, k;

	printf("\n初始向量x0:\n");
	VectorInput(x, n);

	k = 0;
	do{
		for (i = 1; i <= n; ++i)
		{
			printf("\nx[%d]=%f", i, x[i]);
			x_former[i] = x[i];
		}
		printf("\n");
		for (i = 1; i <= n; ++i)
		{
			x[i] = A[i][n + 1];
			for (j = 1; j <= n; ++j){
				if (j != i){
					x[i] -= A[i][j] * x[j];
				}
				if (fabs(A[i][i]) > PRECISION){
					x[i] /= A[i][i];
				}
				else{
					return 1;
				}
			}
		}
		++k;
	} while (IsSatisfyPricision(x, x_former, n) && k<MAX_Number);

	if (k >= MAX_Number){
		return 1;
	}
	else
	{
		printf("\nGauss-Seidel迭代次数为%d 次", k);
		return 0;
	}
}

//使用Gauss-Seidel迭代算法解方程组
void Call_GaussSeide()
{
	cout << "您正在执行的迭代算法是: GaussSeide" << endl;

	int n;
	float A[MAX_n][MAX_n], x[MAX_n];

	printf("\n方阵维数n=");
	scanf_s("%d", &n);
	if (n >= MAX_n - 1)
	{
		printf("\n\007n must < %d!", MAX_n);
		exit(0);
	}

	MatrixInput(A, n, n + 1);

	if (GaussSeidel(A, x, n)){
		printf("\nGauss-Seidel迭代失败!");
	}
	else
	{
		printf("\n结果:");
		VectorOutput(x, n);
	}
	
	return;
}

/*
 * SOR迭代过程中需要的变量
 */
float **A;    /*存放A矩阵*/
float *B;    /*存放b矩阵*/
float *X;    /*存放x矩阵*/
float p;    /*精确度*/
float w;    /*松弛因子*/
int n;     /*未知数个数*/
int c;     /*最大迭代次数*/
int k = 1;    /*实际迭代次数*/

//超松弛迭代SOR算法具体实现(C风格实现)
void SOR(float xk[])
{
	int i, j;
	float t = 0.0;
	float tt = 0.0;
	float *xl;
	xl = (float *)malloc(sizeof(float)*(n + 1));
	for (i = 1; i<n + 1; i++)
	{
		t = 0.0;
		tt = 0.0;
		for (j = 1; j<i; j++)
			t = t + A[i][j] * xl[j];
		for (j = i; j<n + 1; j++)
			tt = tt + A[i][j] * xk[j];
		xl[i] = xk[i] + w*(B[i] - t - tt) / A[i][i];
	}

	t = 0.0;
	for (i = 1; i<n + 1; i++)
	{
		tt = fabs(xl[i] - xk[i]);
		tt = tt*tt;
		t += tt;
	}
	t = sqrt(t);

	for (i = 1; i<n + 1; i++)
		xk[i] = xl[i];

	if (k + 1 <= c&&t>p)
	{
		k++;
		SOR(xk);
	}
}

//使用SOR迭代算法解方程组（w = 1, 1.25, 1.5）
void Call_SOR(){
	cout << "您正在执行的迭代算法是: SOR" << endl;

	int i, j;
	printf("输入矩阵维数N:\n");
	scanf_s("%d", &n);
	A = (float **)malloc(sizeof(float)*(n + 1));
	for (i = 0; i < n + 1; i++){
		A[i] = (float*)malloc(sizeof(float)*(n + 1));
	}

	printf("输入矩阵A:\n");
	for (i = 1; i < n + 1; i++){
		for (j = 1; j < n + 1; j++){
			scanf_s("%f", &A[i][j]);
		}
	}
	for (i = 1; i < n + 1; i++){
		for (j = 1; j < n; j++){
			if (A[i][j] == 0){
				printf("a[%d][%d]不能为0\n", i, j);
			}
		}
	}
	B = (float *)malloc(sizeof(float)*(n + 1));
	printf("输入矩阵b:\n");
	for (i = 1; i < n + 1; i++){
		scanf_s("%f", &B[i]);
	}
	X = (float *)malloc(sizeof(float)*(n + 1));

	printf("输入矩阵x:\n");
	for (i = 1; i < n + 1; i++){
		scanf_s("%f", &X[i]);
	}
	printf("输入精确值:\n");
	scanf_s("%f", &p);
	printf("输入最大迭代次数:\n");
	scanf_s("%d", &c);
	printf("输入松弛因子 w(0<w<2):\n");
	scanf_s("%f", &w);
	
	SOR(X);

	printf("\n结果:");
	VectorOutput(X, n);

	return;
}

int _tmain(int argc, _TCHAR* argv[])
{
	string method;  //迭代算法
	string keepGoing;  //程序执行标识
	do{
		cout << "请输入迭代算法，Jacobi、Gauss-Seidel或者SOR（输入以回车键结束）：";
		cin >> method;
		cout << endl;
		if (method.find("Jacobi") != string::npos){
			Call_Jacobi();
		}
		else if (method.find("Gauss-Seide") != string::npos){
			Call_GaussSeide();
		}
		else if (method.find("SOR") != string::npos){
			Call_SOR();
		}
		else{
			cout << "你的输入是：" << method << " ，程序无法理解您当前的输入..." << endl;
		}
		cout << endl;
		cout << "请输入'quit'退出程序或者输入其他任意字符继续执行其他迭代算法" << endl;
		cin >> keepGoing;
		if (keepGoing.find("quit") != string::npos){
			break;
		}
	} while (true);
	
	printf("\n\n\007 迭代执行完毕，按任意键退出程序!\n");

	system("PAUSE");

	return 0;
}
