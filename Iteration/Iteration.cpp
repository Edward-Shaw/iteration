//
// Iteration.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>  
#include <iomanip>  
#include <string>
#include <vector>
#include <stdio.h>
#include <math.h>

using namespace std;

#define PRECISION	0.000001  //误差精度e
#define MAX_NUMBER	100000  //最大迭代次数
#define N 12 //系数矩阵维数

//获得矩阵A
void Init_Matrix_A(double A[N][N]){

	cout << " 矩阵A：" << endl;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			A[i][j] = 1.0 / (i + j + 1);
			cout << A[i][i] << " ";
		}
	}
	cout << endl;
}

//获得值向量B
void Init_Vector_B(double B[N], double A[N][N]){
	cout << endl;
	cout << " 值向量B：" << endl;

	for (int i = 0; i <= (N - 1); i++)
	{
		for (int j = 0; j <= (N - 1); j++)
		{
			B[i] = B[i] + A[i][j];
		}
		cout << B[i] << " ";
	}

	cout << endl;
}

int _tmain(int argc, _TCHAR* argv[])
{
	//存放系数矩阵及值向量
	double A[N][N], B[N];

	//存放不同迭代算法计算后得到的值的变量
	double Jacobi[N], Jacobi_Temp[N], GaussSeidel[N], SOR_1[N], SOR_2[N], SOR_3[N];

	//存放迭代误差的变量
	double Jacobi_e[N], GaussSeidel_e[N], SOR1_e[N], SOR2_e[N], SOR3_e[N];
	
	int Iteration_Count = 0;  //迭代次数
	double a, a1;  //存放矩阵值的临时变量
	bool stop = false;
	
	cout << endl;
	cout << "当前系数矩阵维数为，N = " << N << endl;
	cout << endl;
	//计算并获取系数矩阵A
	Init_Matrix_A(A);
	
	//对变量进行初始化
	for (int i = 0; i < N; i++)
	{
		Jacobi[i] = 0;
		Jacobi_Temp[i] = 0;
		GaussSeidel[i] = 0;
		SOR_1[i] = 0;
		SOR_2[i] = 0;
		SOR_3[i] = 0;
		B[i] = 0;
	}

	//计算并获取值向量B
	Init_Vector_B(B, A);

	//Jacobi迭代
	while (true)
	{
		for (int u = 0; u < N; u++)
		{
			a = 0;
			for (int j = 0; j < N; j++)
			{
				if (u == j)
				{
					continue;
				}
				a = a + A[u][j] * Jacobi[j];
			}
			Jacobi_Temp[u] = (B[u] - a) / A[u][u];
		}

		Iteration_Count = Iteration_Count + 1;

		for (int i = 0; i < N; i++)
		{
			Jacobi_e[i] = Jacobi[i] - Jacobi_Temp[i];
			Jacobi_e[i] = fabs(Jacobi_e[i]);
		}

		for (int i = 0; i < N; i++)
		{
			if (Jacobi_e[i] >= PRECISION)
			{
				stop = false;
				break;
			}
			else{
				stop = true;
			}
		}

		for (int i = 0; i < N; i++){
			Jacobi[i] = Jacobi_Temp[i];
		}

		if (stop || (Iteration_Count >= MAX_NUMBER)){
			break;
		}
	}

	//显示输出
	cout << endl;
	cout << "****************************** Output ******************************" << endl;
	cout << endl;

	//Jacobi输出
	if (Iteration_Count >= MAX_NUMBER)
	{
		cout << " Jacobi 迭代结果：不收敛" << endl;
	}
	else
	{
		cout << " Jacobi 迭代结果：" << endl;
		for (int i = 0; i < N; i++)
		{
			cout << Jacobi[i] << endl;
		}
		cout << endl;
		cout << "迭代次数：" << Iteration_Count << endl;
	}
	cout << "----------------------------------" << endl;
	//重置迭代次数为0
	Iteration_Count = 0;

	//Gauss-Seidel迭代
	while (true)
	{
		for (int i = 0; i < N; i++)
		{
			a = 0;
			a1 = 0;
			for (int j = 0; j < N ; j++)
			{
				if (i == j){
					continue;
				}
				a = a + A[i][j] * GaussSeidel[j];
			}
			a1 = (B[i] - a) / A[i][i];
			GaussSeidel_e[i] = a1 - GaussSeidel[i];
			GaussSeidel_e[i] = fabs(GaussSeidel_e[i]);
			GaussSeidel[i] = a1;
		}

		Iteration_Count = Iteration_Count + 1;

		for (int u = 0; u < N; u++)
		{
			if (GaussSeidel_e[u] >= PRECISION)
			{
				stop = false;
				break;
			}
			else{
				stop = true;
			}
		}

		if (stop || (Iteration_Count >= MAX_NUMBER)){
			break;
		}
	}

	//Gauss-Seidel输出
	if (Iteration_Count >= MAX_NUMBER)
	{
		cout << " Gauss-Seidel 迭代结果：不收敛" << endl;
	}
	else
	{
		cout << endl;
		cout << " Gauss-Seidel 迭代结果：" << endl;
		for (int i = 0; i < N; i++)
		{
			cout << GaussSeidel[i] << endl;
		}
		cout << endl;
		cout << "迭代次数：" << Iteration_Count << endl;
	}
	cout << "----------------------------------" << endl;
	//重置迭代次数为0
	Iteration_Count = 0;

	//SQR迭代，w为1
	while (true)
	{
		for (int i = 0; i < N; i++)
		{
			a = 0;
			a1 = 0;
			for (int j = 0; j < N; j++)
			{
				if (i == j){
					continue;
				}
				a = a + A[i][j] * SOR_1[j];
			}
			a1 = (B[i] - a) / A[i][i];
			SOR1_e[i] = a1 - SOR_1[i];
			SOR1_e[i] = fabs(SOR1_e[i]);
			SOR_1[i] = a1;
		}

		Iteration_Count = Iteration_Count + 1;

		for (int u = 0; u < N; u++)
		{
			if (SOR1_e[u] >= PRECISION)
			{
				stop = false;
				break;
			}
			else{
				stop = true;
			}
		}

		if (stop || (Iteration_Count >= MAX_NUMBER)){
			break;
		}

	}

	//SQR输出,w为1
	if (Iteration_Count >= MAX_NUMBER)
	{
		cout << " SQR 迭代, w 为 1 结果：不收敛" << endl;
	}
	else
	{
		cout << endl;
		cout << " SQR 迭代, w 为 1 结果：" << endl;
		for (int i = 0; i < N; i++)
		{
			cout << SOR_1[i] << endl;;
		}
		cout << endl;
		cout << "迭代次数：" << Iteration_Count << endl;
	}
	cout << "----------------------------------" << endl;
	//重置迭代次数为0
	Iteration_Count = 0;

	//SQR迭代，w为1.25
	while (true)
	{
		for (int i = 0; i < N; i++)
		{
			a = 0;
			a1 = 0;
			for (int j = 0; j < N; j++)
			{
				if (i == j){
					continue;
				}
				a = a + A[i][j] * SOR_2[j];
			}
			a1 = 1.25*(B[i] - a) - 0.25*SOR_2[i] * A[i][i];
			SOR2_e[i] = SOR_2[i] - a1 / A[i][i];
			SOR2_e[i] = fabs(SOR2_e[i]);
			SOR_2[i] = a1 / A[i][i];
		}

		Iteration_Count = Iteration_Count + 1;

		for (int u = 0; u < N; u++)
		{
			if (SOR2_e[u] > (PRECISION))
			{
				stop = false;
				break;
			}
			else{
				stop = true;
			}
		}

		if (stop || (Iteration_Count >= MAX_NUMBER)){
			break;
		}

	}

	//SQR输出,w为1.25
	if (Iteration_Count >= MAX_NUMBER)
	{
		cout << " SQR 迭代, w 为 1.25 结果：不收敛" << endl;
	}
	else
	{
		cout << endl;
		cout << " SQR 迭代, w 为 1.25 结果：" << endl;
		for (int i = 0; i < N; i++)
		{
			cout << SOR_2[i] << endl;
		}
		cout << endl;
		cout << "迭代次数：" << Iteration_Count << endl;
	}
	cout << "----------------------------------" << endl;
	//重置迭代次数为0
	Iteration_Count = 0;

	//SQR迭代，w为1.5
	while (true)
	{
		for (int i = 0; i < N; i++)
		{
			a = 0;
			a1 = 0;
			for (int j = 0; j < N; j++)
			{
				if (i == j){
					continue;
				}
				a = a + A[i][j] * SOR_3[j];
			}
			a1 = (1 - 1.5)*SOR_3[i] + 1.000 / A[i][i] * 1.5*(B[i] - a);
			SOR3_e[i] = a1 - SOR_3[i];
			SOR3_e[i] = fabs(SOR3_e[i]);
			SOR_3[i] = a1;
		}

		Iteration_Count = Iteration_Count + 1;

		for (int u = 0; u < N; u++)
		{
			if (SOR3_e[u] >= PRECISION)
			{
				stop = false;
				break;
			}
			else{
				stop = true;
			}
		}

		if (stop || (Iteration_Count >= MAX_NUMBER)){
			break;
		}

	}

	//SQR输出,w为1.5
	if (Iteration_Count >= MAX_NUMBER)
	{
		cout << " SQR 迭代, w 为 1.5 结果：不收敛" << endl;
	}
	else
	{
		cout << endl;
		cout << " SQR 迭代, w 为 1.5 结果：" << endl;
		for (int i = 0; i < N; i++)
		{
			cout << SOR_3[i] << endl;
		}
		cout << endl;
		cout << "迭代次数：" << Iteration_Count << endl;
	}
	cout << "----------------------------------" << endl;
	cout << endl;
	cout <<"程序迭代执行完毕！";
	system("pause");
}