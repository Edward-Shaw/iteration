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

//�������е����ֵ  
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

//�ſɱȵ�����ʽ�㷨ʵ��  
void Jacobi(vector<vector<double> > A, vector<double> B, int n)
{
	vector<double> X(n, 0);
	vector<double> Y(n, 0);
	vector<double> D(n, 0);
	int k = 0; //��¼ѭ������  
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
			cout << "����ʧ�ܣ��������Ǻ�����������" << endl;
			return;
		}

		for (int a = 0; a<n; a++){
			D[a] = X[a] - Y[a];
		}
	} while (MaxOfList(D)>0.000001 || MaxOfList(D)<-0.000001);

	return;
}

//ʹ���ſ˱ȵ����㷨�ⷽ����
void Call_Jacobi()
{
	cout << "������ִ�еĵ����㷨��: Jacobi" << endl;
	
	int n;
	cout << "�����뷽����δ֪���ĸ���n��";
	cin >> n;
	cout << endl;

	vector<vector<double> >A(n, vector<double>(n, 0));
	vector<double>B(n, 0);

	cout << "�����뷽�����ϵ������" << endl;
	for (int i = 0; i<n; i++){
		for (int j = 0; j<n; j++){
			cin >> A[i][j];
		}
	}
	cout << endl;

	cout << "�����뷽�����ֵ������" << endl;
	for (int k = 0; k<n; k++){
		cin >> B[k];
	}
	cout << endl;

	cout << "������ķ�����Ϊ��" << endl;
	for (int a = 0; a<n; a++){
		for (int b = 0; b<n; b++){
			cout << A[a][b] << " ";
		}
		cout << "    " << B[a] << endl;
	}
	cout << endl;
	cout << "���ſɱȵ�����ʽ��ķ�����Ľ�Ϊ��" << endl;
	Jacobi(A, B, n);

	return;
}

//�����ʼ����
void VectorInput(float x[], int n)
{
	int i;

	for (i = 1; i <= n; ++i)
	{
		printf("x[%d]=", i);
		scanf_s("%f", &x[i]);
	}
}

//�����������
void MatrixInput(float A[][MAX_n], int m, int n)
{
	int  i, j;
	printf("\n����ϵ������\n");
	for (i = 1; i <= m; ++i)
	{
		printf("�����������%d : ", i);
		for (j = 1; j <= n; ++j){
			scanf_s("%f", &A[i][j]);
		}
	}
}

//�������
void VectorOutput(float x[], int n)
{
	int i;
	for (i = 1; i <= n; ++i){
		printf("\nx[%d]=%f", i, x[i]);
	}
}

//�ж��Ƿ��ڹ涨������
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

//Gauss-Seidel�����㷨ʵ��
int GaussSeidel(float A[][MAX_n], float x[], int n)
{
	float x_former[MAX_n];
	int i, j, k;

	printf("\n��ʼ����x0:\n");
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
		printf("\nGauss-Seidel��������Ϊ%d ��", k);
		return 0;
	}
}

//ʹ��Gauss-Seidel�����㷨�ⷽ����
void Call_GaussSeide()
{
	cout << "������ִ�еĵ����㷨��: GaussSeide" << endl;

	int n;
	float A[MAX_n][MAX_n], x[MAX_n];

	printf("\n����ά��n=");
	scanf_s("%d", &n);
	if (n >= MAX_n - 1)
	{
		printf("\n\007n must < %d!", MAX_n);
		exit(0);
	}

	MatrixInput(A, n, n + 1);

	if (GaussSeidel(A, x, n)){
		printf("\nGauss-Seidel����ʧ��!");
	}
	else
	{
		printf("\n���:");
		VectorOutput(x, n);
	}
	
	return;
}

/*
 * SOR������������Ҫ�ı���
 */
float **A;    /*���A����*/
float *B;    /*���b����*/
float *X;    /*���x����*/
float p;    /*��ȷ��*/
float w;    /*�ɳ�����*/
int n;     /*δ֪������*/
int c;     /*����������*/
int k = 1;    /*ʵ�ʵ�������*/

//���ɳڵ���SOR�㷨����ʵ��(C���ʵ��)
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

//ʹ��SOR�����㷨�ⷽ���飨w = 1, 1.25, 1.5��
void Call_SOR(){
	cout << "������ִ�еĵ����㷨��: SOR" << endl;

	int i, j;
	printf("�������ά��N:\n");
	scanf_s("%d", &n);
	A = (float **)malloc(sizeof(float)*(n + 1));
	for (i = 0; i < n + 1; i++){
		A[i] = (float*)malloc(sizeof(float)*(n + 1));
	}

	printf("�������A:\n");
	for (i = 1; i < n + 1; i++){
		for (j = 1; j < n + 1; j++){
			scanf_s("%f", &A[i][j]);
		}
	}
	for (i = 1; i < n + 1; i++){
		for (j = 1; j < n; j++){
			if (A[i][j] == 0){
				printf("a[%d][%d]����Ϊ0\n", i, j);
			}
		}
	}
	B = (float *)malloc(sizeof(float)*(n + 1));
	printf("�������b:\n");
	for (i = 1; i < n + 1; i++){
		scanf_s("%f", &B[i]);
	}
	X = (float *)malloc(sizeof(float)*(n + 1));

	printf("�������x:\n");
	for (i = 1; i < n + 1; i++){
		scanf_s("%f", &X[i]);
	}
	printf("���뾫ȷֵ:\n");
	scanf_s("%f", &p);
	printf("��������������:\n");
	scanf_s("%d", &c);
	printf("�����ɳ����� w(0<w<2):\n");
	scanf_s("%f", &w);
	
	SOR(X);

	printf("\n���:");
	VectorOutput(X, n);

	return;
}

int _tmain(int argc, _TCHAR* argv[])
{
	string method;  //�����㷨
	string keepGoing;  //����ִ�б�ʶ
	do{
		cout << "����������㷨��Jacobi��Gauss-Seidel����SOR�������Իس�����������";
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
			cout << "��������ǣ�" << method << " �������޷��������ǰ������..." << endl;
		}
		cout << endl;
		cout << "������'quit'�˳���������������������ַ�����ִ�����������㷨" << endl;
		cin >> keepGoing;
		if (keepGoing.find("quit") != string::npos){
			break;
		}
	} while (true);
	
	printf("\n\n\007 ����ִ����ϣ���������˳�����!\n");

	system("PAUSE");

	return 0;
}
