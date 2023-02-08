#pragma once
#include "RJSplineToolKit.h"
#include "b_spline_interpolate.h"
RJVds MakSpiralCurve()
{
	double fR, fb, fstep,fxita;
	int N, i;
	RJVds mCurvVds;
	POINT3D theP;

	N = 25;
	fR = 10;
	fb = 4;
	fstep = PI/4.0;
	for (i = 0; i < N; i++)
	{
		fxita = i*fstep;
		theP.x = fR*cos(fxita);
		theP.y = fR*sin(fxita);
		theP.z = -fb*(N - 1) / 2 + fb*i;
		////
		mCurvVds.push_back(theP);
	}

	return mCurvVds;
}

int ConverDataArray(RJVds mCurvVds,double PointAarray[][3])
{
	int N, i;

	N = (int)mCurvVds.size();
	for (i = 0; i < N; i++)
	{
		PointAarray[i][0] = mCurvVds[i].x;
		PointAarray[i][1] = mCurvVds[i].y;
		PointAarray[i][2] = mCurvVds[i].z;
	}

	return N;
}

RJVds FitSplineCurve(RJVds mCurvVds)
{
	double PMat[50][3];
	int PNum;
	int PPNum[2];
	RJVds PsOfCurve;

	PNum = ConverDataArray(mCurvVds, PMat);
	PPNum[0] = PNum;
	PPNum[1] = 1000;
	PsOfCurve = RJOpenPara3DegSplCur(PMat, PPNum);
	return PsOfCurve;
}

RJVds FitSplineCurv2(RJVds mCurvVds)
{
	double PMat[50][3];
	int PNum;
	int PPNum[2];
	RJVds PsOfCurve;

	PNum = ConverDataArray(mCurvVds, PMat);
	PPNum[0] = PNum;
	PPNum[1] = 1000;
	PsOfCurve = RJClosedPara3DegSplCur(PMat, PPNum);
	return PsOfCurve;
}

RJVds MakNPolygonCurve()
{
	double fR, fstep, fxita;
	int N, i;
	RJVds mCurvVds;
	POINT3D theP;

	N = 7;
	fR = 10;
	fstep = 2*PI / 7.0;
	for (i = 0; i <= N; i++)
	{
		fxita = i*fstep;
		theP.x = fR*cos(fxita);
		theP.y = fR*sin(fxita);
		theP.z =0;
		////
		mCurvVds.push_back(theP);
	}

	return mCurvVds;
}

RJVds MakEllispeCurve()
{
	double fstep, fxita,fa,fb;
	int N, i;
	RJVds mCurvVds;
	POINT3D theP;

	N = 70;
	fa = 10;
	fb = 7;
	fstep = 2 * PI / 70.0;
	for (i = 0; i <= N; i++)
	{
		fxita = i*fstep;
		theP.x = fa*cos(fxita);
		theP.y = fb*sin(fxita);
		theP.z = 0;
		////
		mCurvVds.push_back(theP);
	}

	return mCurvVds;
}

RJVds OffsetPolygonCurve(RJDiscCurv m_TheCurveObj,double fOffSetValue)
{
	RJVdsSet OffSetCurs;
	RJVds PsOfCurve;

	OffSetCurs = m_TheCurveObj.RJOffsetSetCurv(fOffSetValue, true);
	PsOfCurve = OffSetCurs[0];

	return PsOfCurve;
}

void MakInterSplinSurf(RJVdsSet TheCurSet, double PPx[][151], double PPy[][151], double PPz[][151], int BPNum[][2])
{
	double PMat[310][3];
	int i, j, M, N, K;

	K = 10;//每条曲线上采集10个样点
	N = (int)TheCurSet.size();
	BPNum[0][0] = N; BPNum[0][1] = K;
	BPNum[1][0] = 3; BPNum[1][1] = 3;
	BPNum[2][0] = 100; BPNum[2][1] = 100;
	////
	for (i = 0; i < N; i++)
	{
		M = (int)TheCurSet[i].size();//第i个曲线的控制点个数
		for (j = 0; j < M; j++)
		{
			PMat[j][0] = TheCurSet[i][j].x;//控制点P[i][j]分量赋值
			PMat[j][1] = TheCurSet[i][j].y;
			PMat[j][2] = TheCurSet[i][j].z;
		}
		//
		K = ResampleTheCurv2(PMat, M, K);
			for (j = 0; j < K; j++)
		{
			PPx[i][j] = PMat[j][0];//表示第i个曲线上，采样的第j个数据点的x坐标(j=0,1,2...K-1)i=0,...N-1
			PPy[i][j] = PMat[j][1];
			PPz[i][j] = PMat[j][2];
		}
	}///for (i = 0; i < N; i++)
	////
	//MakInterSurf(PPx, PPy, PPz, BPNum);

	//(1)求解曲面信息
	double* Q = (double*)malloc(N * K * 3 * sizeof(double));
	for (i = 0; i < N; i++) {
		for (j = 0; j < K; j++) {
			Q[i * K * 3 + j * 3 + 0] = PPx[i][j];
			Q[i * K * 3 + j * 3 + 1] = PPy[i][j];
		}   Q[i * K * 3 + j * 3 + 2] = PPz[i][j];
	}

	int Q_dim[3] = { N,K,3 };

	//曲面信息的定义
	double* P = (double*)malloc(N * K * 3 * sizeof(double));
	int p = 3;
	int q = 3;
	double* U = (double*)malloc((N + p + 1) * sizeof(double));
	double* V = (double*)malloc((K + q + 1) * sizeof(double));
	interpolate_grid_points(Q, Q_dim, U, V, p, q, P);

	//(2)曲面两个方向计算离散点 100
	int P_dim[3] = { N,K,3 };
	SURFACE surf = { P,U,V,p,q,P_dim };

	double u0, v0;
	double S[3];
	for (i = 0; i <= 100; i++) {
		for (j = 0; j <= 100; j++) {
			u0 = 0.01 * i;
			v0 = 0.01 * j;
			calc_b_spline_surf(&surf, u0, v0, S);

			//复制到PP数组
			PPx[i][j] = S[0];
			PPy[i][j] = S[1];
			PPz[i][j] = S[2];
		}
	}

	free(Q);
	free(P);
	free(U);
	free(V);
}

void MakApproSplinSurf(RJVdsSet TheCurSet, double PPx[][151], double PPy[][151], double PPz[][151], int BPNum[][2])
{
	double PMat[310][3];
	int i, j, M, N, K;
	K = 10;//每条曲线上采集10个样点
	N = (int)TheCurSet.size();
	BPNum[0][0] = N; BPNum[0][1] = K;
	BPNum[1][0] = 3; BPNum[1][1] = 3;
	BPNum[2][0] = 100; BPNum[2][1] = 100;
	for (i = 0; i < N; i++)
	{
		M = (int)TheCurSet[i].size();
		for (j = 0; j < M; j++)
		{
			PMat[j][0] = TheCurSet[i][j].x;
			PMat[j][1] = TheCurSet[i][j].y;
			PMat[j][2] = TheCurSet[i][j].z;
		}
		K = ResampleTheCurv2(PMat, M, K);
		for (j = 0; j < K; j++)
		{
			PPx[i][j] = PMat[j][0];
			PPy[i][j] = PMat[j][1];
			PPz[i][j] = PMat[j][2];
		}
	}
	MakBSplinSurf(PPx, PPy, PPz, BPNum);
}