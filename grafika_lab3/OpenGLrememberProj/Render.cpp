

#include "Render.h"

#include <Windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>
#include <cmath>


double mX, mY, mZ, mT = 0;
double rX = 0, rY = 0, rZ = 0;

void animation() {
	double x, y, z = 0;
	double t;
	double A[] = { 1,3,-5 };
	double B[] = { 2,-7,4 };
	double C[] = { 9,5,1 };
	double D[] = { -4,3,6 };
	if (mT > 1) {
		mT = 0;
		rX = rY = rZ = 0;
	}

	t = mT;
	mX = pow((1 - t), 3) * A[0] + 3 * t*(pow((1 - t), 2))*B[0] + 3 * pow(t, 2)*(1 - t)*C[0] + pow(t, 3)*D[0];
	mY = pow((1 - t), 3)* A[1] + 3 * t*(pow((1 - t), 2))*B[1] + 3 * pow(t, 2)*(1 - t)*C[1] + pow(t, 3)*D[1];
	mZ = pow((1 - t), 3)* A[2] + 3 * t*(pow((1 - t), 2))*B[2] + 3 * pow(t, 2)*(1 - t)*C[2] + pow(t, 3)*D[2];

	t += 0.01;

	if (t <= 1) {
		x = pow((1 - t), 3) * A[0] + 3 * t*(pow((1 - t), 2))*B[0] + 3 * pow(t, 2)*(1 - t)*C[0] + pow(t, 3)*D[0];
		y = pow((1 - t), 3)* A[1] + 3 * t*(pow((1 - t), 2))*B[1] + 3 * pow(t, 2)*(1 - t)*C[1] + pow(t, 3)*D[1];
		z = pow((1 - t), 3)* A[2] + 3 * t*(pow((1 - t), 2))*B[2] + 3 * pow(t, 2)*(1 - t)*C[2] + pow(t, 3)*D[2];

		rX = atan2(y - mY, z - mZ) * 180 / 3.14;
		rY = atan2(x - mX, z - mZ) * 180 / 3.14;
		rZ = atan2(y - mY, x - mX) * 180 / 3.14;
	}

	mT += 0.01;
}

void Render()
{
	double x, y, z = 0;
	double A[] = {1,3,-5};
	double B[] = {2,-7,4};
	double C[] = {9,5,1};
	double D[] = {-4,3,6};
	animation();
	glColor3d(0.5, 0.4, 0);
	glBegin(GL_LINES);
	for (double t = 0; t <= 1; t += 0.001)
	{
		x = pow((1 - t), 3) * A[0] + 3 * t*(pow((1 - t), 2))*B[0] + 3 * pow(t, 2)*(1 - t)*C[0] + pow(t, 3)*D[0];
		y = pow((1 - t), 3)* A[1] + 3 * t*(pow((1 - t), 2))*B[1] + 3 * pow(t, 2)*(1 - t)*C[1] + pow(t, 3)*D[1];
		z = pow((1 - t), 3)* A[2] + 3 * t*(pow((1 - t), 2))*B[2] + 3 * pow(t, 2)*(1 - t)*C[2] + pow(t, 3)*D[2];
		glVertex3d(x, y, z);
	}
	glEnd();

	glTranslated(mX, mY, mZ);
	glRotatef(rX, 1, 0, 0);
	glRotatef(rY, 0, 1, 0);
	glRotatef(rZ, 0, 0, 1);

	glBegin(GL_TRIANGLES);
	glColor3d(0.0, 0.1, 0.2);
	glVertex3d(1, 0, 0);
	glVertex3d(0, 1, 0);
	glVertex3d(0, 0, 2);

	glColor3d(0.3, 0.4, 0.5);
	glVertex3d(-1, 0, 0);
	glVertex3d(0, 1, 0);
	glVertex3d(0, 0, 2);

	glColor3d(0.5, 0.4, 0.3);
	glVertex3d(1, 0, 0);
	glVertex3d(0, -1, 0);
	glVertex3d(0, 0, 2);

	glColor3d(0.2, 0.1, 0.0);
	glVertex3d(-1, 0, 0);
	glVertex3d(0, -1, 0);
	glVertex3d(0, 0, 2);

	glColor3d(0.5, 0.5, 0.2);
	glVertex3d(1, 0, 0);
	glVertex3d(0, 1, 0);
	glVertex3d(0, 0, -2);

	glColor3d(0.3, 0.4, 0);
	glVertex3d(-1, 0, 0);
	glVertex3d(0, 1, 0);
	glVertex3d(0, 0, -2);

	glColor3d(0, 0.4, 0.5);
	glVertex3d(1, 0, 0);
	glVertex3d(0, -1, 0);
	glVertex3d(0, 0, -2);

	glColor3d(0.7, 0.1, 0);
	glVertex3d(-1, 0, 0);
	glVertex3d(0, -1, 0);
	glVertex3d(0, 0, -2);

	glEnd();

	/*double A1[3] = { 2,1,3 };
	double B1[3] = { 4,4,7 };
	double AK[3] = { -1,-5,-3 };
	double BK[3] = { 1,7,-3 };
	glNormal3d(AK[0], AK[1], AK[2]);
	glNormal3d(BK[0], BK[1], BK[2]);
	glBegin(GL_LINES);
	glVertex3d(A1[0], A1[1], A1[2]);
	glVertex3d(A1[0] + AK[0], A1[1] + AK[1], A1[2] + AK[2]);
	glVertex3d(B1[0], B1[1], B1[2]);
	glVertex3d(B1[0] - BK[0], B1[1] - BK[1], B1[2] - BK[2]);
	double t;
	for (t = 0; t <= 1; t += 0.001) {
		x = A1[0] * (2 * pow(t, 3) - 3 * pow(t, 2) + 1) + B1[0] * (-2 * pow(t, 3) + 3 * pow(t, 2)) + AK[0] * (2 * pow(t, 3) - 3 * pow(t, 2) + t) + BK[0] * (pow(t, 3) - pow(t, 2));
		y = A1[1] * (2 * pow(t, 3) - 3 * pow(t, 2) + 1) + B1[1] * (-2 * pow(t, 3) + 3 * pow(t, 2)) + AK[1] * (2 * pow(t, 3) - 3 * pow(t, 2) + t) + BK[1] * (pow(t, 3) - pow(t, 2));
		z = A1[2] * (2 * pow(t, 3) - 3 * pow(t, 2) + 1) + B1[2] * (-2 * pow(t, 3) + 3 * pow(t, 2)) + AK[2] * (2 * pow(t, 3) - 3 * pow(t, 2) + t) + BK[2] * (pow(t, 3) - pow(t, 2));
		glVertex3d(x, y, z);
	}
	glEnd();

	double T[3] = { -3,1,-3 };
	double F[3] = { 4,-5,3 };
	double TK[3] = { 5,-3,3 };
	double FK[3] = { 1,5,89 };
	glNormal3d(TK[0], TK[1], TK[2]);
	glNormal3d(FK[0], FK[1], FK[2]);
	glBegin(GL_LINES);
	glVertex3d(T[0], T[1], T[2]);
	glVertex3d(T[0] + TK[0], T[1] + TK[1], T[2] + TK[2]);
	glVertex3d(F[0], F[1], F[2]);
	glVertex3d(F[0] - FK[0], F[1] - FK[1], F[2] - FK[2]);
	for (t = 0; t <= 1; t += 0.001) {
		x = T[0] * (2 * pow(t, 3) - 3 * pow(t, 2) + 1) + F[0] * (-2 * pow(t, 3) + 3 * pow(t, 2)) + TK[0] * (2 * pow(t, 3) - 3 * pow(t, 2) + t) + FK[0] * (pow(t, 3) - pow(t, 2));
		y = T[1] * (2 * pow(t, 3) - 3 * pow(t, 2) + 1) + F[1] * (-2 * pow(t, 3) + 3 * pow(t, 2)) + TK[1] * (2 * pow(t, 3) - 3 * pow(t, 2) + t) + FK[1] * (pow(t, 3) - pow(t, 2));
		z = T[2] * (2 * pow(t, 3) - 3 * pow(t, 2) + 1) + F[2] * (-2 * pow(t, 3) + 3 * pow(t, 2)) + TK[2] * (2 * pow(t, 3) - 3 * pow(t, 2) + t) + FK[2] * (pow(t, 3) - pow(t, 2));
		glVertex3d(x, y, z);
	}
	glEnd();


	double L[3] = { 2,2,3 };
	double M[3] = { 4,4,2 };
	double C1[3] = { 2,5,5 };
	double D1[3] = { 4,7,7 };

	glBegin(GL_LINES);
	glVertex3d(L[0], L[1], L[2]);
	glVertex3d(M[0], M[1], M[2]);
	glVertex3d(M[0], M[1], M[2]);
	glVertex3d(C1[0], C1[1], C1[2]);
	glVertex3d(C1[0], C1[1], C1[2]);
	glVertex3d(D1[0], D1[1], D1[2]);
	glVertex3d(D1[0], D1[1], D1[2]);
	glVertex3d(L[0], L[1], L[2]);
	for (t = 0; t <= 1; t += 0.01) {
		x = pow((1 - t), 3)*L[0] + 3 * t*pow((1 - t), 2)*M[0] + 3 * pow(t, 2)*(1 - t)*C1[0] + pow(t, 3)*D1[0];
		y = pow((1 - t), 3)*L[1] + 3 * t*pow((1 - t), 2)*M[1] + 3 * pow(t, 2)*(1 - t)*C1[1] + pow(t, 3)*D1[1];
		z = pow((1 - t), 3)*L[2] + 3 * t*pow((1 - t), 2)*M[2] + 3 * pow(t, 2)*(1 - t)*C1[2] + pow(t, 3)*D1[2];
		glVertex3d(x, y, z);
	}
	glEnd();

	double P[3] = { -2,4,3 };
	double O[3] = { 1,2,3 };
	double U[3] = { 5,5,5 };
	double Y[3] = { -5,7,-3 };

	glBegin(GL_LINES);
	glVertex3d(P[0], P[1], P[2]);
	glVertex3d(O[0], O[1], O[2]);
	glVertex3d(O[0], O[1], O[2]);
	glVertex3d(U[0], U[1], U[2]);
	glVertex3d(U[0], U[1], U[2]);
	glVertex3d(Y[0], Y[1], Y[2]);
	glVertex3d(Y[0], Y[1], Y[2]);
	glVertex3d(P[0], P[1], P[2]);
	for (t = 0; t <= 1; t += 0.01) {
		x = pow((1 - t), 3)*P[0] + 3 * t*pow((1 - t), 2)*O[0] + 3 * pow(t, 2)*(1 - t)*U[0] + pow(t, 3)*Y[0];
		y = pow((1 - t), 3)*P[1] + 3 * t*pow((1 - t), 2)*O[1] + 3 * pow(t, 2)*(1 - t)*U[1] + pow(t, 3)*Y[1];
		z = pow((1 - t), 3)*P[2] + 3 * t*pow((1 - t), 2)*O[2] + 3 * pow(t, 2)*(1 - t)*U[2] + pow(t, 3)*Y[2];
		glVertex3d(x, y, z);
	}
	glEnd();*/
}