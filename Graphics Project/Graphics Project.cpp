#include "stdafx.h"
#include "Graphics Project.h"
#include <Commdlg.h>
#include <windows.h>
#include <string>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

//To-Do:
//Complete load: rest of algorithms
//Splines (make sure that it correct)

#define MAX_LOADSTRING 100

// Global Variables:
HINSTANCE hInst;                                // current instance
WCHAR szTitle[MAX_LOADSTRING];                  // The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];            // the main window class name

// Forward declarations of functions included in this code module:
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);

// Menu choice:
string choice = "";

// For changing colors:
COLORREF backgroundColor = RGB(255, 255, 255);
COLORREF drawColor = RGB(0, 0, 0);

// DRAWING ALGORITHMS
//Line
void DrawLine_MidPoint(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {

		int dx, dy, pk, x, y, d = 1;
		dx = x2 - x1;
		dy = y2 - y1;

		if (abs(dx) < abs(dy)) { // |slope| < 1
			if (y1 > y2) {
				swap(y2, y1);
				swap(x2, x1);
				dx *= -1;
				dy *= -1;
			}
			if (x2 < x1) {
				dx *= -1;
				d = -1;
			}
			pk = 2 * dx - dy;
			x = x1;
			y = y1;
			while (y < y2) {
				if (pk <= 0)
					pk += 2 * dx;
				else {
					pk += (2 * dx - 2 * dy);
					x += d;
				}
				SetPixel(hdc, x, y, drawColor);
				++y;
			}
			SetPixel(hdc, x2, y2, drawColor);
		}
		else {
			if (x1 > x2) {
				swap(x2, x1);
				swap(y2, y1);
				dx *= -1;
				dy *= -1;
			}
			if (y1 > y2) {
				dy *= -1;
				d = -1;
			}
			pk = 2 * dy - dx;
			x = x1;
			y = y1;
			SetPixel(hdc, x, y, drawColor);

			while (x < x2) {
				if (pk < 0) {
					pk += 2 * dy;
				}
				else {
					pk += (2 * dy - 2 * dx);
					y += d;
				}

				SetPixel(hdc, x, y, drawColor);
				++x;
			}
			SetPixel(hdc, x2, y2, drawColor);
		}
}

void DrawLine_Parametric(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color) {
	int x, y;
	double dt = 1.0 / 1000;
	for (double t = 0; t <= 1; t += dt) {
		x = (1 - t)*xs + t*xe;
		y = (1 - t)*ys + t*ye;
		SetPixel(hdc, x, y, color);
	}
}

void DrawLine_DDA(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color) {
	int dx = xe - xs, dy = ye - ys;
	SetPixel(hdc, xs, ys, color);
	double slope;
	//slope < 1
	if (abs(dx) >= abs(dy)) {
		slope = (double)dy / dx;

		int x = xs, xincr = dx > 0 ? 1 : -1; // we increase x 1 unit at a time
		double y = ys, yincr = slope*xincr;
		while (x != xe)
		{
			x += xincr; y += yincr;
			SetPixel(hdc, x, round(y), color);
		}
	}

	else {

		slope = (double)dx / dy;
		int y = ys, yincr = dy > 0 ? 1 : -1;
		double x = (double)xs, xincr = slope*yincr;
		while (y != ye)
		{
			x += xincr; y += yincr;
			SetPixel(hdc, round(x), y, color);
		}
	}
}

//Circle
void Draw8Points(HDC hdc, int xc, int yc, int a, int b, COLORREF color)
{
	SetPixel(hdc, xc + a, yc + b, color);
	SetPixel(hdc, xc - a, yc + b, color);
	SetPixel(hdc, xc - a, yc - b, color);
	SetPixel(hdc, xc + a, yc - b, color);
	SetPixel(hdc, xc + b, yc + a, color);
	SetPixel(hdc, xc - b, yc + a, color);
	SetPixel(hdc, xc - b, yc - a, color);
	SetPixel(hdc, xc + b, yc - a, color);
}

void DrawCircle_Polar(HDC hdc, int xc, int yc, int R, COLORREF color)
{
	int x = R, y = 0;
	double theta = 0, dtheta = 1.0 / R;
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x>y)
	{
		theta += dtheta;
		x = round(R*cos(theta));
		y = round(R*sin(theta));
		Draw8Points(hdc, xc, yc, x, y, color);
	}
}

void DrawCircle_IterativePolar(HDC hdc, int xc, int yc, int R, COLORREF color)
{
	double x = R, y = 0;
	double dtheta = 1.0 / R;
	double cdtheta = cos(dtheta), sdtheta = sin(dtheta);
	Draw8Points(hdc, xc, yc, R, 0, color);
	while (x>y)
	{
		double x1 = x*cdtheta - y*sdtheta;
		y = x*sdtheta + y*cdtheta;
		x = x1;
		Draw8Points(hdc, xc, yc, round(x), round(y), color);
	}
}

void DrawCircle_MidPoint(HDC hdc, int xc, int yc, int R, COLORREF color)
{
	int x = 0, y = R;
	int d = 1 - R;
	int c1 = 3, c2 = 5 - 2 * R;
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x<y)
	{
		if (d<0)
		{
			d += c1;
			c2 += 2;
		}
		else
		{
			d += c2;
			c2 += 4;
			y--;
		}
		c1 += 2;
		x++;
		Draw8Points(hdc, xc, yc, x, y, color);
	}
}

void DrawCircle_Cartesian(HDC hdc, int xc, int yc, int R, COLORREF color)
{
	int x = 0, y = R, R2 = R*R;
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x < y)
	{
		x++;
		y = round(sqrt((double)(R2 - x*x)));
		Draw8Points(hdc, xc, yc, x, y, color);
	}
}

//Curves
void DrawCurve_1stDegree(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
	int x, y;
	double dt = 1.0 / 1000;
	for (double t = 0; t <= 1; t += dt) {
		x = (1 - t)*x1 + t*x2;
		y = (1 - t)*y1 + t*y2;
		SetPixel(hdc, x, y, color);
	}
}

void DrawCurve_2ndDegree(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
	int x, y;
	double dt = 1.0 / 1000,
		slope = 500, // why? idk
		a1 = x2 - x1 - slope,
		a2 = y2 - y1 - slope,
		b1 = slope,
		b2 = slope,
		g1 = x1,
		g2 = y1;

	for (double t = 0; t <= 1; t += dt) {
		x = a1*t*t + b1*t + g1;
		y = a2*t*t + b2*t + g2;
		SetPixel(hdc, round(x), round(y), color);
	}
}

void mul(int A[4][4], double B[4], double C[4]) { //Multiply Hermite Matrix By Vector
	for (int i = 0; i < 4; i++) {
		C[i] = 0;
		for (int j = 0; j < 4; j++)
			C[i] += A[i][j] * B[j];
	}
}

double dot(double A[4], double B[4]) {
	double sum = 0;
	for (int i = 0; i < 4; i++)
		sum += A[i] * B[i];
	return sum;
}

void DrawCurve_Hermite(HDC hdc, POINT p1, POINT t1, POINT p2, POINT t2, COLORREF color) {
	double x, y;
	double dt = 1.0 / 10000;

	int basisMatrix[4][4] = { { 2, 1, -2, 1 },{ -3, -2, 3, -1 },{ 0, 1, 0, 0 },{ 1, 0, 0, 0 } };
	double vx[4] = { p1.x, t1.x, p2.x, t2.x };
	double vy[4] = { p1.y, t1.y, p2.y, t2.y };
	double gx[4], gy[4];

	mul(basisMatrix, vx, gx);
	mul(basisMatrix, vy, gy);

	for (double t = 0; t <= 1; t += dt) {
		double vt[4] = { t*t*t, t*t, t, 1 };
		x = dot(gx, vt);
		y = dot(gy, vt);
		SetPixel(hdc, round(x), round(y), color);
	}
}

void DrawCurve_Bezier(HDC hdc, POINT p0, POINT p1, POINT p2, POINT p3, COLORREF color) {
	// this algorithm is like hermite but it uses different slopes; the slope between p0, p1 and p2, p3
	POINT t1, t2;
	t1.x = 3 * (p1.x - p0.x);
	t1.y = 3 * (p1.y - p0.y);
	t2.x = 3 * (p3.x - p2.x);
	t2.y = 3 * (p3.y - p2.y);
	DrawCurve_Hermite(hdc, p0, t1, p3, t2, color);
}

void DrawCurve_Splines(HDC hdc, vector<POINT> P, int n, double c, COLORREF color) {
	double c1 = 1 - c;
	POINT T1, T2;
	T1.x = c1*(P[2].x - P[0].x);
	T1.y = c1*(P[2].y - P[0].y);
	
	for (int i = 2; i < n - 1; i++)
	{
		T2.x = c1*(P[i + 1].x - P[i - 1].x);
		T2.y = c1*(P[i + 1].y - P[i - 1].y);
		DrawCurve_Hermite(hdc, P[i - 1], T1, P[i], T2, color);
	}
}

//Filling
struct Entry
{
	int xmin, xmax;
};

void InitEntries(Entry table[])
{
	for (int i = 0; i<600; i++)
	{
		table[i].xmin = MAXINT;
		table[i].xmax = MININT;
	}
}

void ScanEdge(POINT v1, POINT v2, Entry table[])
{
	if (v1.y == v2.y)return;
	if (v1.y>v2.y)swap(v1, v2);
	double minv = (double)(v2.x - v1.x) / (v2.y - v1.y);
	double x = v1.x;
	int y = v1.y;
	while (y<v2.y)
	{
		if (x<table[y].xmin)table[y].xmin = (int)ceil(x);
		if (x>table[y].xmax)table[y].xmax = (int)floor(x);
		y++;
		x += minv;
	}
}

void DrawScanLines(HDC hdc, Entry table[], COLORREF color)
{
	for (int y = 0; y<600; y++)
		if (table[y].xmin<table[y].xmax)
			for (int x = table[y].xmin; x <= table[y].xmax; x++)
				SetPixel(hdc, x, y, color);
}

void ConvexFill(HDC hdc, vector <POINT> p, int n, COLORREF color)
{
	Entry *table = new Entry[600];
	InitEntries(table);
	POINT v1 = p[n - 1];
	for (int i = 0; i<n; i++)
	{
		POINT v2 = p[i];
		ScanEdge(v1, v2, table);
		v1 = p[i];
	}
	DrawScanLines(hdc, table, color);
	delete table;
}


//Clipping
void PointClipping_rect(HDC hdc, int x, int y, int xleft, int ytop, int xright, int ybottom, COLORREF color)
{
	if (x >= xleft && x <= xright && y >= ytop && y <= ybottom)
		SetPixel(hdc, x, y, color);
}

union OutCode
{
	unsigned All : 4;
	struct { unsigned left : 1, top : 1, right : 1, bottom : 1; };
};

OutCode GetOutCode(double x, double y, int xleft, int ytop, int xright, int ybottom)
{
	OutCode out;
	out.All = 0;
	if (x < xleft)
		out.left = 1;
	else if (x > xright)
		out.right = 1;
	if (y < ytop)
		out.top = 1;
	else if (y > ybottom)
		out.bottom = 1;
	return out;
}

void VIntersect(double xs, double ys, double xe, double ye, int x, double *xi, double *yi)
{
	*xi = x;
	*yi = ys + (x - xs)*(ye - ys) / (xe - xs);
}

void HIntersect(double xs, double ys, double xe, double ye, int y, double *xi, double *yi)
{
	*yi = y;
	*xi = xs + (y - ys)*(xe - xs) / (ye - ys);
}

void LineClipping_Rect(HDC hdc, int xs, int ys, int xe, int ye, int xleft, int ytop, int xright, int ybottom)
{
	double x1 = xs, y1 = ys, x2 = xe, y2 = ye;
	OutCode out1 = GetOutCode(x1, y1, xleft, ytop, xright, ybottom);
	OutCode out2 = GetOutCode(x2, y2, xleft, ytop, xright, ybottom);
	while ((out1.All || out2.All) && !(out1.All & out2.All))
	{
		double xi, yi;
		if (out1.All)
		{
			if (out1.left)
				VIntersect(x1, y1, x2, y2, xleft, &xi, &yi);
			else if (out1.top)
				HIntersect(x1, y1, x2, y2, ytop, &xi, &yi);
			else if (out1.right)
				VIntersect(x1, y1, x2, y2, xright, &xi, &yi);
			else 
				HIntersect(x1, y1, x2, y2, ybottom, &xi, &yi);
			x1 = xi;
			y1 = yi;
			out1 = GetOutCode(x1, y1, xleft, ytop, xright, ybottom);
		}
		else
		{
			if (out2.left)
				VIntersect(x1, y1, x2, y2, xleft, &xi, &yi);
			else if (out2.top)
				HIntersect(x1, y1, x2, y2, ytop, &xi, &yi);
			else if (out2.right)
				VIntersect(x1, y1, x2, y2, xright, &xi, &yi);
			else 
				HIntersect(x1, y1, x2, y2, ybottom, &xi, &yi);
			x2 = xi;
			y2 = yi;
			out2 = GetOutCode(x2, y2, xleft, ytop, xright, ybottom);
		}
	}

	if (!out1.All && !out2.All)
	{
		MoveToEx(hdc, round(x1), round(y1), NULL);
		LineTo(hdc, round(x2), round(y2));
	}
}

void PointClipping_Circle(HDC hdc, int x, int y, int xc, int yc, int r) {
	double delta = sqrt((x - xc)*(x - xc) + (y - yc)*(y - yc));
	if (delta < r) {
		SetPixel(hdc, x, y, drawColor);
	}
}

void LineClipping_Circle(HDC hdc, int x1, int y1, int x2, int y2, int xc, int yc, int r) {

	if (x1 > x2) {
		swap(x1, x2);  swap(y1, y2);
	}
	double xInter1, yInter1, xInter2, yInter2, delta, slope, d1, d2;

	slope = (y2 - y1)*1.0 / (x2 - x1)*1.0;

	xInter1 = -(sqrt((r*r - xc*xc + 2 * x1*xc - x1*x1)*slope*slope + ((2 * xc - 2 * x1)*yc - 2 * y1*xc + 2 * y1*x1)*slope + r*r - yc*yc + 2 * y1*yc - y1*y1) - x1*slope*slope + (y1 - yc)*slope - xc) / (double)(slope*slope + 1);
	yInter1 = y1 + (xInter1 - x1)*slope;

	xInter2 = (sqrt((r*r - xc*xc + 2 * x1*xc - x1*x1)*slope*slope + ((2 * xc - 2 * x1)*yc - 2 * y1*xc + 2 * y1*x1)*slope + r*r - yc*yc + 2 * y1*yc - y1*y1) + x1*slope*slope + (yc - y1)*slope + xc) / (double)(slope*slope + 1);
	yInter2 = y1 + (xInter2 - x1)*slope;

	d1 = sqrt((x1 - xc)*(x1 - xc) + (y1 - yc)*(y1 - yc));
	d2 = sqrt((x2 - xc)*(x2 - xc) + (y2 - yc)*(y2 - yc));

	if (d1 < r && d2 < r) {
		DrawLine_DDA(hdc, x1, y1, x2, y2, drawColor);
	}
	else if (d1 < r) {
		DrawLine_DDA(hdc, x1, y1, xInter2, yInter2, drawColor);
	}
	else if (d2 < r) {
		DrawLine_DDA(hdc, xInter1, yInter1, x2, y2, drawColor);
	}
	else {
		DrawLine_DDA(hdc, xInter1, yInter1, xInter2, yInter2, drawColor);
	}

}

//end of algorithms
/////////////////////////////////////////////////////////////////////////////////////////////

int calcR(POINT p1, POINT p2) { return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y)); }

string toString(string name, vector<POINT> v, COLORREF clr)
{
	string str = ""; 

	if (name == "convexfilling" || name == "curve_splines")
		str += name + "\n" + to_string(v.size());
	else
	str += name; 

	for (int i = 0; i < v.size(); i++)
		str += "\n" + to_string(v[i].x) + "\n" + to_string(v[i].y);

	int r = GetRValue(clr), g = GetGValue(clr), b = GetBValue(clr);
	str += "\n" + to_string(r) + "\n" + to_string(g) + "\n" + to_string(b) + "\n";
	
	return str;
}

COLORREF ShowColorDialog(HWND hwnd) {

	CHOOSECOLOR cc;
	static COLORREF crCustClr[16];

	ZeroMemory(&cc, sizeof(cc));
	cc.lStructSize = sizeof(cc);
	cc.hwndOwner = hwnd;
	cc.lpCustColors = (LPDWORD)crCustClr;
	cc.rgbResult = RGB(0, 255, 0);
	cc.Flags = CC_FULLOPEN | CC_RGBINIT;
	ChooseColor(&cc);

	return cc.rgbResult;
}

void save(string in) {

	OPENFILENAME ofn;

	char szFileName[MAX_PATH] = "";
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = NULL;
	ofn.lpstrFilter = (LPCWSTR)L"Text Files (*.txt)\0*.txt\0All Files (*.*)\0*.*\0";
	ofn.lpstrFile = (LPWSTR)szFileName;
	ofn.nMaxFile = MAX_PATH;
	ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
	ofn.lpstrDefExt = (LPCWSTR)L"txt";

	GetSaveFileName(&ofn);//open save as dialog

	wstring filePath = ofn.lpstrFile;

	ofstream outfile;

	outfile.open(filePath, ios::app);
	int r = GetRValue(backgroundColor), g = GetGValue(backgroundColor), b = GetBValue(backgroundColor);
	outfile << to_string(r) << endl << to_string(g) << endl << to_string(b) << endl; 
	outfile << in;

	outfile.close();
}

void load(HWND hwnd, HDC hdc) {

	OPENFILENAME ofn;

	TCHAR szFile[MAX_PATH];
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.lpstrFile = szFile;
	ofn.lpstrFile[0] = '\0';
	ofn.hwndOwner = hwnd;
	ofn.nMaxFile = sizeof(szFile);
	ofn.lpstrFilter = (LPCWSTR)L"Text Files (*.txt)\0*.txt\0All Files (*.*)\0*.*\0";
	ofn.nFilterIndex = 1;
	ofn.lpstrInitialDir = NULL;
	ofn.lpstrFileTitle = NULL;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

	GetOpenFileName(&ofn);

	wstring filePath = ofn.lpstrFile;
	string line;
	POINT p[10]; 
	int dia, r, g, b;
	ifstream infile;
	infile.open(filePath, ios::in);

	getline(infile, line);
	r = stoi(line);
	getline(infile, line);
	g = stoi(line);
	getline(infile, line);
	b = stoi(line);

	backgroundColor = RGB(r, g, b);

	RECT rc;
	GetClientRect(hwnd, &rc);
	HBRUSH brush = CreateSolidBrush(backgroundColor);
	FillRect(hdc, &rc, brush);

	while (getline(infile, line)) {
	
		if (line == "line_midpoint") {

			for (int i = 0; i < 2; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}
			
			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);

			DrawLine_MidPoint(hdc, p[0].x, p[0].y, p[1].x, p[1].y, RGB(r, g, b));
			continue;
		}

		else if (line == "line_dda") {
			
			for (int i = 0; i < 2; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);

			DrawLine_DDA(hdc, p[0].x, p[0].y, p[1].x, p[1].y, RGB(r, g, b));
			continue;
		}

		else if (line == "line_parametric") {

			for (int i = 0; i < 2; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);
			DrawLine_Parametric(hdc, p[0].x, p[0].y, p[1].x, p[1].y, RGB(r, g, b));
			continue;

		}

		else if (line == "circle_polar") {
			
			for (int i = 0; i < 2; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);
			dia = calcR(p[0], p[1]);
			DrawCircle_Polar(hdc, p[0].x, p[0].y, dia, RGB(r, g, b));
			continue;
		}

		else if (line == "circle_iterpolar") {
			
			for (int i = 0; i < 2; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);
			dia = calcR(p[0], p[1]);
			DrawCircle_IterativePolar(hdc, p[0].x, p[0].y, dia, RGB(r, g, b));
			continue;
		}

		else if (line == "circle_midpoint") {
			
			for (int i = 0; i < 2; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);
			dia = calcR(p[0], p[1]);
			DrawCircle_MidPoint(hdc, p[0].x, p[0].y, dia, RGB(r, g, b));
			continue;

		}

		else if (line == "circle_cartesian") {

			for (int i = 0; i < 2; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);
			dia = calcR(p[0], p[1]);
			DrawCircle_Cartesian(hdc, p[0].x, p[0].y, dia, RGB(r, g, b));
			continue;

		}

		else if (line == "curve_1st") {
			
			for (int i = 0; i < 2; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);

			DrawCurve_1stDegree(hdc, p[0].x, p[0].y, p[1].x, p[1].y, RGB(r, g, b));
			continue;
		}
		else if (line == "curve_2nd") {
			
			for (int i = 0; i < 2; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);

			DrawCurve_2ndDegree(hdc, p[0].x, p[0].y, p[1].x, p[1].y, RGB(r, g, b));
			continue;
		}
		else if (line == "curve_3rdh") {

			for (int i = 0; i < 2; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}
			
			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);
			p[2].x = 400; p[2].y = 400; p[3].x = 400; p[3].y = 400;
			DrawCurve_Hermite(hdc, p[0], p[2], p[1], p[3], RGB(r, g, b));
			continue;
		}
		else if (line == "curve_3rdb") {
			
			for (int i = 0; i < 4; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);

			DrawCurve_Bezier(hdc, p[0], p[1], p[2], p[3], RGB(r, g, b));
			continue;

		}
		else if (line == "curve_splines") {
			vector <POINT> pp;
			int sz;
			POINT p;

			getline(infile, line);
			sz = stoi(line);

			for (int i = 0; i < sz ; i++)
			{
				getline(infile, line);
				p.x = stoi(line);
				getline(infile, line);
				p.y = stoi(line);

				pp.push_back(p);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);

			DrawCurve_Splines(hdc, pp, pp.size(), 0.25, RGB(r, g, b));
			continue;

		}

		else if (line == "convexfilling") {

			vector <POINT> pp;
			POINT po;
			int sz;

			getline(infile, line);
			sz = stoi(line);

			for (int i = 0; i < sz; i++)
			{
				getline(infile, line);
				po.x = stoi(line);
				getline(infile, line);
				po.y = stoi(line);

				pp.push_back(po);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);

		ConvexFill(hdc, pp, pp.size(), RGB(r, g, b));

			continue;
		}

		else if (line == "lineClipping_rect") {

			for (int i = 0; i < 4; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);

			Rectangle(hdc, p[0].x, p[0].y, p[1].x, p[1].y);
			LineClipping_Rect(hdc, p[2].x, p[2].y, p[3].x, p[3].y, p[0].x, p[1].y, p[1].x, p[0].y);

			continue;
		}

		else if (line == "pointClipping_rect") {

			for (int i = 0; i < 3; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			Rectangle(hdc, p[0].x, p[0].y, p[1].x, p[1].y);

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);

			Rectangle(hdc, p[0].x, p[0].y, p[1].x, p[1].y);
			PointClipping_rect(hdc, p[2].x, p[2].y, p[0].x, p[1].y, p[1].x, p[0].y, RGB(r,g,b));
			continue;
		}

		else if (line == "lineClipping_circle") {

			for (int i = 0; i < 4; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);

			int r = calcR(p[0], p[1]);
			DrawCircle_MidPoint(hdc, p[0].x, p[0].y, r, RGB(r,g,b));
			LineClipping_Circle(hdc, p[2].x, p[2].y, p[3].x, p[3].y, p[0].x, p[0].y, r);
			continue;

		}
		else if (line == "pointClipping_circle") {

			for (int i = 0; i < 3; i++)
			{
				getline(infile, line);
				p[i].x = stoi(line);
				getline(infile, line);
				p[i].y = stoi(line);
			}

			getline(infile, line);
			r = stoi(line);
			getline(infile, line);
			g = stoi(line);
			getline(infile, line);
			b = stoi(line);

			int r = calcR(p[0], p[1]);
			DrawCircle_MidPoint(hdc, p[0].x, p[0].y, r, RGB(r, g, b));
			PointClipping_Circle(hdc, p[2].x, p[2].y, p[0].x, p[0].y, r);
			continue;

		}

	}
}

/////////////////////////////////////////////////////////////////////////////////
int APIENTRY wWinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE hPrevInstance, _In_ LPWSTR lpCmdLine, _In_ int nCmdShow)

{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: Place code here.

    // Initialize global strings
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_GRAPHICSPROJECT, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    // Perform application initialization:
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_GRAPHICSPROJECT));

    MSG msg;

    // Main message loop:
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    return (int) msg.wParam;
}

//
//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_GRAPHICSPROJECT));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW + 1);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_GRAPHICSPROJECT);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   hInst = hInstance; // Store instance handle in our global variable

   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
      CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);

   return TRUE;
}

//
//  FUNCTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  PURPOSE:  Processes messages for the main window.
//
//  WM_COMMAND  - process the application menu
//  WM_PAINT    - Paint the main window
//  WM_DESTROY  - post a quit message and return
//
//

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	
	HDC hdc;
	hdc = GetDC(hWnd);
	static string input = "";
	static int r, count = 0, numOfPoints = 0;
	static POINT p0, p1, p2, p3, p4, p5, t1, t2;
	static vector<POINT> p;
	
	switch (message)
	{

	case WM_RBUTTONDOWN: 

		if (choice == "convexfilling") {
			hdc = GetDC(hWnd);
			ConvexFill(hdc, p, numOfPoints, drawColor);
			input += toString("convexfilling", p, drawColor);
			numOfPoints = 0;
			p.clear();
			ReleaseDC(hWnd, hdc);

		}

		else if (choice == "curve_splines") {
			hdc = GetDC(hWnd);
			DrawCurve_Splines(hdc, p, p.size(), 0.5, drawColor);
			input += toString("curve_splines", p, drawColor);
			numOfPoints = 0;
			p.clear();
			ReleaseDC(hWnd, hdc);

		}
		
		break;

	case WM_LBUTTONDOWN:
	{
		
		if (choice == "line_midpoint") {
			if (count == 0) {
				p0.x = LOWORD(lParam);
				p0.y = HIWORD(lParam);
				p.push_back(p0);
				count++;
			}
			else if (count == 1) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				count = 0;
				p.push_back(p1);
				DrawLine_MidPoint(hdc, p0.x, p0.y, p1.x, p1.y, drawColor);
				input += toString("line_midpoint", p, drawColor);
			
				p.clear();
				ReleaseDC(hWnd, hdc);
			}
		}
		else if (choice == "line_dda") {
			if (count == 0) {
				p0.x = LOWORD(lParam);
				p0.y = HIWORD(lParam);
				p.push_back(p0);
				count++;
			}
			else if (count == 1) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				count = 0;
				p.push_back(p1);
				DrawLine_DDA(hdc, p0.x, p0.y, p1.x, p1.y, drawColor);
				input += toString("line_dda", p, drawColor);
				
				p.clear();
				ReleaseDC(hWnd, hdc);
			}
		}
		else if (choice == "line_parametric") {
			if (count == 0) {
				p0.x = LOWORD(lParam);
				p0.y = HIWORD(lParam);
				p.push_back(p0);
				count++;
			}
			else if (count == 1) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				count = 0;
				p.push_back(p1);
				DrawLine_Parametric(hdc, p0.x, p0.y, p1.x, p1.y, drawColor);
				input += toString("line_parametric", p, drawColor);
				p.clear();
				ReleaseDC(hWnd, hdc);
			}
		}
		else if (choice == "circle_polar") {
			if (count == 0) {
				p0.x = LOWORD(lParam);
				p0.y = HIWORD(lParam);
				p.push_back(p0);
				count++;
			}
			else if (count == 1) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				r = calcR(p0, p1);
				count = 0;
				p.push_back(p1);
				DrawCircle_Polar(hdc, p0.x, p0.y, r, drawColor);
				input += toString("circle_polar", p, drawColor);
				p.clear();
				ReleaseDC(hWnd, hdc);
			}
		}
		else if (choice == "circle_iterpolar") {
			if (count == 0) {
				p0.x = LOWORD(lParam);
				p0.y = HIWORD(lParam);
				p.push_back(p0);
				count++;
			}
			else if (count == 1) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				r = calcR(p0, p1);
				count = 0;
				p.push_back(p1);
				DrawCircle_IterativePolar(hdc, p0.x, p0.y, r, drawColor);
				input += toString("circle_iterpolar", p, drawColor);
				p.clear();
				ReleaseDC(hWnd, hdc);
			}
		}
		else if (choice == "circle_midpoint") {
			if (count == 0) {
				p0.x = LOWORD(lParam);
				p0.y = HIWORD(lParam);
				p.push_back(p0);
				count++;
			}
			else if (count == 1) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				r = calcR(p0, p1);
				count = 0;
				p.push_back(p1);
				DrawCircle_MidPoint(hdc, p0.x, p0.y, r, drawColor);
				input += toString("circle_midpoint", p, drawColor);
				p.clear();
				ReleaseDC(hWnd, hdc);
			}
		}
		else if (choice == "circle_cartesian") {
			if (count == 0) {
				p0.x = LOWORD(lParam);
				p0.y = HIWORD(lParam);
				p.push_back(p0);
				count++;
			}
			else if (count == 1) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				r = calcR(p0, p1);
				count = 0;
				p.push_back(p1);
				DrawCircle_Cartesian(hdc, p0.x, p0.y, r, drawColor);
				input += toString("circle_cartesian", p, drawColor);
				p.clear();
				ReleaseDC(hWnd, hdc);
			}
		}
		
		else if (choice == "curve_1st") {
			if (count == 0) {
				p0.x = LOWORD(lParam);
				p0.y = HIWORD(lParam);
				p.push_back(p0);
				count++;
			}
			else if (count == 1) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				count = 0;
				p.push_back(p1);
				DrawCurve_1stDegree(hdc, p0.x, p0.y, p1.x, p1.y, drawColor);
				input += toString("curve_1st", p, drawColor);
				p.clear();
				ReleaseDC(hWnd, hdc);
			}
		}
		else if (choice == "curve_2nd") {
			if (count == 0) {
				p0.x = LOWORD(lParam);
				p0.y = HIWORD(lParam);
				p.push_back(p0);
				count++;
			}
			else if (count == 1) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				count = 0;
				p.push_back(p1);
				DrawCurve_2ndDegree(hdc, p0.x, p0.y, p1.x, p1.y, drawColor);
				input += toString("curve_2nd", p, drawColor);
				p.clear();
				ReleaseDC(hWnd, hdc);
			}
		}
		else if (choice == "curve_3rdh") {
			if (count == 0) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				p.push_back(p1);
				count++;
			}
			else if (count == 1) {
				p2.x = LOWORD(lParam);
				p2.y = HIWORD(lParam);
				count = 0;
				t1.x = 400; t1.y = 400; t2.x = 400; t2.y = 400;
				DrawCurve_Hermite(hdc, p1, t1, p2, t2, drawColor);
				p.push_back(p2);
				input += toString("curve_3rdh", p, drawColor);
				p.clear();
				ReleaseDC(hWnd, hdc);
			}

		}
		else if (choice == "curve_3rdb") {
			if (count == 0) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				p.push_back(p1);
				count++;
			}
			else if (count == 1) {
				p2.x = LOWORD(lParam);
				p2.y = HIWORD(lParam);
				p.push_back(p2);
				count++;
			}
			else if (count == 2) {
				p3.x = LOWORD(lParam);
				p3.y = HIWORD(lParam);
				p.push_back(p3);
				count++;
			}
			else if (count == 3) {
				p4.x = LOWORD(lParam);
				p4.y = HIWORD(lParam);
				p.push_back(p4);
				count = 0;
				DrawCurve_Bezier(hdc, p1, p2, p3, p4, drawColor);
				input += toString("curve_3rdb", p, drawColor);
				p.clear();
				ReleaseDC(hWnd, hdc);
			}
		}
		else if (choice == "curve_splines") {

			numOfPoints++;
			p0.x = LOWORD(lParam);
			p0.y = HIWORD(lParam);
			p.push_back(p0);
		}
		else if (choice == "convexfilling") {

			numOfPoints++;
			p0.x = LOWORD(lParam);
			p0.y = HIWORD(lParam);
			p.push_back(p0);

		}
		else if (choice == "lineClipping_rect") {
			if (count == 0) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				p.push_back(p1);
				count++;
			}
			else if (count == 1) {
				hdc = GetDC(hWnd);
				p2.x = LOWORD(lParam);
				p2.y = HIWORD(lParam);
				p.push_back(p2);
				Rectangle(hdc, p1.x, p1.y, p2.x, p2.y);
				count++;
			}
			else if (count == 2) {
				p3.x = LOWORD(lParam);
				p3.y = HIWORD(lParam);
				p.push_back(p3);
				count++;
			}
			else if (count == 3) {
				p4.x = LOWORD(lParam);
				p4.y = HIWORD(lParam);
				p.push_back(p4);
				
				LineClipping_Rect(hdc, p3.x, p3.y, p4.x, p4.y, p1.x, p2.y, p2.x, p1.y);
				
				input += toString("lineClipping_rect", p, drawColor);
				p.clear();
				count = 0;
				ReleaseDC(hWnd, hdc);
			}
		}
		else if (choice == "pointClipping_rect") {
			if (count == 0) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				p.push_back(p1);
				count++;
			}
			else if (count == 1) {
				hdc = GetDC(hWnd);
				p2.x = LOWORD(lParam);
				p2.y = HIWORD(lParam);
				p.push_back(p2);
				Rectangle(hdc, p1.x, p1.y, p2.x, p2.y);
				count++;
			}
			else if (count == 2) {
				p3.x = LOWORD(lParam);
				p3.y = HIWORD(lParam);
				p.push_back(p3);
				PointClipping_rect(hdc, p3.x, p3.y, p1.x, p2.y, p2.x, p1.y, drawColor);
				input += toString("pointClipping_rect", p, drawColor);
				p.clear();
				count = 0;
				ReleaseDC(hWnd, hdc);
			}

		}
		else if (choice == "lineClipping_circle") {
			if (count == 0) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				p.push_back(p1);
				count++;
			}
			else if (count == 1) {
				p2.x = LOWORD(lParam);
				p2.y = HIWORD(lParam);
				p.push_back(p2);
				r = calcR(p1, p2);
				DrawCircle_MidPoint(hdc, p1.x, p1.y, r, drawColor);
				count++;
			}
			else if (count == 2) {
				p3.x = LOWORD(lParam);
				p3.y = HIWORD(lParam);
				p.push_back(p3);
				count++;
			}
			else if (count == 3) {
				p4.x = LOWORD(lParam);
				p4.y = HIWORD(lParam);
				p.push_back(p4);
				LineClipping_Circle(hdc, p3.x, p3.y, p4.x, p4.y, p1.x, p1.y, r);
				input += toString("lineClipping_circle", p, drawColor);
				p.clear();
				count = 0;
				ReleaseDC(hWnd, hdc);
			}
		}
		else if (choice == "pointClipping_circle") {
			if (count == 0) {
				p1.x = LOWORD(lParam);
				p1.y = HIWORD(lParam);
				p.push_back(p1);
				count++;
			}
			else if (count == 1) {
				p2.x = LOWORD(lParam);
				p2.y = HIWORD(lParam);
				p.push_back(p2);
				r = calcR(p1, p2);
				DrawCircle_MidPoint(hdc, p1.x, p1.y, r, drawColor);
				count++;
			}
			else if (count == 2) {
				p3.x = LOWORD(lParam);
				p3.y = HIWORD(lParam);
				p.push_back(p3);
				PointClipping_Circle(hdc, p3.x, p3.y, p1.x, p1.y, r);
				input += toString("pointClipping_circle", p, drawColor);
				p.clear();
				count = 0;
				ReleaseDC(hWnd, hdc);
			}
		}
	}

		break;
	
    case WM_COMMAND:
	{
		int wmId = LOWORD(wParam);
		// Parse the menu selections:
		switch (wmId)
		{
		case IDM_ABOUT:
			DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
			break;
		case IDM_EXIT:
			DestroyWindow(hWnd);
			break;
		case IDM_LINE_MIDPOINT:
			choice = "line_midpoint";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case IDM_LINE_DDA:
			choice = "line_dda";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case  IDM_LINE_PARAMETRIC:
			choice = "line_parametric";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case  IDM_CIRCLE_MIDPOINT:
			choice = "circle_midpoint";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case  IDM_CIRCLE_CARTESIAN:
			choice = "circle_cartesian";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case  IDM_CIRCLE_POLAR:
			choice = "circle_polar";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case  IDM_CIRCLE_ITERPOLAR:
			choice = "circle_iterpolar";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case  IDM_CURVE_1ST:
			choice = "curve_1st";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case  IDM_CURVE_2ND:
			choice = "curve_2nd";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case  IDM_CURVE_3RD_HERMITE:
			choice = "curve_3rdh";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case  IDM_CURVE_3RD_BEZIER:
			choice = "curve_3rdb";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case  IDM_CURVE_SPLINES:
			choice = "curve_splines";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case  IDM_CONVEX_FILLING:
			choice = "convexfilling";
			/* InvalidateRect(hWnd, 0, TRUE); */
			break;
		case IDM_POINTCLIP_RECT:
			choice = "pointClipping_rect";
			break;
		case IDM_LINECLIP_RECT:
			choice = "lineClipping_rect";
			break;
		case IDM_LINECLIP_CIRCLE:
			choice = "lineClipping_circle";
			break;
		case IDM_POINTCLIP_CIRCLE:
			choice = "pointClipping_circle";
			break;

		case IDM_CHANGE_BACKGROUND:
			backgroundColor = ShowColorDialog(hWnd);
			InvalidateRect(hWnd, NULL, TRUE);
			break;
		case IDM_CHANGE_DRAWCOLOR:
			drawColor = ShowColorDialog(hWnd);
			break;

		case IDM_SAVE:
			save(input);
			break;
		case IDM_LOAD:
			load(hWnd, hdc);
			break;

		default:
			return DefWindowProc(hWnd, message, wParam, lParam);
		}
	}
        break;
	case WM_ERASEBKGND:
	{
		RECT rc;
		GetClientRect(hWnd, &rc);
		HBRUSH brush = CreateSolidBrush(backgroundColor);
		FillRect(hdc, &rc, brush);
		return TRUE;

	}
    case WM_PAINT:
        {
            PAINTSTRUCT ps;
            HDC hdc = BeginPaint(hWnd, &ps);
            // TODO: Add any drawing code that uses hdc here...
            EndPaint(hWnd, &ps);
        }
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }

    return 0;
}

// Message handler for about box.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
