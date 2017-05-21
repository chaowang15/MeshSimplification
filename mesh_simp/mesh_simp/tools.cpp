#include <string>
#include "tools.h"


// Compute the distance between two VERTICES
double getLineLength(Vertex *v0, Vertex *v1)
{
	return sqrt((v0->coord[0] - v1->coord[0])*(v0->coord[0] - v1->coord[0])
		+ (v0->coord[1] - v1->coord[1])*(v0->coord[1] - v1->coord[1])
		+ (v0->coord[2] - v1->coord[2])*(v0->coord[2] - v1->coord[2]));
}

// Compute the distance between two VERTICES
double getLineLength(double ax,double ay,double az,double bx,double by,double bz)
{
	return sqrt((ax-bx)*(ax-bx)+(ay-by)*(ay-by)+(az-bz)*(az-bz));
}


// Cross product
void crossProduct(Vertex *v0, Vertex *v1, Vertex *v2, double coord[3])
{
	double ax = v1->coord[0] - v0->coord[0], ay = v1->coord[1] - v0->coord[1], 
		az = v1->coord[2] - v0->coord[2];
	double bx = v2->coord[0] - v0->coord[0], by = v2->coord[1] - v0->coord[1], 
		bz = v2->coord[2] - v0->coord[2];
	coord[0] = ay * bz - az * by;
	coord[1] = az * bx - ax * bz;
	coord[2] = ax * by - ay * bx;
}

// Cross product
void crossProduct(Vector *vec1, Vector *vec2, double coord[3])
{
	double ax = vec1->coord[0], ay = vec1->coord[1], az = vec1->coord[2];
	double bx = vec2->coord[0], by = vec2->coord[1], bz = vec2->coord[2];
	coord[0] = ay * bz - az * by;
	coord[1] = az * bx - ax * bz;
	coord[2] = ax * by - ay * bx;	
}

// Dot Product
double dot3(Vector *vec1, Vector *vec2)
{
	return (vec1->coord[0]*vec2->coord[0] + vec1->coord[1]*vec2->coord[1]
	+ vec1->coord[2]*vec2->coord[2]);
}

// Compute the length of a vector
double vectorLength(double v[3])
{
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// Add a triangle's normal onto a vertex's normal
void addTriNorm2VtxNorm(Vertex *v, Triangle *t)
{
	v->norm[0] += t->norm[0];
	v->norm[1] += t->norm[1];
	v->norm[2] += t->norm[2];
}

bool IsEdgeInTriangle(int v0, int v1, Triangle *t)
{
	if (v0 == v1)
		return false;
	if (isIncluded(v0, t->vtx[0], t->vtx[1], t->vtx[2]) != -1 && 
		isIncluded(v1, t->vtx[0], t->vtx[1], t->vtx[2]) != -1)
		return true;
	else
		return false;
}


// Compute the angle (in degree) between two input vectors
double computeAngleOfVectors(Vector *vec1, Vector *vec2)
{
	double len1 = sqrt(vec1->coord[0]*vec1->coord[0] + 
		vec1->coord[1]*vec1->coord[1] + vec1->coord[2]*vec1->coord[2]);
	double len2 = sqrt(vec2->coord[0]*vec2->coord[0] + 
		vec2->coord[1]*vec2->coord[1] + vec2->coord[2]*vec2->coord[2]);
	double vecDot = dot3(vec1, vec2);
	return 180 * acos(vecDot/(len1*len2)) / (double)PI;
}

// Compute the angle (in degree) between two input vectors
double computeAngleOfVectors(double vec1[3], double vec2[3])
{
	double len1 = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
	double len2 = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1] + vec2[2]*vec2[2]);
	double vecDot = dot3(vec1, vec2);
	return 180 * acos(vecDot/(len1*len2)) / (double)PI;
}

// Compute the distance between an input point and an input
// line determined by two input VERTICES
double computeDisBtwnPointAndLine(double p[3], double v0[3], double v1[3])
{
	double m[3], n[3];
	for (int i = 0; i < 3; ++i)
	{
		m[i] = v1[i] - v0[i];
		n[i] = p[i] - v0[i];
	}
	double a = dot3(m, n);
	double mLen = vectorLength(m);
	double nLen = vectorLength(n);
	return sqrt(nLen*nLen-a*a/(mLen*mLen));
}

// Given the full path of a file, find its file name 
void getFileName(char *filePath, char *fileName)
{
	int i = 0, j = 0;
	int len = strlen(filePath);
	for (i = len-1; i >= 0; i--)
	{
		if (filePath[i] == '\\' || filePath[i] == '/')
		{
			break;
		}
	}
	for (j = 0; j<len-i-1; ++j)
	{
		fileName[j] = filePath[i+j+1];
	}
	fileName[j] = '\0';
}

void getFileSuffix(char *filePath, char *fileSuffix)
{
	int i,j;
	int len = strlen(filePath);
	for (i = len-1; i >= 0; i--)
	{
		if (filePath[i] == '.')
		{
			break;
		}
	}
	for (j = 0; j<len-i; ++j)
	{
		fileSuffix[j] = filePath[i+j];
	}
	fileSuffix[j] = '\0';
}

int get3rdVtrOfTriangle(Triangle *t, int a, int b)
{
	int m = isIncluded(a, t->vtx[0], t->vtx[1], t->vtx[2]);
	int n = isIncluded(b, t->vtx[0], t->vtx[1], t->vtx[2]);
	if (m == -1 || n == -1)
		return -1;
	else
		return 3-(m+n);
}

void sort3(unsigned int &v0, unsigned int &v1, unsigned int &v2)
{
	unsigned int a = (v0 < v1) ? v0 : v1;
	a = (a < v2) ? a : v2;
	unsigned int b = (v0 > v1) ? v0 : v1;
	b = (b > v2) ? b : v2;
	unsigned int c = v0 + v1 + v2;
	v0 = a;
	v1 = c-a-b;
	v2 = b;
}


// 格式转换：string转QString
QString s2q(const string &s) 
{  
	return QString(QString::fromLocal8Bit(s.c_str()));  
}  

// 格式转换：QString转string
string q2s(const QString &s) 
{  
	return string((const char *)s.toLocal8Bit());  
}

// 格式转换：wstring转string
std::string ws2s(const std::wstring& ws)
{
	std::string curLocale = setlocale(LC_ALL, NULL);        // curLocale = "C";
	setlocale(LC_ALL, "chs");
	const wchar_t* _Source = ws.c_str();
	size_t _Dsize = 2 * ws.size() + 1;
	char *_Dest = new char[_Dsize];
	memset(_Dest,0,_Dsize);
	wcstombs(_Dest,_Source,_Dsize);
	std::string result = _Dest;
	delete []_Dest;
	setlocale(LC_ALL, curLocale.c_str());
	return result;
}

// 格式转换：string转wstring
std::wstring s2ws(const std::string& s)
{
	setlocale(LC_ALL, "chs"); 
	const char* _Source = s.c_str();
	size_t _Dsize = s.size() + 1;
	wchar_t *_Dest = new wchar_t[_Dsize];
	wmemset(_Dest, 0, _Dsize);
	mbstowcs(_Dest,_Source,_Dsize);
	std::wstring result = _Dest;
	delete []_Dest;
	setlocale(LC_ALL, "C");
	return result;
}

QString int2QString(const int &n)
{
	return QString::number(n);
}

char* cstyleString(QString& str)
{
	char* s = NULL;
	int i;
	s = (char*)malloc((1+str.length()) * sizeof(char));
	for(i = 0 ; i < str.length() ; i++)
		s[i] = str.at(i).toLatin1();
	s[i] = 0;
	return s;
}