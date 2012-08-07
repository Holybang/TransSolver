/** @file 
* @brief  数学库
* @author 赵洪德
* @date 2008-10-29
* @version 0.1
*
* 阐述数学库中的类
*/






#ifndef _mathematics_h_
#define _mathematics_h_

#include <cmath>
#include <stdlib.h> 
#include <iostream>
#include <iomanip>

//using namespace std;

// ------------------------------------------------------------
// constants
// ------------------------------------------------------------

// PI & Co.
#ifndef M_PI
const double M_PI				= 3.14159265358979323846f;
#endif
const double M_DEG_TO_RAD		= 0.01745329251994329547f;
const double M_RAD_TO_DEG		= 57.29577951308232286465f;
const double M_PI_DIV_2			= 1.570796326794896558f;
const double M_PI_DIV_4			= 0.785398163397448279f;

// tolerance values

const double M_LOW_TOLERANCE		= 0.000001f;
const double M_TOLERANCE			= 0.000010f;
const double M_HIGH_TOLERANCE	= 0.000100f;

// ------------------------------------------------------------
// basic functions
// ------------------------------------------------------------

/**
*@note 求取f,min, max中位于中间的数   
*@param f 参数一
*@param min 参数二
*@param max 参数三
*@return 中间值
*/
inline double Clamp(double f, double min, double max)
{
	return min > f ? min : (max > f ? f : max);
}

inline double Lerp(double f, double from, double to)
{
	return from + (to - from) * f;
}
/**
*@note 最小值函数 
*@param f1 参数一
*@param f2 参数二
*@return 最小值
*/
inline double Min(double f1, double f2)
{
	return f1 < f2 ? f1 : f2;
}
/**
*@note 最大值函数 
*@param f1 参数一
*@param f2 参数二
*@return 最大值
*/
inline double Max(double f1, double f2)
{
	return f1 > f2 ? f1 : f2;
}
/**
*@note 求角度f的正弦值                              
*@param f 角度
*@return 正弦值
*/
inline double SinD(double f)
{
	return (double)sin(f * M_DEG_TO_RAD);  //变量f为角度值，f * M_DEG_TO_RAD转化为弧度值
}
/**
*@note 求数值的反正弦值                                      
*@param fValue 数值
*@return 反正弦值
*/
inline double ASinD(double fValue)
{
	if (-1.0 < fValue)
	{
		if ( fValue < 1.0 )
			return (asin(fValue)*M_RAD_TO_DEG);     //返回角度，asin(fValue)*M_RAD_TO_DEG是转化为角度值
		else
			return -90.0f;                          //exception
	}
	else
	{
		return 90.0f;                                //exception
	}
}
/**
*@note 求角度f的余弦值                              
*@param f 角度
*@return 余弦值
*/
inline double CosD(double f)
{
	return (double)cos(f * M_DEG_TO_RAD);            //f为角度，f*M_DEG_TO_RAD转化为弧度值
}
/**
*@note 求数值的反余弦值                                      
*@param fValue 数值
*@return 反余弦值
*/
inline double ACosD(double fValue)
{
	if ( -1.0 < fValue )
	{
		if ( fValue < 1.0 )
			return (acos(fValue)*M_RAD_TO_DEG);       //返回角度，acos(fValue)*M_RAD_TO_DEG是转化为角度值
		else
			return 0.0;                              //exception
	}
	else
	{
		return 180.0;     ;                                 //exception
	}
}

/**
*@note 求角度f的正切值                              
*@param f 角度
*@return 正切值
*/
inline double TanD(double f)
{
	return (double)tan(f * M_DEG_TO_RAD);
}

/**
*@note 求数值的反正切值                                      
*@param fValue 数值
*@return 反正切值
*/
inline double ATanD(double fValue)
{
	return (atan(fValue) * M_RAD_TO_DEG);
}
/**
*@note 求数值f1/f2的反正切值
*@param f1 
*@param f2
*@return 正切值
*/
inline double ATanD2 (double f1, double f2)
{
	return (atan2( f1, f2) * M_RAD_TO_DEG);
}

/**
*@note 将角度f转化为弧度 
*@param f 角度
*@return 弧度值
*/
inline double Degree(double f)
{
	return f * M_DEG_TO_RAD;
}
/**
*@note 将弧度f转化为角度
*@param f 弧度值
*@return 角度值
*/
inline double Radians(double f)
{
	return f * M_RAD_TO_DEG;
}

inline bool IsZeroL(double f)
{
	return fabs(f) > M_LOW_TOLERANCE;
}

inline bool IsZero(double f)
{
	return fabs(f) > M_TOLERANCE;
}

inline bool IsZeroH(double f)
{
	return fabs(f) > M_HIGH_TOLERANCE;
}
/**
*@note 生成随机数   
*@return 随机数
*/
inline double Rand(double max)
{
	return (double(rand()) / double(RAND_MAX)) * max;    //    0<=rand<=RAND_MAX
}
/**
*@note 生成（min, max）之间的一个随机数
*@param min 最小数
*@param max 最大数
*@return 随机数
*/
inline double Rand(double min, double max)
{
	return Rand(max - min) + min;                       //  生成（min, max）之间的一个随机数   
}

/**
*@note 取d的绝对值
*@return d的绝对值
*/
inline double AbsDouble(double d)
{
	if(d < 0){
		return -d;
	} else{
		return d;
	}
}

enum EulerAnglesSeq{
	SEQ_ZXY,
	SEQ_XYZ,

	SEQ_ERROR
};


class Vector2;
class Matrix3;
class Vector3;
class Plane;
class Matrix4;
class Quaternion;



class Vector2
{
public:

	// Members

	double x, y;

public:

	// Constructors

	Vector2();
	Vector2(const double* other);
	Vector2(double _x, double _y);
	Vector2(const Vector2& other);

	// Assign

	Vector2& operator=(const Vector2& other);
	void Assign(double _x, double _y);
	void Assign(const Vector2& other);
	void Assign(double* other);

	// Swap

	void Swap(Vector2& other);

	// Casting

	operator double*();
	operator const double*() const;

	// Math

	// this += other
	void Add(const Vector2& other);
	// this = v1 + v2
	void Add(const Vector2& v1, const Vector2& v2);
	// this -= other
	void Sub(const Vector2& other);
	// this = v1 - v2
	void Sub(const Vector2& v1, const Vector2& v2);
	// this *= f
	void Mul(double f);
	// this = other * f
	void Mul(const Vector2& other, double f);
	// this /= f
	void Div(double f);
	// this = other / f
	void Div(const Vector2& other, double f);
	// this = -this
	void Invert();
	// this = -other
	void Invert(const Vector2& other);

	// Linear combination

	// this += other * s
	void Combine(const Vector2& other, double s);
	// this = this * s + other * t
	void Combine(double s, const Vector2& other, double t);
	// this = v1 * s + v2 * t
	void Combine(const Vector2& v1, double s, const Vector2& v2, double t);

	// operators

	Vector2& operator+=(const Vector2& other);
	Vector2& operator-=(const Vector2& other);
	Vector2& operator*=(double f);
	Vector2& operator/=(double f);
	Vector2& operator+=(double f);
	Vector2& operator-=(double f);
	Vector2 operator-();
	Vector2 operator+(const Vector2& other) const;
	Vector2 operator-(const Vector2& other) const;
	Vector2 operator*(double f) const;
	Vector2 operator/(double f) const;
	Vector2 operator+(double f) const;
	Vector2 operator-(double f) const;

	bool operator==(const Vector2& other);
	bool operator!=(const Vector2& other);

	// Linear Algebra

	double Length() const;
	double LengthSqr() const;
	double Distance(const Vector2& other) const;
	double DistanceSqr(const Vector2& other) const;
	void Normalize();
	void Normalize(const Vector2& other);

	void Rotate(double angle);
	void Rotate(Vector2 other, double angle);

	// Others

	void Clamp(double min, double max);
	// this = v1 + (v2 - v1) * s
	void Lerp(const Vector2& v1, const Vector2& v2, double s);

	// Actions with Matrix3

	void Transform(const Matrix3& mat, const Vector2& other);
	void Transform(const Matrix3& mat);

	double Dot(const Vector2& vec) const;

	// friend operators
	friend std::ostream& operator<<(std::ostream& os, Vector2& v);

};


class Matrix3
{
public:

	// members

	union
	{
		double M[9];
		double m[3][3];
	};

	static double Identity_Matrix[];

public:

	// Constructors

	Matrix3();
	Matrix3(const double* other);
	Matrix3(const Matrix3& other);
	Matrix3(const double f0, const double f1, const double f2,
			const double f3, const double f4, const double f5,
			const double f6, const double f7, const double f8);

	// Assignment

	Matrix3& operator=(const Matrix3& other);
	void Assign(const double* other);
	void Assign(const Matrix3& other);

	// Casting

	operator double*();
	operator const double*() const;

	// Simple methods

	Matrix3& SetIdentity();
	Matrix3& SetZero();

	// this *= other
	Matrix3& Mul(const Matrix3& other);
	// this = m1 * m2
	Matrix3& Mul(const Matrix3& m1, const Matrix3& m2);
	// this *= f
	Matrix3& Mul(double f);
	// this = other * f
	Matrix3& Mul(const Matrix3& other, double f);

	/************************************************************************/
	/* Added by Jin Bingwen                                                 */
	/************************************************************************/
	Vector3 Mul(const Vector3 v);

	// this /= f
	Matrix3& Div(double f);
	// this = other / f
	Matrix3& Div(const Matrix3& other, double f);

	void Transpose();

	/************************************************************************/
	/* Add by Liu Lu                                                        */
	/************************************************************************/

	Matrix3& Add(Matrix3& other);
	Matrix3& Add( Matrix3 &M1, Matrix3 &M2 );
	Matrix3& Subtract(Matrix3& other);
	Matrix3& Subtract( Matrix3 &M1, Matrix3 &M2 );

	// Dst = this * v
	void PostMultiply( Vector3 &v, Vector3 &Dst );
	// Dst = v*this
	void PreMultiply( Vector3 &v, Vector3 &Dst );

	bool IsTranspose( Matrix3 &other );

	void SetFromOuterProduct( Vector3 &v1, Vector3 &v2 );

	Matrix3& Invert(const Matrix3& other);

	Matrix3& CofactorMat(const Matrix3& other);

	double Determinant(const Matrix3& other);

	/************************************************************************/
	/* End                                                                  */
	/************************************************************************/

	// Operators

	Matrix3& operator*=(const Matrix3& other);
	Matrix3 operator*(const Matrix3& other) const;

	// Transformation methods

	Matrix3& Translate(double tx, double ty);
	Matrix3& SetTranslate(double tx, double ty);
	Matrix3& Rotate(double r);
	Matrix3& SetRotate(double r);
	Matrix3& Scale(double sx, double sy);
	Matrix3& SetScale(double sx, double sy);

	// Other
	void FromEulerAnglesZYX (const double& fYAngle, const double& fPAngle, const double& fRAngle);
	void FromEulerAnglesZXY (const Vector3& vec);
	Vector3 ToEulerAnglesZXY (void);

	// friends
	friend std::ostream& operator<<(std::ostream& os, Matrix3 &m);

};


class Vector3
{
public:



	// Member variables

	double x, y, z;

	// constant
	static const Vector3 ZERO;
	static const Vector3 UNIT_SCALE;
	static const Vector3 UNIT_X;
	static const Vector3 UNIT_Y;
	static const Vector3 UNIT_Z;
	static const Vector3 NEGATIVE_UNIT_X;
	static const Vector3 NEGATIVE_UNIT_Y;
	static const Vector3 NEGATIVE_UNIT_Z;

public:

	// Constructors

	Vector3();
	Vector3(double _x, double _y, double _z);
	Vector3(const Vector3& other);
	Vector3(const double* other);

	// assignemt

	Vector3& operator=(const Vector3& other);
	void Assign(double _x, double _y, double _z);
	void Assign(const Vector3& other);
	void Assign(double* other);

	// Swap

	void Swap(Vector3& other);

	// Casting

	operator double*();
	operator const double*() const;

	// Math

	// this += other
	void Add(const Vector3& other);
	// this = v1 + v2
	void Add(const Vector3& v1, const Vector3& v2);
	// this -= other
	void Sub(const Vector3& other);
	// this = v1 - v2
	void Sub(const Vector3& v1, const Vector3& v2);
	// this *= f
	void Mul(double f);
	// this = other * f
	void Mul(const Vector3& other, double f);
	// this /= f
	void Div(double f);
	// this = other / f
	void Div(const Vector3& other, double f);
	// this = -this
	void Invert();
	// this = -other
	void Invert(const Vector3& other);

	// Linear combination

	// this += other * s
	void Combine(const Vector3& other, double s);
	// this = this * s + other * t
	void Combine(double s, const Vector3& other, double t);
	// this = v1 * s + v2 * t
	void Combine(const Vector3& v1, double s, const Vector3& v2, double t);

	// EvilOne's beloved overloads

	Vector3& operator+=(const Vector3& other);
	Vector3& operator-=(const Vector3& other);
	Vector3& operator*=(double f);
	Vector3& operator/=(double f);
	Vector3& operator+=(double f);
	Vector3& operator-=(double f);
	Vector3 operator-();
	Vector3 operator+(const Vector3& other) const;
	Vector3 operator-(const Vector3& other) const;
	Vector3 operator*(double f) const;
	Vector3 operator/(double f) const;
	Vector3 operator+(double f) const;
	Vector3 operator-(double f) const;

	bool operator==(const Vector3& other)const;
	bool operator!=(const Vector3& other)const;
	// Linear Algebra

	double Dot(const Vector3& v) const;
	// this = v1 x v2
	void Cross(const Vector3& v);
	void Cross(const Vector3& v1, const Vector3& v2);
	double Length() const;
	double LengthSqr() const;
	double Distance(const Vector3& other) const;
	double DistanceSqr(const Vector3& other) const;
	void Normalize();
	void Normalize(const Vector3& other);

	// Other usefull funtions

	void Clamp(double min, double max);
	// this = v1 + (v2 - v1) * s
	void Lerp(const Vector3& v1, const Vector3& v2, double s);

	// Operations with a plane

	void Intersect(const Vector3& start, double s, const Vector3& end, double e);
	void Intersect(const Plane& p, const Vector3& start, const Vector3& end);
	void Reflect(const Plane& p, const Vector3& v);

	// Operations with Matrix

	void Transform(const Matrix4& mat, const Vector3& other);
	void Transform(const Matrix4& mat);
	void TransformTranspose(const Matrix4& mat, const Vector3& other);
	void TransformTranspose(const Matrix4& mat);

	// Other 

	Vector3 Scaling(const Vector3& v);
	// friend operators
	friend std::ostream& operator<<(std::ostream& os, Vector3& v);

};

// ------------------------------------------------------------
// Plane
// ------------------------------------------------------------

class Plane
{
public:

	// members 

	union
	{
		struct
		{
			double a, b, c, d;
		};

		struct
		{
			Vector3 n;
			double d;
		};
	};

public:

	// Constructors

	Plane();
	Plane(const Vector3& _n, double _d);
	Plane(double _a, double _b, double _c, double _d);
	Plane(double* other);
	Plane(const Plane& other);

	// Assignment

	Plane& operator=(const Plane& other);
	void Assign(const Vector3& _n, double _d);
	void Assign(double _a, double _b, double _c, double _d);
	void Assign(double* other);
	void Assign(const Plane& other);

	// Functions

	double Distance(const Vector3& v) const;

	// Interaction with matrix

	void Transform(const Matrix4 &mat);
	void Transform(const Matrix4 &mat, const Plane& other);

	// friend operators
	friend std::ostream& operator<<(std::ostream& os, Plane& p);

};

// ------------------------------------------------------------
// Matrix4
// ------------------------------------------------------------

class Matrix4
{
public:

	// members

	union
	{
		double M[16];
		double m[4][4];
	};

	static double Identity_Matrix[];
	static double Orientation_Switch_Matrix[];
	static double Perspective_Matrix[];

public:

	// Constructors

	Matrix4();
	Matrix4(const double* other);
	Matrix4(const Matrix4& other);
	Matrix4(const Quaternion& other);

	// Assignment

	Matrix4& operator=(const Matrix4& other);
	Matrix4& operator=(const Quaternion& other);
	void Assign(const Quaternion& other);
	void Assign(const double* other);
	void Assign(const Matrix4& other);

	// Casting

	operator double*();
	operator const double*() const;

	// Simple methods

	Matrix4& SetIdentity();
	Matrix4& SetZero();
	Matrix4& SetPerspective();
	Matrix4& SetSwitchOrientation();

	// this *= other
	Matrix4& Mul(const Matrix4& other);
	// this = m1 * m2
	Matrix4& Mul(const Matrix4& m1, const Matrix4& m2);
	// this *= f
	Matrix4& Mul(double f);
	// this = other * f
	Matrix4& Mul(const Matrix4& other, double f);
	// this /= f
	Matrix4& Div(double f);
	// this = other / f
	Matrix4& Div(const Matrix4& other, double f);

	//add by frankmao, matrix * vector
	Vector3 Mul(const Vector3& v);  

	// Operators

	Matrix4& operator*=(const Matrix4& other);
	Matrix4 operator*(const Matrix4& other) const;

	// Transformation methods

	Matrix4& SetTranslate(double tx, double ty, double tz);
	Matrix4& Translate(double tx, double ty, double tz);
	Matrix4& SetScale(double sx, double sy, double sz);
	Matrix4& Scale(double sx, double sy, double sz);

	// rotation around three euler-angles

	Matrix4& SetRotate(const Vector3& r);
	Matrix4& SetRotate(double rx, double ry, double rz);
	Matrix4& Rotate(const Vector3& r);
	Matrix4& Rotate(double rx, double ry, double rz);

	// rotation euler-angle around axis

	Matrix4& SetRotate(double angle, const Vector3& r);
	Matrix4& SetRotate(double angle, double x, double y, double z);
	Matrix4& Rotate(double angle, const Vector3& r);
	Matrix4& Rotate(double angle, double x, double y, double z);

	// Invert/Transpose

	Matrix4& Adjoint();
	Matrix4& Adjoint(const Matrix4& other);

	Matrix4& Transpose();
	Matrix4& Transpose(const Matrix4& other);

	Matrix4& Invert();
	Matrix4& Invert(const Matrix4& other);

	double Det() const;

	// Perpsective

	Matrix4& SetProjection(double fov, double aspect, double znear, double zfar);
	Matrix4& SetOthogonal(double znear, double zfar);

	// static

	static double Det2x2(double a, double b, double c, double d);
	static double Det3x3(double a1, double a2, double a3,
		double b1, double b2, double b3,
		double c1, double c2, double c3);

	// friends
	friend std::ostream& operator<<(std::ostream& os, Matrix4 &m);

};

// ------------------------------------------------------------
// Quaternion
// ------------------------------------------------------------

class Quaternion
{
public:

	// constant
	static const Quaternion ZERO;
	static const Quaternion IDENTITY;

	// Members
	double x, y, z, w;

public:

	// Constructors

	Quaternion();
	Quaternion(const Quaternion& other);
	Quaternion(const double* other);
	Quaternion(double _x, double _y, double _z, double _w);

	Quaternion(const Vector3& xAxis, const Vector3& yAxis, const Vector3& zAxis);
	// Assignment

	Quaternion& operator=(const Quaternion& other);
	void Assign(const Quaternion& other);
	void Assign(const double* other);
	void Assign(double _x, double _y, double _z, double _w);

	// Simple math

	void Invert();
	void Invert(const Quaternion& other);
	Quaternion operator- () const;
	Quaternion& Mul(const Quaternion& other);
	Quaternion& Mul(const Quaternion& q1, const Quaternion& q2);
	Quaternion operator*(const Quaternion& other) const;
	Quaternion& operator*=(const Quaternion& other);


	double Length() const;
	double LengthSqr() const;

	void Normalize();
	void Normalize(const Quaternion& other);

	// Rotation around euler-angle and axis

	Quaternion& SetRotate(double angle, const Vector3& v);
	Quaternion& Rotate(double angle, const Vector3& v);
	Quaternion& SetRotate(double angle, double _x, double _y, double _z);
	Quaternion& Rotate(double angle, double _x, double _y, double _z);

	


	// rotation around three euler-angles

	Quaternion& SetRotate(const Vector3& v);
	Quaternion& Rotate(const Vector3& v);
	Quaternion& SetRotate(double rx, double ry, double rz);
	Quaternion& Rotate(double rx, double ry, double rz);

	Quaternion& FromEulerAngleZXY(double rz, double rx, double ry);
	Vector3 ToEulerAngleZXY(void)const;
	Quaternion& FromEulerAngleXYZ(double rx, double ry, double rz);
	Vector3 ToEulerAngleXYZ(void)const;

	Quaternion& FromEulerAngle(double rx, double ry, double rz, EulerAnglesSeq seq = SEQ_ZXY);
	Vector3 ToEulerAngle(EulerAnglesSeq seq = SEQ_ZXY) const;

	void FromAngleAxis(const double& fAngle, const Vector3& kAxis) ;
	void ToAngleAxis (double& rfAngle, Vector3& rkAxis) const;
	void ToRotationMatrix (Matrix3& kRot) const;

	void FromRotationMatrix(const Matrix3& kRot);
	void FromAxes (const Vector3& xaxis, const Vector3& yaxis, const Vector3& zaxis);


	double Dot (const Quaternion& rkQ) const;
	double Norm () const;
	Quaternion operator* (double fScalar) const;
	Quaternion Inverse () const;
	static Quaternion nlerp(double fT, const Quaternion& rkP,const Quaternion& rkQ, bool shortestPath);
	Quaternion operator- (const Quaternion& rkQ) const;
	Quaternion TransitionFromCurrentTo(const Quaternion rkQ);
	Quaternion operator+ (const Quaternion& rkQ) const;
	Quaternion Log () const;

	// Interpolation

	void Lerp(const Quaternion& q0, const Quaternion& q1, double t);
	void Slerp(const Quaternion& q0, const Quaternion& q1, double t);

	inline bool operator== (const Quaternion& rhs) const
	{
		return (rhs.x == x) && (rhs.y == y) &&
			(rhs.z == z) && (rhs.w == w);
	}

	// operator with vector

	//Quaternion Mul (const Quaternion& q, const Vector3& v) const;
	//Quaternion Mul (const Vector3& v, const Quaternion& q) const;
	Vector3 operator*(const Vector3& v) const;

	//Vector3 QVRotation(const Quaternion& q, const Vector3& v) const;

	// friend operators
	friend std::ostream& operator<<(std::ostream& os, Quaternion& v);

};

// ------------------------------------------------------------
// Include implementations
// ------------------------------------------------------------

#include "Vector2.h"
#include "Matrix3.h"
#include "Vector3.h"
#include "Plane.h"
#include "Matrix4.h"
#include "Quaternion.h"

#endif
