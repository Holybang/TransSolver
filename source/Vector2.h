/** @file 
* @brief  数学库
* @author 赵洪德
* @date 2008-10-29
* @version 0.1
*
* 二维向量
*/

#ifndef _mathematics_h_
#error "Dont include Vector2 directly!!! Use Mathematics.h instead."
#endif
/**
*@constructor
*/
inline Vector2::Vector2()
: x(0)
, y(0)
{}
/**
*@note constructor
*@param other 指针
*/

inline Vector2::Vector2(const double* other)
{
	x = other[0]; y = other[1];
}
/**
*@note constructor
*@param _x
*@param _y
*/

inline Vector2::Vector2(double _x, double _y)
{
	x = _x; y = _y;
}
/**
*@constructor
*@param other vector2
*/
inline Vector2::Vector2(const Vector2& other)
{
	x = other.x; y = other.y;
}
/**
*@note operator    
*/
inline Vector2& Vector2::operator=(const Vector2& other)
{
	x = other.x; y = other.y;
	return *this;
}
/**
*@note Assign 
*/
inline void Vector2::Assign(double _x, double _y)
{
	x = _x; y = _y;
}


inline void Vector2::Assign(const Vector2& other)
{
	x = other.x; y = other.y;
}

inline void Vector2::Assign(double* other)
{
	x = other[0]; y = other[1];
}

/*
*@note 将向量other和this相互交换  
*@param other Vector2
*@return void 
*/

inline void Vector2::Swap(Vector2& other)
{
	Vector2 tmp = *this;
	Assign(other);
	other.Assign(tmp);
}
/*
*@note operator
*/
inline Vector2::operator double*()
{
	return &x;
}
/**
*@note operator
**/
inline Vector2::operator const double*() const
{
	return &x;
}

/**
*@note 加法1: this+=other
*@param other 相加的二维向量
*@return void
*/

inline void Vector2::Add(const Vector2& other)
{
	x += other.x; y += other.y;
}
/**
*@note 加法2: v1+v2
*@param v1 参加相加的二维向量
*@param v2 参加相加的二维向量
*/
inline void Vector2::Add(const Vector2& v1, const Vector2& v2)
{
	x = v1.x + v2.x; y = v1.y + v2.y;
}

/**
*@note 减法1: this-=other
*@param other 减数二维向量
*/
inline void Vector2::Sub(const Vector2& other)
{
	x -= other.x; y -= other.y;
}
/**
*@note 减法2: this=v1-v2        
*@param v1 被减数二维向量
*@param v2 减数二维向量
*/
inline void Vector2::Sub(const Vector2& v1, const Vector2& v2)
{
	x = v1.x - v2.x; y = v1.y - v2.y;
}
/**
*@note 乘法1：this*=f
*@param f 与向量相乘的数
*/
inline void Vector2::Mul(double f)
{
	x *= f; y *= f;
}
/**
*@note 乘法2：this=other*f 
*@param other 相乘的二维向量
*@param f 相乘的数
*/

inline void Vector2::Mul(const Vector2& other, double f)
{
	x = other.x * f; y = other.y * f;
}
/**
*@note 除法1：this/=f
*@param f 与向量相除的数
*/
inline void Vector2::Div(double f)
{
	f = 1.0f / f;
	x *= f; y *= f;
}
/**
*@note  除法2：this=other/f
*param other 二维向量
*param f 与向量相除的数
*/
inline void Vector2::Div(const Vector2& other, double f)
{
	f = 1.0f / f;
	x = other.x * f; y = other.y * f;
}
/**
*@note v=-v
*/
inline void Vector2::Invert()
{
	x = -x; y = -y;
}
/**
*@note v=-other
*/
inline void Vector2::Invert(const Vector2& other)
{
	x = -other.x; y = -other.y;
}

/**
*@note v=v+other*s 
*@param other 二维向量
*@param s 参数变换的浮点数
*/
inline void Vector2::Combine(const Vector2& other, double s)
{
	x += other.x * s; y += other.y * s;
}
/**
*@note v=v*s+other*t
*@param s 参与变换的浮点数
*@param other 二维向量
*@param t 参与变换的浮点数
*/
inline void Vector2::Combine(double s, const Vector2& other, double t)
{
	x = x * s + other.x * t; y = y * s + other.y * t;
}
/**
*@note v=v*s+other*t 
*@param v1 参与变换的二维向量
*@param v2 参与变换的二维向量
*@param s 参与变换的浮点数
*@param t 参与变换的浮点数
*/
inline void Vector2::Combine(const Vector2& v1, double s, const Vector2& v2, double t)
{
	x = v1.x * s + v2.x * t; y = v1.y * s + v2.y * t;
}

/**
*@note operator 
*/
inline Vector2& Vector2::operator+=(const Vector2& other)
{
	x += other.x; y += other.y;
	return *this;
}
/**
*@note operator 
*/
inline Vector2& Vector2::operator-=(const Vector2& other)
{
	x -= other.x; y -= other.y;
	return *this;
}
/**
*@note operator 
*/
inline Vector2& Vector2::operator*=(double f)
{
	x *= f; y *= f;
	return *this;
}
/**
*@note operator 
*/
inline Vector2& Vector2::operator-=(double f)
{
	x -= f; y -= f;
	return *this;
}
/**
*@note operator 
*/
inline Vector2& Vector2::operator+=(double f)
{
	x += f; y += f;
	return *this;
}
/**
*@note operator 
*/
inline Vector2& Vector2::operator/=(double f)
{
	f = 1.0f / f;
	x *= f; y *= f;
	return *this;
}
/**
*@note operator 
*/
inline Vector2 Vector2::operator-()
{	
	return Vector2(-x, -y);
}
/**
*@note operator 
*/
inline Vector2 Vector2::operator+(const Vector2& other) const
{
	return Vector2(x + other.x, y + other.y);
}
/**
*@note operator 
*/
inline Vector2 Vector2::operator-(const Vector2& other) const
{
	return Vector2(x - other.x, y - other.y);
}

//inline Vector2 Vector2::operator+(double f) const
//{
//	return Vector2(x + f, y + f);
//}
/**
*@note operator 
*/
inline Vector2 Vector2::operator-(double f) const
{
	return Vector2(x - f, y - f);
}
/**
*@note operator 
*/
inline Vector2 Vector2::operator*(double f) const
{
	return Vector2(x * f, y * f);
}
/**
*@note operator 
*/
inline Vector2 Vector2::operator/(double f) const
{
	f = 1.0f / f;
	return Vector2(x * f, y * f);
}

inline bool Vector2::operator==( const Vector2& other )
{
	return x == other.x && y == other.y;
}

inline bool Vector2::operator!=(const Vector2& other )
{
	return x != other.x || y != other.y;
}

/**
*@note 求向量的模
*@return 向量的模
*/
inline double Vector2::Length() const
{
	return (double)sqrt(x * x + y * y);           //返回向量的模
}
/**
*@note 求向量的模的平方
*@return 向量模的平方
*/
inline double Vector2::LengthSqr() const
{
	return x * x + y * y;                         //返回向量的模的平方  
}
/**
*@note 求距离
*@param other 二维向量
*@return 两个二维向量之间的距离
*/
inline double Vector2::Distance(const Vector2& other) const
{
	Vector2 tmp = other - *this;
	return tmp.Length();                        //返回sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
/**
*@note 求距离的平方 
*@param other 二维向量
*@return 距离的平方
*/
inline double Vector2::DistanceSqr(const Vector2& other) const
{
	Vector2 tmp = other - *this;
	return tmp.LengthSqr();                     //返回((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
/**
*@note 将向量单位化（方向不变，模为1）
*/
inline void Vector2::Normalize()
{
	double f = 1.0f / Length();
	x = x * f; y = y * f;
}
/**
*@note 将other向量单位化（方向不变，模为1）
*@param other 二维向量
*/
inline void Vector2::Normalize(const Vector2& other)
{
	double f = 1.0f / other.Length();
	x = other.x * f; y = other.y * f;
}
/**
*@note 将向量顺时针旋转angle弧度
*@param angle 旋转的角度
*/
inline void Vector2::Rotate(double angle)
{
	/*
	    |x|  |  CosD(angle) SinD(angle) |  |x|
		|y|= | -SinD(angle) CosD(angle) | *|y|
		                                       */
	double _x = x;
	x =  _x * CosD(angle) + y * SinD(angle);
	y = -_x * SinD(angle) + y * CosD(angle);
}
/**
*@note 将other向量顺时针旋转angle弧度 
*@param other 二维向量
*@param angle 旋转的角度
*/
inline void Vector2::Rotate(Vector2 other, double angle)
{
	/*
	  |x|  | CosD(angle) SinD(angle) |  |x|
	  |y|= | SinD(angle) CosD(angle) | *|y|
	                                        */
	double _x = other.x;
	double _y = other.y;

	x =  _x * CosD(angle) + _y * SinD(angle);
	y = -_x * SinD(angle) + _y * CosD(angle);
}
inline void Vector2::Clamp(double min, double max)
{
	x = ::Clamp(x, min, max); y = ::Clamp(y, min, max);
}

inline void Vector2::Lerp(const Vector2& v1, const Vector2& v2, double s)
{
	x = ::Lerp(s, v1.x, v2.x); y = ::Lerp(s, v1.y, v2.y);
}
/**
*@note  Matrix3矩阵与Vector2相乘
*@param mat 三维矩阵
*@param other 二维向量
*/
inline void Vector2::Transform(const Matrix3& mat, const Vector2& other)
{
	const double* m = mat.M;
	double fx = other.x; double fy = other.y;
	x = m[0] * fx + m[3] * fy + m[6];  // w assumed to be 1.0F
	y = m[1] * fx + m[4] * fy + m[7];
}

inline void Vector2::Transform(const Matrix3& mat)
{
	const double* m = mat.M;
	double fx = x; double fy = y;
	x = m[0] * fx + m[3] * fy + m[6];  // w assumed to be 1.0F
	y = m[1] * fx + m[4] * fy + m[7];
}
/**
*@note 两个向量相乘                                          
*/
inline double Vector2::Dot(const Vector2& vec) const
{
     /*return	| x1  y1|* |x2|
		                   |y2| 
						         */
	return x * vec.x + y * vec.y;
}

