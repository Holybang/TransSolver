/** @file 
* @brief  ��ѧ��
* @author �Ժ��
* @date 2008-10-29
* @version 0.1
*
* ��ά����
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
*@param other ָ��
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
*@note ������other��this�໥����  
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
*@note �ӷ�1: this+=other
*@param other ��ӵĶ�ά����
*@return void
*/

inline void Vector2::Add(const Vector2& other)
{
	x += other.x; y += other.y;
}
/**
*@note �ӷ�2: v1+v2
*@param v1 �μ���ӵĶ�ά����
*@param v2 �μ���ӵĶ�ά����
*/
inline void Vector2::Add(const Vector2& v1, const Vector2& v2)
{
	x = v1.x + v2.x; y = v1.y + v2.y;
}

/**
*@note ����1: this-=other
*@param other ������ά����
*/
inline void Vector2::Sub(const Vector2& other)
{
	x -= other.x; y -= other.y;
}
/**
*@note ����2: this=v1-v2        
*@param v1 ��������ά����
*@param v2 ������ά����
*/
inline void Vector2::Sub(const Vector2& v1, const Vector2& v2)
{
	x = v1.x - v2.x; y = v1.y - v2.y;
}
/**
*@note �˷�1��this*=f
*@param f ��������˵���
*/
inline void Vector2::Mul(double f)
{
	x *= f; y *= f;
}
/**
*@note �˷�2��this=other*f 
*@param other ��˵Ķ�ά����
*@param f ��˵���
*/

inline void Vector2::Mul(const Vector2& other, double f)
{
	x = other.x * f; y = other.y * f;
}
/**
*@note ����1��this/=f
*@param f �������������
*/
inline void Vector2::Div(double f)
{
	f = 1.0f / f;
	x *= f; y *= f;
}
/**
*@note  ����2��this=other/f
*param other ��ά����
*param f �������������
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
*@param other ��ά����
*@param s �����任�ĸ�����
*/
inline void Vector2::Combine(const Vector2& other, double s)
{
	x += other.x * s; y += other.y * s;
}
/**
*@note v=v*s+other*t
*@param s ����任�ĸ�����
*@param other ��ά����
*@param t ����任�ĸ�����
*/
inline void Vector2::Combine(double s, const Vector2& other, double t)
{
	x = x * s + other.x * t; y = y * s + other.y * t;
}
/**
*@note v=v*s+other*t 
*@param v1 ����任�Ķ�ά����
*@param v2 ����任�Ķ�ά����
*@param s ����任�ĸ�����
*@param t ����任�ĸ�����
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
*@note ��������ģ
*@return ������ģ
*/
inline double Vector2::Length() const
{
	return (double)sqrt(x * x + y * y);           //����������ģ
}
/**
*@note ��������ģ��ƽ��
*@return ����ģ��ƽ��
*/
inline double Vector2::LengthSqr() const
{
	return x * x + y * y;                         //����������ģ��ƽ��  
}
/**
*@note �����
*@param other ��ά����
*@return ������ά����֮��ľ���
*/
inline double Vector2::Distance(const Vector2& other) const
{
	Vector2 tmp = other - *this;
	return tmp.Length();                        //����sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
/**
*@note ������ƽ�� 
*@param other ��ά����
*@return �����ƽ��
*/
inline double Vector2::DistanceSqr(const Vector2& other) const
{
	Vector2 tmp = other - *this;
	return tmp.LengthSqr();                     //����((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
/**
*@note ��������λ�������򲻱䣬ģΪ1��
*/
inline void Vector2::Normalize()
{
	double f = 1.0f / Length();
	x = x * f; y = y * f;
}
/**
*@note ��other������λ�������򲻱䣬ģΪ1��
*@param other ��ά����
*/
inline void Vector2::Normalize(const Vector2& other)
{
	double f = 1.0f / other.Length();
	x = other.x * f; y = other.y * f;
}
/**
*@note ������˳ʱ����תangle����
*@param angle ��ת�ĽǶ�
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
*@note ��other����˳ʱ����תangle���� 
*@param other ��ά����
*@param angle ��ת�ĽǶ�
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
*@note  Matrix3������Vector2���
*@param mat ��ά����
*@param other ��ά����
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
*@note �����������                                          
*/
inline double Vector2::Dot(const Vector2& vec) const
{
     /*return	| x1  y1|* |x2|
		                   |y2| 
						         */
	return x * vec.x + y * vec.y;
}

