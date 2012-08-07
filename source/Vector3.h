/** @file 
* @brief  ��ѧ��
* @author �Ժ��
* @date 2008-10-29
* @version 0.1
*
* ��ά����
*/

#ifndef _mathematics_h_
#error "Dont include Vector3 directly!!! Use Mathematics.h instead."
#endif

/**
*@note Constructors                                    
*/


inline Vector3::Vector3()
{x = 0; y = 0; z = 0;}
/**
*@note Constructors                                    
*/
inline Vector3::Vector3(double _x, double _y, double _z)
{
	x = _x; y = _y; z = _z;
}
/**
*@note Constructors                                    
*/
inline Vector3::Vector3(const Vector3& other)
{
	x = other.x; y = other.y; z = other.z;
}
/**
*@note Constructors                                    
*/
inline Vector3::Vector3(const double* other)
{
	x = other[0]; y = other[1]; z = other[2];
}

/**
*@note operator  
*/
inline Vector3& Vector3::operator=(const Vector3& other)
{
	x = other.x; y = other.y; z = other.z;
	return *this;
}
/**
*@note assign                                                
*/
inline void Vector3::Assign(double _x, double _y, double _z)
{
	x = _x; y = _y; z = _z;
}
/**
*@note assign                                                
*/
inline void Vector3::Assign(const Vector3& other)
{
	x = other.x; y = other.y; z = other.z;
}
/**
*@note assign                                                
*/
inline void Vector3::Assign(double* other)
{
	x = other[0]; y = other[1]; z = other[2];
}
/**
*@note ������������  
*@param other ��Ҫ��������ά����
*/
inline void Vector3::Swap(Vector3& other)
{
	double _x = x, _y = y, _z = z;
	x = other.x; y = other.y; z = other.z;
	other.x = _x; other.y = _y; other.z = _z;
}
/**
*@note operator                                             
*/
inline Vector3::operator double*()
{
	return &x;
}
/**
*@note operator                                             
*/
inline Vector3::operator const double*() const
{
	return &x;
}
/**
*@note �ӷ�1: this+=other 
*@param other ����
*/
inline void Vector3::Add(const Vector3& other)
{
	x += other.x; y += other.y; z += other.z;
}

/**
*@�ӷ�2: this=v1+v2   
*@param v1 ������ά����
*@param v2 ������ά����
*/
inline void Vector3::Add(const Vector3& v1, const Vector3& v2)
{
	x = v1.x + v2.x; y = v1.y + v2.y; z = v1.z + v2.z;
}
/**
*@note ����1: this-=other
*@param other ������ά����
*/
inline void Vector3::Sub(const Vector3& other)
{
	x -= other.x; y -= other.y; z -= other.z;
}
/**
*@note ����2: this=v1-v2
*@param v1 ������
*@param v2 ����
*/
inline void Vector3::Sub(const Vector3& v1, const Vector3& v2)
{
	x = v1.x - v2.x; y = v1.y - v2.y; z = v1.z - v2.z;
}
/**
*@note �˷�1��this*=f 
*param f ��������˵���
*/

inline void Vector3::Mul(double f)
{
	x *= f; y *= f; z *= f;
}
/**
*@note �˷�2��this=other*f
*@param other ��ά����
*@param f ��������˵���
*/

inline void Vector3::Mul(const Vector3& other, double f)
{
	x = other.x * f; y = other.y * f; z = other.z * f;
}
/**
*@note ����1��this/=f
*@param f �������������
*/
inline void Vector3::Div(double f)
{
	f = 1.0f / f;
	x *= f; y *= f; z *= f;
}
/**
*@note ����2��this=other/f 
*@param other ��λ����
*@param f �������������
*/

inline void Vector3::Div(const Vector3& other, double f)
{
	f = 1.0f / f;
	x = other.x * f; y = other.y * f; z = other.z * f;
}
/**
*@note v=-v 
*/
inline void Vector3::Invert()
{
	x = -x; y = -y; z = -z;
}
/**
*@note v=-other                                           
*/
inline void Vector3::Invert(const Vector3& other)
{
	x = -other.x; y = -other.y; z = -other.z;
}

/**
*@note  v=v+other*s                                            
*/
inline void Vector3::Combine(const Vector3& other, double s)
{
	x += other.x * s; y += other.y * s; z += other.z * s;
}
/**
*@note v=v*s+other*t                                          
*/

inline void Vector3::Combine(double s, const Vector3& other, double t)
{
	x = x * s + other.x * t;
	y = y * s + other.y * t;
	z = z * s + other.z * t;
}

/**
*@note v=v*s+other*t                                          
*/


inline void Vector3::Combine(const Vector3& v1, double s, const Vector3& v2, double t)
{
	x = v1.x * s + v2.x * t;
	y = v1.y * s + v2.y * t;
	z = v1.z * s + v2.z * t;
}

/**
*@note perator                                            
*/

inline Vector3& Vector3::operator+=(const Vector3& other)
{
	x += other.x; y += other.y; z += other.z;
	return *this;
}
/**
*@note perator                                            
*/

inline Vector3& Vector3::operator-=(const Vector3& other)
{
	x -= other.x; y -= other.y; z -= other.z;
	return *this;
}
/**
*@note perator                                            
*/

inline Vector3& Vector3::operator*=(double f)
{
	x *= f; y *= f; z *= f;
	return *this;
}
/**
*@note perator                                            
*/

inline Vector3& Vector3::operator-=(double f)
{
	x -= f; y -= f; z -= f;
	return *this;
}
/**
*@note perator                                            
*/

inline Vector3& Vector3::operator+=(double f)
{
	x += f; y += f; z += f;
	return *this;
}
/**
*@note perator                                            
*/

inline Vector3& Vector3::operator/=(double f)
{
	f = 1.0f / f;
	x *= f; y *= f; z *= f;
	return *this;
}
/**
*@note perator                                            
*/

inline Vector3 Vector3::operator-()
{
	return Vector3(-x, -y, -z);
}
/**
*@note perator                                            
*/

inline Vector3 Vector3::operator+(const Vector3& other) const
{
	return Vector3(x + other.x, y + other.y, z + other.z);
}
/**
*@note perator                                            
*/

inline Vector3 Vector3::operator-(const Vector3& other) const
{
	return Vector3(x - other.x, y - other.y, z - other.z);
}
/**
*@note perator                                            
*/

inline Vector3 Vector3::operator+(double f) const
{
	return Vector3(x + f, y + f, z + f);
}
/**
*@note perator                                            
*/

inline Vector3 Vector3::operator-(double f) const
{
	return Vector3(x - f, y - f, z - f);
}
/**
*@note perator                                            
*/

inline Vector3 Vector3::operator*(double f) const
{
	return Vector3(x * f, y * f, z * f);
}
/**
*@note perator                                            
*/

inline Vector3 Vector3::operator/(double f) const
{
	f = 1.0f / f;
	return Vector3(x * f, y * f, z * f);
}

/**
*@note perator                                            
*/

inline bool Vector3::operator==(const Vector3& other)const
{
	return x==other.x && y == other.y && z ==other.z;
}
/**
*@note perator                                            
*/

inline bool Vector3::operator!=(const Vector3& other)const
{
	return x != other.x || y != other.y || z != other.z;
}


/**
*@note ���������   
*@param v ��˵���ά����
*@return ��˵Ľ��
*/

inline double Vector3::Dot(const Vector3& v) const
{
	return x * v.x + y * v.y + z * v.z;
}
/**
*@note v=v*Matrix3  
*@param v ��ά����
*/

inline void Vector3::Cross(const Vector3& v)
{
	/*
	                             | 0    -v.z    v.y |
	        |x  y  z|=|x  y  z| *|v.z    0      -v.x|
								 |-v.y   v.x      0 |
								                      */
	double _x = y * v.z - z * v.y;
	double _y = z * v.x - x * v.z;
	double _z = x * v.y - y * v.x;
	x = _x; y = _y; z = _z;
}

inline void Vector3::Cross(const Vector3& v1, const Vector3& v2)
{
	/*
	                             | 0    -v.z    v.y |
	        |x  y  z|=|x  y  z| *|v.z    0      -v.x|
	                             |-v.y   v.x      0 |
                                            	       */
	double _x = v1.y * v2.z - v1.z * v2.y;
	double _y = v1.z * v2.x - v1.x * v2.z;
	double _z = v1.x * v2.y - v1.y * v2.x;
	x = _x; y = _y; z = _z;
}
/**
*@note ��������ģ  
*/

inline double Vector3::Length() const
{
	return (double)sqrt(x * x + y * y + z * z);  // ����������ģ
}
/**
*@note ��������ģ��ƽ��                                     
*/

inline double Vector3::LengthSqr() const
{
	return x * x + y * y + z * z;                //����������ģ��ƽ��
}
/**
*@note distance  
*param ��ά����
*/

inline double Vector3::Distance(const Vector3& other) const
{
	Vector3 tmp = other - *this;                //temp=other-this
	return tmp.Length();                        //����temp��ģ
}
/**
*@note distance ��ƽ��  
*@param other ��ά����
*/
inline double Vector3::DistanceSqr(const Vector3& other) const
{
	Vector3 tmp = other - *this;
	return tmp.LengthSqr();
}
/**
*@note ��������λ��                                       
*/

inline void Vector3::Normalize()
{
	if(Length()==0)
		/*throw std::runtime_error("zero vector");*/
		return;
	double f = 1.0f / Length();
	x = x * f; y = y * f; z = z * f;
}
/**
*@note ��������λ��         
*@param other ��Ҫ��λ��������
*/
inline void Vector3::Normalize(const Vector3& other)
{
	double f = 1.0f / other.Length();
	x = other.x * f; y = other.y * f; z = other.z * f;
}

inline void Vector3::Clamp(double min, double max)
{
	x = ::Clamp(x, min, max);
	y = ::Clamp(y, min, max);
	z = ::Clamp(z, min, max);
}

inline void Vector3::Lerp(const Vector3& v1, const Vector3& v2, double s)
{
	x = ::Lerp(s, v1.x, v2.x);
	y = ::Lerp(s, v1.y, v2.y);
	z = ::Lerp(s, v1.z, v2.z);
}

// Operations with a plane

inline void Vector3::Intersect(const Vector3& start, double s, const Vector3& end, double e)
{
	//double ws = s / (e - s);
	//double we = e / (e - s);
	//x = start.x * ws + end.x * we;
	//y = start.y * ws + end.y * we;
	//z = start.z * ws + end.z * we;
	//edited by jinbingwen 2009-7-13 20:34:04
	double ws = s / (e + s);
	double we = e / (e + s);
	x = start.x * we + end.x * ws;
	y = start.y * we + end.y * ws;
	z = start.z * we + end.z * ws;
}

inline void Vector3::Intersect(const Plane& p, const Vector3& start, const Vector3& end)
{
	double s = p.Distance(start);
	double e = p.Distance(end);
	Intersect(start, s, end, e);
}

inline void Vector3::Reflect(const Plane& p, const Vector3& v)
{
	Vector3 tmp = p.n * 2 * v.Dot(p.n) - v;
	x = tmp.x; y = tmp.y; z = tmp.z;
}

// Operations with a matrix

inline void Vector3::Transform(const Matrix4& mat)
{
	const double* m = mat.M;
	double fx = x; double fy = y; double fz = z;
	x = m[0] * fx + m[4] * fy + m[8] * fz + m[12];  // w assumed to be 1.0F
	y = m[1] * fx + m[5] * fy + m[9] * fz + m[13];
	z = m[2] * fx + m[6] * fy + m[10] * fz + m[14];
}

inline void Vector3::Transform(const Matrix4& mat, const Vector3& other)
{
	const double* m = mat.M;
	double fx = other.x; double fy = other.y; double fz = other.z;
	x = m[0] * fx + m[4] * fy + m[8] * fz + m[12];  // w assumed to be 1.0F
	y = m[1] * fx + m[5] * fy + m[9] * fz + m[13];
	z = m[2] * fx + m[6] * fy + m[10] * fz + m[14];
}

inline void Vector3::TransformTranspose(const Matrix4& mat)
{
	const double* m = mat.M;
	double fx = x; double fy = y; double fz = z;
	x = m[0] * fx + m[1] * fy + m[2] * fz + m[3];  // w assumed to be 1.0F
	y = m[4] * fx + m[5] * fy + m[6] * fz + m[7];
	z = m[8] * fx + m[9] * fy + m[10] * fz + m[11];
}

inline void Vector3::TransformTranspose(const Matrix4& mat, const Vector3& other)
{
	const double* m = mat.M;
	double fx = other.x; double fy = other.y; double fz = other.z;
	x = m[0] * fx + m[1] * fy + m[2] * fz + m[3];  // w assumed to be 1.0F
	y = m[4] * fx + m[5] * fy + m[6] * fz + m[7];
	z = m[8] * fx + m[9] * fy + m[10] * fz + m[11];
}

// Other

inline Vector3 Vector3::Scaling(const Vector3& v)
{
	Vector3 vec;
	vec.x = x * v.x;
	vec.y = y * v.y;
	vec.z = z * v.z;
	return vec;
}
