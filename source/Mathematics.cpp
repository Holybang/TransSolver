/** @file 
* @brief  数学库
* @author 赵洪德
* @date 2008-10-29
* @version 0.1
*
* 数学库的方法
*/


#include "Mathematics.h"

const Quaternion Quaternion::ZERO(0.0f, 0.0f, 0.0f, 0.0f);
const Quaternion Quaternion::IDENTITY(0.0f, 0.0f, 0.0f, 1.0f);
const Vector3 Vector3::ZERO(0.0f, 0.0f, 0.0f);
const Vector3 Vector3::UNIT_SCALE(1.0f, 1.0f, 1.0f);
const Vector3 Vector3::UNIT_X(1.0f, 0, 0);
const Vector3 Vector3::UNIT_Y(0, 1.0f, 0);
const Vector3 Vector3::UNIT_Z(0, 0, 1.0f);
const Vector3 Vector3::NEGATIVE_UNIT_X(-1.0f, 0, 0);
const Vector3 Vector3::NEGATIVE_UNIT_Y(0, -1.0f, 0);
const Vector3 Vector3::NEGATIVE_UNIT_Z(0, 0, -1.0f);
std::ostream& operator<<(std::ostream& os, Vector2& v)
{
	os << "(" << std::setprecision(4) << v.x << ", "  << v.y << ")";
	return os;
}


std::ostream& operator<<(std::ostream& os, Vector3& v)
{
	os << "(" << std::setprecision(4) << v.x << ", "  << v.y << ", "  << v.z << ")";
	return os;
}


std::ostream& operator<<(std::ostream& os, Plane& p)
{
	os << "(" << std::setprecision(4) << p.a << ", " 
		<< p.b << ", "  << p.c << ", " << p.d << ")";
	return os;
}
/**
*@note  定义Identity_Matrix为单位矩阵                       
*/


double Matrix3::Identity_Matrix[] =
{
	1.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 1.0f
};

void Matrix3::Transpose()
{
	Matrix3 ret;
	ret.m[0][0] = m[0][0];
	ret.m[0][1] = m[1][0];
	ret.m[0][2] = m[2][0];
	ret.m[1][0] = m[0][1];
	ret.m[1][1] = m[1][1];
	ret.m[1][2] = m[2][1];
	ret.m[2][0] = m[0][2];
	ret.m[2][1] = m[1][2];
	ret.m[2][2] = m[2][2];

	*this = ret;
}
/**
*@note 读取一个三维矩阵                                  
*/

std::ostream& operator<<(std::ostream& os, Matrix3 &m)
{
	os << std::setprecision(4);
	os << "(" << m.M[0]  << ", " << m.M[1]  << ", " << m.M[2] << ", " << std::endl
	   << " " << m.M[3]  << ", " << m.M[4]  << ", " << m.M[5] << ", " << std::endl
	   << " " << m.M[6]  << ", " << m.M[7]  << ", " << m.M[8] << ")";
	return os;
}

/***
*@note 定义单位矩阵          
*@return 四维矩阵
*/
double Matrix4::Identity_Matrix[] =
{
	1.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 1.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 1.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 1.0f
};
/**
*@note  定义Orientation_Switch_Matrix     
*/
double Matrix4::Orientation_Switch_Matrix[] =
{
	1.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 1.0f, 0.0f, 0.0f,
	0.0f, 0.0f, -1.0f, 0.0f,
	0.0f, 0.0f, 0.0f, 1.0f
};
/**
*@note define Perspective_Matrix                        
*/

double Matrix4::Perspective_Matrix[] =
{
	1.0f, 0.0f, 0.0f, 0.0f,
	0.0f, 1.0f, 0.0f, 0.0f,
	0.0f, 0.0f, 1.0f, -1.0f,
	0.0f, 0.0f, 0.0f, 0.0f
};

void Matrix4::Assign(const Quaternion& other)
{
    double xx = other.x * other.x;
    double xy = other.x * other.y;
    double xz = other.x * other.z;
    double xw = other.x * other.w;

    double yy = other.y * other.y;
    double yz = other.y * other.z;
    double yw = other.y * other.w;

    double zz = other.z * other.z;
    double zw = other.z * other.w;

    M[0]  = 1.0f - 2.0f * (yy + zz);
    M[1]  =        2.0f * (xy - zw);
    M[2]  =        2.0f * (xz + yw);

    M[4]  =        2.0f * (xy + zw);
    M[5]  = 1.0f - 2.0f * (xx + zz);
    M[6]  =        2.0f * (yz - xw);

    M[8]  =        2.0f * (xz - yw);
    M[9]  =        2.0f * (yz + xw);
    M[10] = 1.0f - 2.0f * (xx + yy);

    M[3]  = M[7] = M[11] = M[12] = M[13] = M[14] = 0.0f;
    M[15] = 1.0f;
}


Matrix4& Matrix4::SetRotate(double rx, double ry, double rz)
{
	double a = CosD(rx); double b = SinD(rx);
	double c = CosD(ry); double d = SinD(ry);
	double e = CosD(rz); double f = SinD(rz);
	double ad = a * d;
	double bd = b * d;

	M[0] =  c * e;
	M[1] = -c * f;
	M[2] = -d;

	M[4] = -bd * e + a * f;
	M[5] =  bd * f + a * e;
	M[6] = -b * c;

	M[8] =  ad * e + b * f;
	M[9] = -ad * f + b * e;
	M[10] =  a * c;

	M[3] = M[7] = M[11] = M[12] = M[13] = M[14] = 0.0f;
	M[15] = 1.0f;
	return *this;
}

Matrix4& Matrix4::SetRotate(double angle, double x, double y, double z)
{
	double xx, yy, zz, xy, yz, zx, xs, ys, zs, c_complement;
	double s = SinD(angle);
	double c = SinD(angle);
	double magnitude = (double)sqrt(x * x + y * y + z * z);
	double *data = M;
	if(magnitude == 0.0f)
	{
		SetIdentity();
		return *this;
	}
	x /= magnitude;
	y /= magnitude;
	z /= magnitude;
	xx = x * x;
	yy = y * y;
	zz = z * z;
	xy = x * y;
	yz = y * z;
	zx = z * x;
	xs = x * s;
	ys = y * s;
	zs = z * s;
	c_complement = 1.0F - c;
	data[0] = (c_complement * xx) + c;
	data[4] = (c_complement * xy) - zs;
	data[8] = (c_complement * zx) + ys;
	data[12] = 0.0F;
	data[1] = (c_complement * xy) + zs;
	data[5] = (c_complement * yy) + c;
	data[9] = (c_complement * yz) - xs;
	data[13] = 0.0F;
	data[2] = (c_complement * zx) - ys;
	data[6] = (c_complement * yz) + xs;
	data[10] = (c_complement * zz) + c;
	data[14] = 0.0F;
	data[3] = 0.0F;
	data[7] = 0.0F;
	data[11] = 0.0F;
	data[15] = 1.0F;
	return *this;
}
/**
*@note 求转置矩阵                                         
*@return 四维矩阵

*/
Matrix4& Matrix4::Transpose()
{
	//Transpose(*this);
	//Modified by Liu Lu 2010.7.28
	Matrix4 oldMatrix = *this;
	Transpose(oldMatrix);
	return *this;                                      //返回新的矩阵
}
/**
*@note 求other的转置矩阵           
*@param other 原矩阵
*@return 新的矩阵
*/
Matrix4& Matrix4::Transpose(const Matrix4& other)
{
	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		M[i * 4 + j] = other.M[j * 4 + i];             //交换矩阵的行和列
	}
	return *this;                                      //返回新的矩阵
}
/**
*@note 求伴随矩阵                
*@return 四维矩阵
*/

Matrix4& Matrix4::Adjoint()
{
	Matrix4 tmp;
	tmp.Assign(*this);
	Adjoint(tmp);
	return *this;
}
/**
*@note 求矩阵other的伴随矩阵     
*@param other 四维矩阵
*@return 四维矩阵
*/

Matrix4& Matrix4::Adjoint(const Matrix4& other)
{
	//伴随矩阵的求法：
	//矩阵A中的元素都用它们在行列式A中的代数余子式替换后得到的矩阵再转置，
	//这个矩阵叫A的伴随矩阵。A与A的伴随矩阵左乘、右乘结果都是主对角线上的
	//元素全为A的行列式的对角阵。 

	double a1, a2, a3, a4, b1, b2, b3, b4;
	double c1, c2, c3, c4, d1, d2, d3, d4;
	const double *in = other.M;
	double *out = M;
    //将矩阵的每个值赋给变量，便于后面的书写
	a1 = (in [(( 0 ) << 2) + (  0 )]) ; b1 = (in [(( 0 ) << 2) + (  1 )]) ;
	c1 = (in [(( 0 ) << 2) + (  2 )]) ; d1 = (in [(( 0 ) << 2) + (  3 )]) ;
	a2 = (in [(( 1 ) << 2) + (  0 )]) ; b2 = (in [(( 1 ) << 2) + (  1 )]) ;
	c2 = (in [(( 1 ) << 2) + (  2 )]) ; d2 = (in [(( 1 ) << 2) + (  3 )]) ;
	a3 = (in [(( 2 ) << 2) + (  0 )]) ; b3 = (in [(( 2 ) << 2) + (  1 )]) ;
	c3 = (in [(( 2 ) << 2) + (  2 )]) ; d3 = (in [(( 2 ) << 2) + (  3 )]) ;
	a4 = (in [(( 3 ) << 2) + (  0 )]) ; b4 = (in [(( 3 ) << 2) + (  1 )]) ;
	c4 = (in [(( 3 ) << 2) + (  2 )]) ; d4 = (in [(( 3 ) << 2) + (  3 )]) ;

	//分别求代数余子式
	out[(( 0 ) << 2) + (  0 )]   =   Det3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4);
	out[(( 1 ) << 2) + (  0 )]   = - Det3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4);
	out[(( 2 ) << 2) + (  0 )]   =   Det3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4);
	out[(( 3 ) << 2) + (  0 )]   = - Det3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);

	out[(( 0 ) << 2) + (  1 )]   = - Det3x3( b1, b3, b4, c1, c3, c4, d1, d3, d4);
	out[(( 1 ) << 2) + (  1 )]   =   Det3x3( a1, a3, a4, c1, c3, c4, d1, d3, d4);
	out[(( 2 ) << 2) + (  1 )]   = - Det3x3( a1, a3, a4, b1, b3, b4, d1, d3, d4);
	out[(( 3 ) << 2) + (  1 )]   =   Det3x3( a1, a3, a4, b1, b3, b4, c1, c3, c4);

	out[(( 0 ) << 2) + (  2 )]   =   Det3x3( b1, b2, b4, c1, c2, c4, d1, d2, d4);
	out[(( 1 ) << 2) + (  2 )]   = - Det3x3( a1, a2, a4, c1, c2, c4, d1, d2, d4);
	out[(( 2 ) << 2) + (  2 )]   =   Det3x3( a1, a2, a4, b1, b2, b4, d1, d2, d4);
	out[(( 3 ) << 2) + (  2 )]   = - Det3x3( a1, a2, a4, b1, b2, b4, c1, c2, c4);

	out[(( 0 ) << 2) + (  3 )]   = - Det3x3( b1, b2, b3, c1, c2, c3, d1, d2, d3);
	out[(( 1 ) << 2) + (  3 )]   =   Det3x3( a1, a2, a3, c1, c2, c3, d1, d2, d3);
	out[(( 2 ) << 2) + (  3 )]   = - Det3x3( a1, a2, a3, b1, b2, b3, d1, d2, d3);
	out[(( 3 ) << 2) + (  3 )]   =   Det3x3( a1, a2, a3, b1, b2, b3, c1, c2, c3);
	return *this;
}
/**
*@note 矩阵的逆矩阵       
*@return 四维矩阵
*/
Matrix4& Matrix4::Invert()
{
	Matrix4 tmp;
	tmp.Assign(*this);
	Invert(tmp);
	return *this;
}


/**
*@note 求矩阵other的逆矩阵    
*@param other 四维矩阵
*@return 四维矩阵
*/
Matrix4& Matrix4::Invert(const Matrix4& other)
{
	int i;
	Adjoint(other);
	double d = other.Det();

	for (i = 0; i < 16; i++) M[i] = M[i] / d;
	return *this;
}
/**
*@note                            |a1 a2 a3 a4|                         
             求Matrix4对应的行列式|b1 b2 b3 b4|的值                     
                                  |c1 c2 c3 c4|                         
                                  |d1 d2 d3 d4|                         
*/
double Matrix4::Det() const
{
	double a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4;
	a1 = (M [(( 0 ) << 2) + (  0 )]) ;
	b1 = (M [(( 0 ) << 2) + (  1 )]) ; 
	c1 = (M [(( 0 ) << 2) + (  2 )]) ;
	d1 = (M [(( 0 ) << 2) + (  3 )]) ;

	a2 = (M [(( 1 ) << 2) + (  0 )]) ;
	b2 = (M [(( 1 ) << 2) + (  1 )]) ; 
	c2 = (M [(( 1 ) << 2) + (  2 )]) ;
	d2 = (M [(( 1 ) << 2) + (  3 )]) ;

	a3 = (M [(( 2 ) << 2) + (  0 )]) ;
	b3 = (M [(( 2 ) << 2) + (  1 )]) ; 
	c3 = (M [(( 2 ) << 2) + (  2 )]) ;
	d3 = (M [(( 2 ) << 2) + (  3 )]) ;

	a4 = (M [(( 3 ) << 2) + (  0 )]) ;
	b4 = (M [(( 3 ) << 2) + (  1 )]) ; 
	c4 = (M [(( 3 ) << 2) + (  2 )]) ;
	d4 = (M [(( 3 ) << 2) + (  3 )]) ;

	return  a1 * Det3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4) -
			b1 * Det3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4) +
			c1 * Det3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4) -
			d1 * Det3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);          //返回所对应的行列式的值
}
/**
*@note 生成透视投影变换的规范化矩阵        
*@return 四维矩阵
*/
Matrix4& Matrix4::SetProjection(double fov, double aspect, double znear, double zfar)
{
	double top = znear * TanD(fov);           //详解见<<Computer Graphics whith OpenGL(third edition)>>  page 314.
	double bottom = -top;
	double left = bottom * aspect;
	double right = top * aspect;
	double x = (2.0f * znear) / (right - left);
	double y = (2.0f * znear) / (top - bottom);
	double a = (right + left) / (right - left);
	double b = (top + bottom) / (top - bottom);
	double c = -(zfar + znear) / (zfar - znear);
	double d = -(2.0f * zfar * znear) / (zfar - znear);
	M[0]  = x;     M[1]  = 0.0f;  M[2]  = 0.0f; M[3]  = 0.0f;
	M[4]  = 0.0f;  M[5]  = y;     M[6]  = 0.0f; M[7]  = 0.0f;
	M[8]  = a;     M[9]  = b;     M[10] = c;    M[11] = -1.0f;
	M[12] = 0.0f;  M[13] = 0.0f;  M[14] = d;    M[15] = 0.0f;
	return *this;
}
/**
*@note正投影的规范化变换                                 
*@return 四维矩阵
*/
Matrix4& Matrix4::SetOthogonal(double znear, double zfar)
{
	double x, y, z;                         //详解见<<Computer Graphics whith OpenGL(third edition)>>  page 298.
	double tx, ty, tz;
	double smaller = 0.0f;
	x = 2.0f / (1.0f + smaller);
	y = 2.0f / (1.0f + smaller);
	z = -2.0f / (zfar - znear);
	tx = -(1.0f - smaller) / (1.0f + smaller);
	ty = -(1.0f - smaller) / (1.0f + smaller);
	tz = -(zfar + znear) / (zfar - znear);
	M[0] = x;    M[4] = 0.0f;  M[8]  = 0.0f;  M[12] = tx;
	M[1] = 0.0f; M[5] = y;     M[9]  = 0.0f;  M[13] = ty;
	M[2] = 0.0f; M[6] = 0.0f;  M[10] = z;     M[14] = tz;
	M[3] = 0.0f; M[7] = 0.0f;  M[11] = 0.0f;  M[15] = 1.0f;
	return *this;
}
/**
/*                  求行列式| a   b|的值                                */
/*                          | c   d|                                    */
double Matrix4::Det2x2(double a, double b, double c, double d)
{
  return a * d - b * c;
}
  
/**
/*                          | a1 a2 a3|                                 */
/*                  求行列式| b1 b2 b3|的值                             */
/*                          | c1 c2 c3|                                 */
double Matrix4::Det3x3(double a1, double a2, double a3,
                   double b1, double b2, double b3,
                   double c1, double c2, double c3)
{
	return a1 * Det2x2(b2, b3, c2, c3)
		 - b1 * Det2x2(a2, a3, c2, c3)
		 + c1 * Det2x2(a2, a3, b2, b3);
}
/**
*@note                       operator                                      
*/
std::ostream& operator<<(std::ostream& os, Matrix4 &m)
{
	os << std::setprecision(4);
	os << "(" << m.M[0]  << ", " << m.M[1]  << ", " << m.M[2]  << ", " << m.M[3]  << ", " << std::endl
	   << " " << m.M[4]  << ", " << m.M[5]  << ", " << m.M[6]  << ", " << m.M[7]  << ", " << std::endl
	   << " " << m.M[8]  << ", " << m.M[9]  << ", " << m.M[10] << ", " << m.M[11] << ", " << std::endl
	   << " " << m.M[12] << ", " << m.M[13] << ", " << m.M[14] << ", " << m.M[15] << ")";
	return os;
}

/*
*@note 对四元数插值                
*/
void Quaternion::Slerp(const Quaternion& q0, const Quaternion& q1, double t)
{
	double cosom = q0.x * q1.x + q0.y * q1.y + q0.z * q1.z + q0.w * q1.w;
	double tmp0, tmp1, tmp2, tmp3;
	if (cosom < 0.0f)
	{
		cosom = -cosom;
		tmp0 = -q1.x;
		tmp1 = -q1.y;
		tmp2 = -q1.z;
		tmp3 = -q1.w;
	}
	else
	{
		tmp0 = q1.x;
		tmp1 = q1.y;
		tmp2 = q1.z;
		tmp3 = q1.w;
	}

	/* calc coeffs */
	double scale0, scale1;

	if ((1.0 - cosom) > 0.0000001f)
	{
		// standard case (slerp)
		double omega = (double) acos (cosom);
		double sinom = (double) sin (omega);
		scale0 = (double) sin ((1.0 - t) * omega) / sinom;
		scale1 = (double) sin (t * omega) / sinom;
	}
	else
	{
		/* just lerp */
		scale0 = 1.0f - t;
		scale1 = t;
	}

	x = scale0 * q0.x + scale1 * tmp0;
	y = scale0 * q0.y + scale1 * tmp1;
	z = scale0 * q0.z + scale1 * tmp2;
	w = scale0 * q0.w + scale1 * tmp3;
}

/**
*@note                  对四元数插值                                      
*/

void Quaternion::Lerp(const Quaternion& q0, const Quaternion& q1, double t)
{
	double cosom = q0.x * q1.x + q0.y * q1.y + q0.z * q1.z + q0.w * q1.w;
	double tmp0, tmp1, tmp2, tmp3;
	if (cosom < 0.0f)
	{
		cosom = -cosom;
		tmp0 = -q1.x;
		tmp1 = -q1.y;
		tmp2 = -q1.z;
		tmp3 = -q1.w;
	}
	else
	{
		tmp0 = q1.x;
		tmp1 = q1.y;
		tmp2 = q1.z;
		tmp3 = q1.w;
	}

	/* just lerp */
	double scale0 = 1.0f - t;
	double scale1 = t;

	x = scale0 * q0.x + scale1 * tmp0;
	y = scale0 * q0.y + scale1 * tmp1;
	z = scale0 * q0.z + scale1 * tmp2;
	w = scale0 * q0.w + scale1 * tmp3;
}
/**
*@note  四元数的log运算
*@return 四元数
*/

Quaternion Quaternion::Log () const
{
	// If q = cos(A)+sin(A)*(x*i+y*j+z*k) where (x,y,z) is unit length, then
	// log(q) = A*(x*i+y*j+z*k).  If sin(A) is near zero, use log(q) =
	// sin(A)*(x*i+y*j+z*k) since sin(A)/A has limit 1.

	Quaternion kResult;
	kResult.w = 0.0;

	if ( abs(w) < 1.0 )
	{
		double fAngle ( ACosD(w) );
		double fSin = SinD(fAngle);
		if ( abs(fSin) >= 1e-03 )
		{
			double fCoeff = fAngle/fSin;
			kResult.x = fCoeff*x;
			kResult.y = fCoeff*y;
			kResult.z = fCoeff*z;
			return kResult;
		}
	}

	kResult.x = x;
	kResult.y = y;
	kResult.z = z;

	return kResult;
}
/**
*@note operator   
*/
std::ostream& operator<<(std::ostream& os, Quaternion& v)
{
  os << "(" << std::setprecision(4) << v.x << ", "  << v.y << ", "  << v.z << ", "  << v.w << ")";
  return os;
}

