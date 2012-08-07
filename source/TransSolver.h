#ifndef TransSolver_h__
#define TransSolver_h__

#include <QMatrix>
#include "Mathematics.h"
#include <vector>
#include <cv.h>
#include <string>
#include <QImage>

using std::vector;
using std::string;

namespace Caca
{
	struct RigidTransData 
	{
		vector< vector<CvMat*> > m_Avec2;
		vector<Vector2> m_pc;
		vector<double> m_Dis;
	};

	class TransSolver
	{
	public:
		static QMatrix similarityTrans(const vector<Vector2>& src,
			const vector<Vector2>& dst);

		static QMatrix affineTrans(const vector<Vector2>& src, 
			const vector<Vector2>& dst);

		static QMatrix rigidTrans(const vector<Vector2>& src, 
			const vector<Vector2>& dst);

		static QMatrix rigidTrans( const vector<Vector2> &src, 
			const vector<Vector2> &dst, 
			Vector2 v, 
			vector< vector<Vector2> >* AVec1 = 0,
			double* dis = 0,
			Vector2* newV = 0,
			int a = 2);
		static void rigidTransFast(vector<Vector2> &src,
			vector<Vector2> &dst,
			Vector2& v,
			vector< vector<Vector2> >* AVec1 /* = 0 */, 
			double* dis /* = 0 */, 
			Vector2* newV /* = 0 */,
			int a = 2);
		static void movingLeastDeform(const QImage& src, QImage& dst, 
			CvMat* mapx, CvMat* mapy,
			vector<Vector2>& srccp, 
			vector<Vector2>& dstcp,
			vector< vector< vector<Vector2> > >* AVec2 = 0,
			vector<double>* disVec = 0);
		static void movingLeastDeform(const vector< vector<Vector2> >& src, 
			vector< vector<Vector2> >& dst,
			vector<Vector2>& srccp, 
			vector<Vector2>& dstcp,
			vector< vector< vector<Vector2> > >* AVec2 = 0,
			vector<double>* disVec = 0,
			int a = 2);

		static void printCVMAT( CvMat* mat, string name );
	protected:
	private:
	};


}

#endif // TransSolver_h__