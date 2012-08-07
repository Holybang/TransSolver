#include "transsolver.h"
#include "..\Caca\Source\Log\LogWriter.h"
#include <highgui.h>

using namespace Caca;

QMatrix TransSolver::affineTrans( const vector<Vector2>& src, const vector<Vector2>& dst )
{
	if(dst.empty() || src.empty()) {
		//printf("input of affineTrans is empty\n");
		//WriteLog("input of affineTrans is empty\n");
		return QMatrix();
	}
	Vector2 pc, qc;
	for(int i = 0; i < src.size() && i < dst.size(); i++) {
		pc += src[i];
		qc += dst[i];
	}
	pc /= src.size() < dst.size() ? src.size() : dst.size();
	qc /= src.size() < dst.size() ? src.size() : dst.size();

	CvMat *pit, *pi, *qi;
	pit = cvCreateMat(2, 1, CV_64FC1);
	pi = cvCreateMat(1, 2, CV_64FC1);
	qi = cvCreateMat(1, 2, CV_64FC1);
	CvMat *spitpi = cvCreateMat(2, 2, CV_64FC1);
	CvMat *pitpi = cvCreateMat(2, 2, CV_64FC1);
	CvMat *pitqi = cvCreateMat(2, 2, CV_64FC1);
	CvMat *spitqi = cvCreateMat(2, 2, CV_64FC1);
	cvZero(spitpi);
	cvZero(spitqi);

	for(int i = 0; i < src.size() && i < dst.size(); i++) {		
		Vector2 qpi = src[i] - pc;
		Vector2 qqi = dst[i] - qc;
		cvmSet(pit, 0, 0, qpi.x);
		cvmSet(pit, 1, 0, qpi.y);
		cvmSet(pi, 0, 0, qpi.x);
		cvmSet(pi, 0, 1, qpi.y);
		cvmSet(qi, 0, 0, qqi.x);
		cvmSet(qi, 0, 1, qqi.y);
		cvMatMul(pit, pi, pitpi);
		cvAdd(pitpi, spitpi, spitpi);
		cvMatMul(pit, qi, pitqi);
		cvAdd(pitqi, spitqi, spitqi);
	}

	//printCVMAT(spitpi, "spitpi");
	//printCVMAT(spitqi, "spitqi");

	CvMat* ispitpi = cvCreateMat(2, 2, CV_64FC1);
	cvInvert(spitpi, ispitpi);
	//printCVMAT(ispitpi, "ispitpi");

	CvMat* M = cvCreateMat(2, 2, CV_64FC1);
	cvMatMul(ispitpi, spitqi, M);
	//printCVMAT(M, "M");

	double m11 = cvmGet(M, 0, 0);
	double m12 = cvmGet(M, 0, 1);
	double m21 = cvmGet(M, 1, 0);
	double m22 = cvmGet(M, 1, 1);
	QMatrix qm(m11, m12, m21, m22, 0, 0);
	QMatrix pcm(1.0, 0, 0, 1.0, -pc.x, -pc.y);
	QMatrix qcm(1.0, 0, 0, 1.0, qc.x, qc.y);
	QMatrix ret = pcm*qm*qcm;

	cvReleaseMat(&pit);
	cvReleaseMat(&pi);
	cvReleaseMat(&qi);
	cvReleaseMat(&spitpi);
	cvReleaseMat(&pitpi);
	cvReleaseMat(&pitqi);
	cvReleaseMat(&spitqi);
	cvReleaseMat(&ispitpi);
	cvReleaseMat(&M);

	return ret;
}

QMatrix Caca::TransSolver::similarityTrans( const vector<Vector2>& src, const vector<Vector2>& dst )
{
	if(dst.empty() || src.empty()) {
		//printf("input of affineTrans is empty\n");
		//WriteLog("input of affineTrans is empty\n");
		return QMatrix();
	}
	Vector2 pc, qc;
	for(int i = 0; i < src.size() && i < dst.size(); i++) {
		pc += src[i];
		qc += dst[i];
	}
	pc /= src.size() < dst.size() ? src.size() : dst.size();
	qc /= src.size() < dst.size() ? src.size() : dst.size();

	CvMat *m0 = cvCreateMat(2, 2, CV_64FC1);
	CvMat *m1 = cvCreateMat(2, 2, CV_64FC1);
	CvMat *m = cvCreateMat(2, 2, CV_64FC1);
	CvMat *M = cvCreateMat(2, 2, CV_64FC1);

	cvSetZero(m0);
	cvSetZero(m1);
	cvSetZero(m);
	cvSetZero(M);

	double us = 0;
	for(int i = 0; i , src.size() && i < dst.size(); i++) {
		Vector2 pi = src[i] - pc;
		Vector2 qi = dst[i] - qc;

		us += pi.Dot(pi);

		Vector2 pip(-pi.y, pi.x);
		Vector2 qip(-qi.y, qi.x);

		cvmSet(m0, 0, 0, pi.x);
		cvmSet(m0, 0, 1, pi.y);
		cvmSet(m0, 1, 0, -pip.x);
		cvmSet(m0, 1, 1, -pip.y);

		cvmSet(m1, 0, 0, qi.x);
		cvmSet(m1, 1, 0, qi.y);
		cvmSet(m1, 0, 1, -qip.x);
		cvmSet(m1, 1, 1, -qip.y);

		cvMatMul(m0, m1, m);

		cvAdd(M, m, M);

// 		printCVMAT(m0, "m0");
// 		printCVMAT(m1, "m1");
// 		printCVMAT(m, "m");
// 		printCVMAT(M, "M");
	}

	cvScale(M, M, 1.0/us);

	double m11 = cvmGet(M, 0, 0);
	double m12 = cvmGet(M, 0, 1);
	double m21 = cvmGet(M, 1, 0);
	double m22 = cvmGet(M, 1, 1);
	QMatrix qm(m11, m12, m21, m22, 0, 0);
	QMatrix pcm(1.0, 0, 0, 1.0, -pc.x, -pc.y);
	QMatrix qcm(1.0, 0, 0, 1.0, qc.x, qc.y);
	QMatrix ret = pcm*qm*qcm;

	cvReleaseMat(&m0);
	cvReleaseMat(&m1);
	cvReleaseMat(&m);
	cvReleaseMat(&M);

	return ret;
}

QMatrix TransSolver::rigidTrans( const vector<Vector2> &src, const vector<Vector2> &dst, 
						 Vector2 v,
						 vector< vector<Vector2> >* AVec1 /*= 0*/,
						 double* dis /*= 0*/,
						 Vector2* newV /*= 0*/,
						 int a /*= 2*/)
{
	if(dst.empty() || src.empty()) {
		return QMatrix();
	}
	

	vector<double> weights(src.size());
	Vector2 pc, qc;
	double wsum = 0;

	for(int i = 0; i < src.size(); i++) {
		double dis = src[i].Distance(v);
		if(dis < M_LOW_TOLERANCE) {
			dis = M_LOW_TOLERANCE;
		}
		double w = 1.0/pow(dis, a);
		weights[i] = w;
		pc += src[i]*w;
		qc += dst[i]*w;
		wsum += w;
	}
	pc /= wsum;
	qc /= wsum;

	if(AVec1) {
		Vector2 frv;
		if(AVec1->empty()) {
			AVec1->resize(src.size(), vector<Vector2>(2));
			for(int i = 0; i < src.size(); i++) {
				vector<Vector2>& Ai = AVec1->at(i);

				Vector2 pi = src[i] - pc;
				Vector2 pi_(pi.y, -pi.x);				
				Vector2 vpc = v - pc;
				Vector2 vpc_(vpc.y, -vpc.x);
				double& wi = weights[i];

				Ai[0].x = pi.Dot(vpc)*wi;
				Ai[0].y = pi_.Dot(vpc)*wi;
				Ai[1].x = pi.Dot(vpc_)*wi;
				Ai[1].y = pi_.Dot(vpc_)*wi; 

				Vector2 qi = dst[i] - qc;
				Vector2 qiAi(qi.Dot(Ai[0]), qi.Dot(Ai[1]));
				frv += qiAi;
			}
		} else {
			for(int i = 0; i < src.size(); i++) {
				vector<Vector2>& Ai = AVec1->at(i);
				Vector2 qi = dst[i] - qc;
				Vector2 qiAi(qi.Dot(Ai[0]), qi.Dot(Ai[1]));
				frv += qiAi;
			}
		}

		frv.Normalize();

		if(dis) {
			if(*dis < 0) {
				*dis = v.Distance(pc);
			}
		}
		*newV = frv*(*dis) + qc;
		
		return QMatrix();
	} else {
		double u = 0;
		double u1, u2;
		u1 = 0;
		u2 = 0;
		for(int i = 0; i < src.size() && i < dst.size(); i++) {		
			Vector2 qpi = src[i] - pc;
			Vector2 qqi = dst[i] - qc;
			Vector2 pi(qpi.x, qpi.y);
			Vector2 qi(qqi.x, qqi.y);
			u1 += pi.Dot(qi)*weights[i];
			Vector2 pi_(pi.y, -pi.x);
			u2 += qi.Dot(pi_)*weights[i];
		}
		u = sqrt(u1*u1 + u2*u2);
		if(u < M_LOW_TOLERANCE) {
			u = M_LOW_TOLERANCE;
		}

		CvMat* R = cvCreateMat(2, 2, CV_64FC1);
		CvMat* r = cvCreateMat(2, 2, CV_64FC1);
		for(int i = 0; i < 4; i++) {
			cvSetReal1D(R, i, 0);
			cvSetReal1D(r, i, 0);
		}
		for(int i = 0; i < src.size() && i < dst.size(); i++) {
			Vector2 qpi = src[i] - pc;
			Vector2 qqi = dst[i] - qc;
			Vector2 pi(qpi.x, qpi.y);
			Vector2 qi(qqi.x, qqi.y);
			Vector2 pi_(pi.y, -pi.x);
			Vector2 qi_(qi.y, -qi.x);		
			cvmSet(r, 0, 0, pi.Dot(qi));
			cvmSet(r, 0, 1, pi.Dot(qi_));
			cvmSet(r, 1, 0, pi_.Dot(qi));
			cvmSet(r, 1, 1, pi_.Dot(qi_));		
			CvScalar;
			cvAXPY(r, weights[i]/u, R, R);
			//cvAdd(R, r, R);
		}

		double m11 = cvmGet(R, 0, 0);
		double m12 = cvmGet(R, 0, 1);
		double m21 = cvmGet(R, 1, 0);
		double m22 = cvmGet(R, 1, 1);
		QMatrix qm(m11, m12, m21, m22, 0, 0);
		QMatrix pcm(1.0, 0, 0, 1.0, -pc.x, -pc.y);
		QMatrix qcm(1.0, 0, 0, 1.0, qc.x, qc.y);
		QMatrix ret = pcm*qm*qcm;

		cvReleaseMat(&r);
		cvReleaseMat(&R);

		return ret;
	}
}


void TransSolver::printCVMAT( CvMat* mat, string name )
{
	printf("%s\n", name.data());
	for(int i = 0; i < mat->rows; i++) {
		for(int j = 0; j < mat->cols; j++) {
			printf("%lf	", cvmGet(mat, i, j));
		}
		printf("\n");
	}
}

void TransSolver::movingLeastDeform( const QImage& src, QImage& dst,
							 CvMat* mapx, CvMat* mapy,
							 vector<Vector2>& srccp, 
							 vector<Vector2>& dstcp,
							vector< vector< vector<Vector2> > > * AVec2 /*= 0*/,
							vector<double>* disVec /*= 0*/)
{
	//static int gen_map_time = 0;
	//static int remap_time = 0;

	//QTime time;
	//time.start();

	if(srccp.size() != dstcp.size()) {
		printf("the numbers (%d %d) of src control points and dst control points must be equal.\n", srccp.size(), dstcp.size());
		WriteLog(QString("the numbers (%d %d) of src control points and dst control points must be equal.").arg(srccp.size()).arg(dstcp.size()));
		return;
	}
	int w = src.width(); 
	int h = src.height();
	if(AVec2) {
		if(AVec2->empty()) {
			AVec2->resize(w*h);
			disVec->resize(w*h, -1.0);
			for(int y = 0; y < h; y++) {
				for(int x = 0; x < w; x++) {
					Vector2 v(x, y);	
					QMatrix trans;
					int idx = x + y*w;
					vector< vector<Vector2> > &AVec1 = AVec2->at(idx);
					double& dis = disVec->at(idx);
					Vector2 newV;
					rigidTrans(dstcp, srccp, v, &AVec1, &dis, &newV);					
					cvmSet(mapx, y, x, newV.x);
					cvmSet(mapy, y, x, newV.y);
				}
			}
		} else {
			for(int y = 0; y < h; y++) {
				for(int x = 0; x < w; x++) {
					Vector2 v(x, y);
					int idx = x + y*w;
					vector< vector<Vector2> > &AVec1 = (*AVec2)[idx];
					double& dis = (*disVec)[idx];
					Vector2 newV;
					rigidTransFast(dstcp, srccp, v, &AVec1, &dis, &newV);					
					cvmSet(mapx, y, x, newV.x);
					cvmSet(mapy, y, x, newV.y);
				}
			}
		}

		
	} else {
		for(int y = 0; y < h; y++) {
			for(int x = 0; x < w; x++) {
				Vector2 v(x, y);	
				QMatrix trans;
				trans = rigidTrans(dstcp, srccp, v);
				QPointF qv(x, y);
				qv = trans.map(qv);
				cvmSet(mapx, y, x, qv.x());
				cvmSet(mapy, y, x, qv.y());	
			}
		}
	}

	//gen_map_time += time.elapsed();
	//time.restart();
	
	//cvSaveImage("mapx.jpg", mapx);
	//cvSaveImage("mapy.jpg", mapy);
	
	CvSize size = cvSize(src.width(), src.height());
	IplImage* srcIplImg = cvCreateImageHeader(size, IPL_DEPTH_8U, 4);
	srcIplImg->imageData = (char*)(src.scanLine(0));
	IplImage* tarIplImg = cvCreateImageHeader(size, IPL_DEPTH_8U, 4);
	tarIplImg->imageData = (char*)(dst.scanLine(0));
	cvRemap(srcIplImg, tarIplImg, mapx, mapy, 
		CV_INTER_LINEAR+CV_WARP_FILL_OUTLIERS, cvScalar(0, 0, 0, 255));
	cvReleaseImageHeader(&tarIplImg);
	cvReleaseImageHeader(&srcIplImg);

	//remap_time += time.elapsed();

	//printf("	gen map time: %d\n", gen_map_time);
	//printf("	remap time: %d\n", remap_time);
}

void TransSolver::movingLeastDeform( const vector< vector<Vector2> >& src, 
							 vector< vector<Vector2> >& dst, 
							 vector<Vector2>& srccp, 
							 vector<Vector2>& dstcp, 
							 vector< vector< vector<Vector2> > >* AVec2 /*= 0*/, 
							 vector<double>* disVec /*= 0*/,
							 int a /*= 2*/)
{
	if(srccp.size() != dstcp.size()) {
		printf("the numbers (%d %d) of src control points and dst control points must be equal.\n", srccp.size(), dstcp.size());
		WriteLog(QString("the numbers (%d %d) of src control points and dst control points must be equal.").arg(srccp.size()).arg(dstcp.size()));
		return;
	}
	int w = src.front().size(); 
	int h = src.size();
	vector<Vector2> vec = src.front();
	dst.resize(src.size(), vec);
	if(AVec2) {
		if(AVec2->empty()) {
			AVec2->resize(w*h);
			disVec->resize(w*h, -1.0);
			for(int y = 0; y < h; y++) {
				for(int x = 0; x < w; x++) {
					Vector2 v = src[y][x];	
					QMatrix trans;
					int idx = x + y*w;
					vector< vector<Vector2> > &AVec1 = AVec2->at(idx);
					double& dis = disVec->at(idx);
					Vector2 newV;
					rigidTrans(srccp, dstcp, v, &AVec1, &dis, &newV, a);	
					dst[y][x] = newV;
				}
			}
		} else {
			for(int y = 0; y < h; y++) {
				for(int x = 0; x < w; x++) {
					Vector2 v = src[y][x];
					int idx = x + y*w;
					vector< vector<Vector2> > &AVec1 = (*AVec2)[idx];
					double& dis = (*disVec)[idx];
					Vector2 newV;
					rigidTransFast(srccp, dstcp, v, &AVec1, &dis, &newV);	
					dst[y][x] = newV;					
				}
			}
		}
	} else {
		for(int y = 0; y < h; y++) {
			for(int x = 0; x < w; x++) {
				Vector2 v = src[y][x];	
				QMatrix trans;
				trans = rigidTrans(srccp, dstcp, v, 0, 0, 0, a);
				QPointF qv(v.x, v.y);
				qv = trans.map(qv);
				dst[y][x] = Vector2(qv.x(), qv.y());
			}
		}
	}
}


void TransSolver::rigidTransFast( vector<Vector2> &src, 
						  vector<Vector2> &dst, 
						  Vector2& v,
						  vector< vector<Vector2> >* AVec1 /* = 0 */, 
						  double* dis /* = 0 */,
						  Vector2* newV /* = 0 */,
						  int a /*= 2*/)
{
	int count = src.size();
	static vector<double> weights(count);
	Vector2 pc;
	Vector2 qc;
	double wsum = 0;
	//static int a = 2;
	//qc.x = 0;
	//qc.y = 0;
	//pc.x = 0;
	//pc.y = 0;
	//wsum = 0;

	double w;
	double d;
	vector<Vector2>::iterator srcIter = src.begin();
	vector<Vector2>::iterator dstIter = dst.begin();
	vector<double>::iterator wIter = weights.begin();
	while(srcIter != src.end()) {
		Vector2& p = *srcIter;
		d = p.Distance(v);
		if(d < M_LOW_TOLERANCE) {
			d = M_LOW_TOLERANCE;
		}
		w = 1.0/(pow(d, a));
		(*wIter) = w;
		pc += p*w;
		qc += (*dstIter)*w;
		wsum += w;
		srcIter++;
		dstIter++;
		wIter++;
	}
	//for(int i = 0; i < count; i++) {
	//	Vector2& p = src[i];
	//	d = p.Distance(v);
	//	if(d < M_LOW_TOLERANCE) {
	//		d = M_LOW_TOLERANCE;
	//	}
	//	w = 1.0/(d*d);
	//	weights[i] = w;
	//	pc += p*w;
	//	qc += dst[i]*w;
	//	wsum += w;
	//}
	pc /= wsum;
	qc /= wsum;

	Vector2 frv;
	//frv.x = 0;
	//frv.y = 0;

	vector< vector<Vector2> >::iterator AIter = AVec1->begin();
	Vector2 qi, qiAi;
	dstIter = dst.begin();
	while(AIter != AVec1->end()) {
		vector<Vector2>& Ai = *AIter;
		qi = *dstIter - qc;
		qiAi.x = qi.Dot(Ai[0]);
		qiAi.y = qi.Dot(Ai[1]);
		frv.x += qiAi.x;
		frv.y += qiAi.y;
		AIter++;
		dstIter++;
	}
	//for(int i = 0; i < count; i++) {
	//	vector<Vector2>& Ai = (*AVec1)[i];
	//	qi = dst[i] - qc;
	//	qiAi.x = qi.Dot(Ai[0]);
	//	qiAi.y = qi.Dot(Ai[1]);
	//	frv.x += qiAi.x;
	//	frv.y += qiAi.y;
	//}

	frv.Normalize();

	*newV = frv*(*dis) + qc;


}

QMatrix Caca::TransSolver::rigidTrans( const vector<Vector2>& src, const vector<Vector2>& dst )
{
	if(dst.empty() || src.empty()) {
		return QMatrix();
	}

	vector<double> weights(src.size(), 1.0/src.size());
	Vector2 pc, qc;
	double wsum = 0;

	for(int i = 0; i < src.size(); i++) {		
		double w = 1.0/src.size();
		weights[i] = w;
		pc += src[i]*w;
		qc += dst[i]*w;
		wsum += w;
	}
	pc /= wsum;
	qc /= wsum;

	double u = 0;
	double u1, u2;
	u1 = 0;
	u2 = 0;
	for(int i = 0; i < src.size() && i < dst.size(); i++) {		
		Vector2 qpi = src[i] - pc;
		Vector2 qqi = dst[i] - qc;
		Vector2 pi(qpi.x, qpi.y);
		Vector2 qi(qqi.x, qqi.y);
		u1 += pi.Dot(qi)*weights[i];
		Vector2 pi_(pi.y, -pi.x);
		u2 += qi.Dot(pi_)*weights[i];
	}
	u = sqrt(u1*u1 + u2*u2);
	if(u < M_LOW_TOLERANCE) {
		u = M_LOW_TOLERANCE;
	}

	CvMat* R = cvCreateMat(2, 2, CV_64FC1);
	CvMat* r = cvCreateMat(2, 2, CV_64FC1);
	for(int i = 0; i < 4; i++) {
		cvSetReal1D(R, i, 0);
		cvSetReal1D(r, i, 0);
	}
	for(int i = 0; i < src.size() && i < dst.size(); i++) {
		Vector2 qpi = src[i] - pc;
		Vector2 qqi = dst[i] - qc;
		Vector2 pi(qpi.x, qpi.y);
		Vector2 qi(qqi.x, qqi.y);
		Vector2 pi_(pi.y, -pi.x);
		Vector2 qi_(qi.y, -qi.x);		
		cvmSet(r, 0, 0, pi.Dot(qi));
		cvmSet(r, 0, 1, pi.Dot(qi_));
		cvmSet(r, 1, 0, pi_.Dot(qi));
		cvmSet(r, 1, 1, pi_.Dot(qi_));		
		CvScalar;
		cvAXPY(r, weights[i]/u, R, R);
		//cvAdd(R, r, R);
	}

	double m11 = cvmGet(R, 0, 0);
	double m12 = cvmGet(R, 0, 1);
	double m21 = cvmGet(R, 1, 0);
	double m22 = cvmGet(R, 1, 1);
	QMatrix qm(m11, m12, m21, m22, 0, 0);
	QMatrix pcm(1.0, 0, 0, 1.0, -pc.x, -pc.y);
	QMatrix qcm(1.0, 0, 0, 1.0, qc.x, qc.y);
	QMatrix ret = pcm*qm*qcm;

	cvReleaseMat(&r);
	cvReleaseMat(&R);

	return ret;
}
