// g2o - General Graph Optimization
// Copyright (C) 2011 H. Strasdat
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Modified by Raúl Mur Artal (2014)
// Added EdgeSE3ProjectXYZ (project using focal_length in x,y directions)
// Modified by Raúl Mur Artal (2016)
// Added EdgeStereoSE3ProjectXYZ (project using focal_length in x,y directions)
// Added EdgeSE3ProjectXYZOnlyPose (unary edge to optimize only the camera pose)
// Added EdgeStereoSE3ProjectXYZOnlyPose (unary edge to optimize only the camera pose)

#ifndef G2O_SIX_DOF_TYPES_EXPMAP
#define G2O_SIX_DOF_TYPES_EXPMAP

#include "../core/base_vertex.h"
#include "../core/base_binary_edge.h"
#include "../core/base_unary_edge.h"
#include "se3_ops.h"
#include "se3quat.h"
#include "types_sba.h"
#include "so3.h"
#include <Eigen/Geometry>
#include <Eigen/StdVector>


namespace g2o {
namespace types_six_dof_expmap {
void init();
}

using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6d;


/**
 * \brief SE3 Vertex parameterized internally with a transformation matrix
 and externally with its exponential map
 */


class  VertexSE3Expmap : public BaseVertex<6, SE3Quat>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:

  VertexSE3Expmap();
  ~VertexSE3Expmap();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  virtual void setToOriginImpl() {
    _estimate = SE3Quat();
  }

  virtual void oplusImpl(const double* update_)  {
    Eigen::Map<const Vector6d> update(update_);
    setEstimate(SE3Quat::exp(update)*estimate());
  }
};


/**
 * \brief SO3 Vertex parameterized internally with a rotation matrix and externally with its exponential map
 */

class VertexSO3Expmap : public BaseVertex<3, SO3>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:

  VertexSO3Expmap();
  bool read(std::istream& is);
  bool write(std::ostream& os) const;

  virtual void setToOriginImpl(){
    _estimate = SO3();
  }

  virtual void oplusImpl(const double* update_){
    Eigen::Map<const Vector3d> update(update_);
    SO3 tmp(SO3::exp(update));
    setEstimate(tmp*estimate());      // exp(w + delta_w ) = exp(delta_w ) * exp(w); left multipline
  }
};

class  EdgeSE3ProjectXYZ: public  BaseBinaryEdge<2, Vector2d, VertexSBAPointXYZ, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:

  EdgeSE3ProjectXYZ();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    Vector2d obs(_measurement);
    _error = obs-cam_project(v1->estimate().map(v2->estimate()));
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    return (v1->estimate().map(v2->estimate()))(2)>0.0;
  }
    

  virtual void linearizeOplus();

  Vector2d cam_project(const Vector3d & trans_xyz) const;

  double fx, fy, cx, cy;
};


class  EdgeStereoSE3ProjectXYZ: public  BaseBinaryEdge<3, Vector3d, VertexSBAPointXYZ, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:

  EdgeStereoSE3ProjectXYZ();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    Vector3d obs(_measurement);
    _error = obs - cam_project(v1->estimate().map(v2->estimate()),bf);
    //std::cout << "EdgeStereoSE3ProjectXYZ _error : " << _error.transpose() << std::endl;
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    return (v1->estimate().map(v2->estimate()))(2)>0.0;
  }


  virtual void linearizeOplus();

  Vector3d cam_project(const Vector3d & trans_xyz, const float &bf) const;

  double fx, fy, cx, cy, bf;
};

class  EdgeSE3ProjectXYZOnlyPose: public  BaseUnaryEdge<2, Vector2d, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:

  EdgeSE3ProjectXYZOnlyPose();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    Vector2d obs(_measurement);
    _error = obs-cam_project(v1->estimate().map(Xw));
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    return (v1->estimate().map(Xw))(2)>0.0;
  }


  virtual void linearizeOplus();

  Vector2d cam_project(const Vector3d & trans_xyz) const;

  Vector3d Xw;
  double fx, fy, cx, cy;
};


class  EdgeStereoSE3ProjectXYZOnlyPose: public  BaseUnaryEdge<3, Vector3d, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:

  EdgeStereoSE3ProjectXYZOnlyPose();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    Vector3d obs(_measurement);
    _error = obs - cam_project(v1->estimate().map(Xw));
    Vector7d est = v1->estimate().toVector();
    //std::cout << "_measurement: " << _measurement.transpose() << ", Xw: " << Xw.transpose()<<
    //", estimate: " << est.transpose() << std::endl;
    //std::cout << "_error = " << _error.transpose() << std::endl;
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    return (v1->estimate().map(Xw))(2)>0.0;
  }


  virtual void linearizeOplus();

  Vector3d cam_project(const Vector3d & trans_xyz) const;

  Vector3d Xw;
  double fx, fy, cx, cy, bf;
};

class EdgeSE3ProjectLineXYZ : public BaseBinaryEdge<2, Vector2d , VertexSBALineXYZ , VertexSE3Expmap >{

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:

    EdgeSE3ProjectLineXYZ();

    bool read(std::istream& is);
    bool write(std::ostream& os) const;

    void computeError(){
      VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
      SE3Quat T(vj->estimate());
      VertexSBALineXYZ *vi = static_cast<VertexSBALineXYZ *>(_vertices[0]);
      Vector6d line_p_q = vi->estimate();
      Vector3d line_p;
      line_p << line_p_q[0] , line_p_q[1] , line_p_q[2];
      Vector3d line_q;
      line_q << line_p_q[3] , line_p_q[4] , line_p_q[5]; 

      //////////////////////////////////////////////////////////////////////////////////////////
      Vector3d h_pu; 
      Vector3d h_qu; 
      h_pu << pu[0] , pu[1] , 1.0;
      h_qu << qu[0] , qu[1] , 1.0;
      //std::cout << "h_pu: " << h_pu.transpose() << ", h_qu: " << h_qu.transpose() << std::endl;
      Vector3d l_pu_qu = h_pu.cross(h_qu);
      double l_pu_qu_n = sqrt(l_pu_qu[0]*l_pu_qu[0]+l_pu_qu[1]*l_pu_qu[1]);
      l_pu_qu = l_pu_qu / l_pu_qu_n;
      //std::cout << "l_pu_qu = " << l_pu_qu.transpose() << std::endl;

      Vector3d line_p_trans = T.map(line_p);
      Vector3d line_q_trans = T.map(line_q);
      if(fabs(line_q_trans[2])<0.01){
        std::cout << "T : " << T.to_homogeneous_matrix() << ", line_q : " << line_q.transpose() << std::endl;
      }

      //std::cout << "line_p_trans : " << line_p_trans.transpose() << ", line_q_trans: " << line_q_trans.transpose() 
      //<<", line_p_q : " << line_p_q.transpose()<< std::endl;
      Vector3d h_pu_hat, h_qu_hat;
      h_pu_hat << fx * line_p_trans[0] / line_p_trans[2] + cx, fx * line_p_trans[1] / line_p_trans[2] + cy , 1.0;
      h_qu_hat << fx * line_q_trans[0] / line_q_trans[2] + cx, fy * line_q_trans[1] / line_q_trans[2] + cy , 1.0;  

      Vector2d obs;
      obs << l_pu_qu.transpose() * h_pu_hat, l_pu_qu.transpose() * h_qu_hat;
      //std::cout << "obs : " << obs.transpose() << std::endl;
      _error = obs;
      std::cout << "EdgeSE3ProjectLineXYZ error : " << _error.transpose() << std::endl; 
    }

    virtual void linearizeOplus();

    Vector3d cam_project(const Vector3d& trans_xyz) const;
    void setParams(const Vector2d& _pu , const Vector2d& _qu);

    double fx , fy , cx , cy;
    Vector2d pu , qu;          // the line's two endpoint measurement in the image

};

class EdgeSE3ProjectLineXYZOnlyPose : public BaseUnaryEdge<2, Vector2d , VertexSE3Expmap >{

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:

    EdgeSE3ProjectLineXYZOnlyPose();

    bool read(std::istream& is);
    bool write(std::ostream& os) const;

    void computeError(){
      VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[0]);
      SE3Quat T(vj->estimate());
      //VertexSBALineXYZ *vi = static_cast<VertexSBALineXYZ *>(_vertices[0]);
      //Vector6d line_p_q = vj->estimate();
      Vector3d line_p;//
      line_p << line_p_q[0] , line_p_q[1] , line_p_q[2];
      Vector3d line_q;
      line_q << line_p_q[3] , line_p_q[4] , line_p_q[5]; 

      //////////////////////////////////////////////////////////////////////////////////////////
      Vector3d h_pu;
      h_pu << (pu[0] , pu[1] , 1.0);
      Vector3d h_qu;
      h_qu << (qu[0] , qu[1] , 1.0);
      Vector3d l_pu_qu = h_pu.cross(h_qu);
      double l_pu_qu_n = sqrt(l_pu_qu[0]*l_pu_qu[0]+l_pu_qu[1]*l_pu_qu[1]);
      l_pu_qu = l_pu_qu / l_pu_qu_n;//l_pu_qu.norm();
      //std::cout << "l_pu_qu = " << l_pu_qu.transpose() << std::endl;

      Vector3d line_p_trans = T.map(line_p);
      Vector3d line_q_trans = T.map(line_q);
      //std::cout << "line_p_trans = " << line_p_trans.transpose() << ", line_q_trans = " << line_q_trans << std::endl;

      Vector3d h_pu_hat, h_qu_hat;
      h_pu_hat << fx * line_p_trans[0] / line_p_trans[2] + cx, fx * line_p_trans[1] / line_p_trans[2] + cy , 1.0;
      h_qu_hat << fx * line_q_trans[0] / line_q_trans[2] + cx, fy * line_q_trans[1] / line_q_trans[2] + cy , 1.0;  

      Vector2d obs;
      obs << l_pu_qu.transpose() * h_pu_hat, l_pu_qu.transpose() * h_qu_hat;
      _error = obs;
      //std::cout << "error = " << obs.transpose() << std::endl;
    }

    virtual void linearizeOplus();

    Vector3d cam_project(const Vector3d& trans_xyz) const;
    void setParams(const Vector2d& _pu , const Vector2d& _qu);
    void setParams(const Vector4f& puqu);
    void setMeasurement3DLine(const Vector6d& _line_p_q);

    double fx , fy , cx , cy;
    Vector2d pu , qu;          // the line's two endpoint measurement in the image
    Vector6d line_p_q;

};



class EdgeSE3ProjectPlaneSO3 : public BaseBinaryEdge<3, Vector3d , VertexSO3Expmap , VertexSE3Expmap> {

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:

    EdgeSE3ProjectPlaneSO3();

    bool read(std::istream& is);
    bool write(std::ostream& os) const;

    void computeError(){

      VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
      SE3Quat T(vj->estimate());
      VertexSO3Expmap *vi = static_cast<VertexSO3Expmap *>(_vertices[0]);
      SO3 plane(vi->estimate());

      Quaterniond pi_x = plane_c.unit_quaternion();      // the observation
      Vector4d pi_x_4d;
      pi_x_4d << pi_x.x() , pi_x.y() , pi_x.z() , pi_x.w();
      
      Vector4d y = T.to_homogeneous_matrix().transpose() * pi_x_4d;
      Vector4d q_y = y / y.norm();
      Quaterniond n_q_y(q_y[3] , q_y[0] , q_y[1] , q_y[2]);
      SO3 q_T(n_q_y);
      SO3 q_u = q_T * plane.inverse();
      Vector3d u = q_u.log();

      Vector3d obs;
      obs = u;
      _error = obs;

    }

    virtual void linearizeOplus();

    void setParams(const SO3& plane_);

    SO3 plane_c;     // the plane coordinates in the current camera frame C     
};

class EdgeSE3ProjectPlaneSO3OnlyPose : public BaseUnaryEdge<3, Vector3d , VertexSE3Expmap> {

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

public:

    EdgeSE3ProjectPlaneSO3OnlyPose();

    bool read(std::istream& is);
    bool write(std::ostream& os) const;

    void computeError(){

      VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[0]);
      SE3Quat T(vj->estimate());
      //VertexSO3Expmap *vi = static_cast<VertexSO3Expmap *>(_vertices[0]);
      SO3 plane = plane_g; //(vi->estimate());

      Quaterniond pi_x = plane_c.unit_quaternion();      // the observation
      Vector4d pi_x_4d;
      pi_x_4d << pi_x.x() , pi_x.y() , pi_x.z() , pi_x.w();
      
      Vector4d y = T.to_homogeneous_matrix().transpose() * pi_x_4d;
      Vector4d q_y = y / y.norm();
      Quaterniond n_q_y(q_y[3] , q_y[0] , q_y[1] , q_y[2]);
      SO3 q_T(n_q_y);
      SO3 q_u = q_T * plane.inverse();
      Vector3d u = q_u.log();

      Vector3d obs;
      obs = u;
      _error = obs;

    }

    virtual void linearizeOplus();

    void setParams(const SO3& plane_, const SO3& _plane_g);

    SO3 plane_c;     // the plane coordinates in the current camera frame C    
    SO3 plane_g; 
};

class EdgeParallel: public BaseBinaryEdge<3, Vector3d , VertexSBAPluckerLine, VertexSBAPluckerLine> {

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
public:
    EdgeParallel() {};

    void computeError() {
        VertexSBAPluckerLine *vi0 = static_cast<VertexSBAPluckerLine *>(_vertices[0]), *vi1 = static_cast<VertexSBAPluckerLine *>(_vertices[1]);
        PluckerLine pl_w0 = vi0->estimate(), pl_w1 = vi1->estimate();
        Vector3d n0 = pl_w0.normal(), n1 = pl_w1.normal(), v0 = pl_w0.direction(), v1 = pl_w1.direction();
        v1 /= v1.norm();
        v0 /= v0.norm();
        _error = v0.cross(v1);
    }

    virtual void linearizeOplus();
    bool read(std::istream& is) {};
    bool write(std::ostream& os) const {};
    
};



class EdgeSE3ProjectPluckerLine : public BaseBinaryEdge<2 , Vector2d , VertexSBAPluckerLine , VertexSE3Expmap>{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  public:
    EdgeSE3ProjectPluckerLine();

    bool read(std::istream& is);
    bool write(std::ostream& os) const;

    void computeError(){

      VertexSE3Expmap* vj = static_cast<VertexSE3Expmap*>(_vertices[1]);
      SE3Quat T(vj->estimate());
      VertexSBAPluckerLine* vi = static_cast<VertexSBAPluckerLine *>(_vertices[0]);
      PluckerLine pl_w = vi->estimate();

      PluckerLine pl_c = pl_w.transform(T.to_homogeneous_matrix());
      //std::cout << "pl_w : " << pl_w << ", pl_c : " << pl_c << std::endl;

      Eigen::Matrix3d K;
      K << fy , 0 , 0 ,
            0 , fx , 0,
            -fy*cx , -fx*cy , fx*fy;
      Vector3d nn = pl_c.normal();
      Vector3d l = K * nn;

      Vector3d _spoint , _epoint;
      _spoint << spoint.x() , spoint.y() , 1.0;
      _epoint << epoint.x() , epoint.y() , 1.0;

      Vector2d obs(_measurement);
      obs << _spoint.transpose() * l , _epoint.transpose() * l;
      double len = sqrt(l.x()*l.x()+l.y()*l.y());
      obs = obs / len;
      _error = obs;
      //std::cout << "_error : " << _error.transpose() << std::endl;

    }

    bool isDepthPositive(){
      VertexSE3Expmap* vj = static_cast<VertexSE3Expmap*>(_vertices[1]);
      SE3Quat T(vj->estimate());
      VertexSBAPluckerLine* vi = static_cast<VertexSBAPluckerLine *>(_vertices[0]);
      PluckerLine pl_w = vi->estimate();

      PluckerLine pl_c = pl_w.transform(T.to_homogeneous_matrix());

      Vector3d pln = pl_c.normal();
      Vector3d plv = pl_c.direction();
      //pln = pln / pln.norm();
      //plv = plv / plv.norm();
      if(plv.norm() < 0.001){
        return false;
      }
/*
      double dd = pln.norm() / plv.norm();
      if(dd < 0.1 || dd > 20){
        return false;
      }

      Vector3d plnv = plv / plv.norm();//pln.cross(plv);
      Vector4d plane_c;
      plane_c.setZero();
      plane_c.segment<3>(0) = plnv;
*/
      Eigen::Matrix4d Lc;
      Lc.setZero();
      Lc.block<3,3>(0,0) = skew(pln);
      Lc.block<3,1>(0,3) = plv;
      Lc.block<1,3>(3,0) = -plv.transpose();

/*
      Vector4d point_zero = Lc * plane_c;

      if(fabs(point_zero[0]) < 0.0001){
        return false;
      }

      Vector3d point_zero_c ;
      point_zero_c << point_zero[0] / point_zero[3] , 
                      point_zero[1] / point_zero[3] ,
                      point_zero[2] / point_zero[3] ;
      if(point_zero_c[2] < 0.3 || point_zero_c[2] > 10){
        return false;
      }

      point_zero_c = point_zero_c / point_zero_c.norm();
      Vector3d nz;
      nz << 0 , 0 , 1;
      const double viewCos = nz.dot(point_zero_c);
      if(viewCos < 0.75 ){
        return false;
      }
*/
      //trimming the endpoints
      Eigen::Vector3d _sp , _ep;
      _sp << spoint.x() , spoint.y() , 1.0;
      _ep << epoint.x() , epoint.y() , 1.0;
      Eigen::Vector3d line_n = _sp.cross(_ep);
      line_n = line_n / line_n.norm();
      Eigen::Vector3d _sp1 , _ep1;
      _sp1 << _sp.x() + 40 * line_n.x() , _sp.y() + 40 * line_n.y() , 1.0;
      _ep1 << _ep.x() + 40 * line_n.x() , _ep.y() + 40 * line_n.y() , 1.0;

      Eigen::Vector3d _sp3d , _ep3d;
      _sp3d << (_sp.x() - cx) / fx , (_sp.y() - cy) / fy , 1.0;
      _ep3d << (_ep.x() - cx) / fx , (_ep.y() - cy) / fy , 1.0;
      Eigen::Vector3d _sp3d1 , _ep3d1;
      _sp3d << (_sp1.x() - cx) / fx , (_sp1.y() - cy) / fy , 1.0;
      _ep3d << (_ep1.x() - cx) / fx , (_ep1.y() - cy) / fy , 1.0;

      Eigen::Vector3d spn , epn;
      spn = _sp3d.cross(_sp3d1);
      spn = spn / spn.norm();
      epn = _ep3d.cross(_ep3d1);
      epn = epn / epn.norm();

      Eigen::Vector4d spl , epl;
      spl.setZero();
      spl.segment<3>(0) = spn;
      epl.setZero();
      epl.segment<3>(0) = epn;

      Eigen::Vector4d sp_ins , ep_ins;
      sp_ins = Lc * spl;
      ep_ins = Lc * epl;

      sp_ins = sp_ins / sp_ins[3];
      ep_ins = ep_ins / ep_ins[3];

      if(sp_ins[2] < 0.2 || ep_ins[2] < 0.2){
        return false;
      }


      return true;
    }

    virtual void linearizeOplus();

    void setParams(const Vector2d& _spoint , const Vector2d& _epoint);
    float fx , fy , cx , cy;
    Vector2d spoint , epoint;

};

class EdgeSE3ProjectPluckerLineOnlyPose : public BaseUnaryEdge <2, Vector2d , VertexSE3Expmap>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  public:
    EdgeSE3ProjectPluckerLineOnlyPose();

    bool read(std::istream& is);
    bool write(std::ostream& os) const;

    void computeError(){
      VertexSE3Expmap* vj = static_cast<VertexSE3Expmap*>(_vertices[0]);
      SE3Quat T(vj->estimate());
      //VertexSBAPluckerLine* vi = static_cast<VertexSBAPluckerLine *>(_vertices[0]);
      //PluckerLine pl_w = vi->estimate();

      PluckerLine pl_c = pl_w.transform(T.to_homogeneous_matrix());

      Eigen::Matrix3d K;
      K << fy , 0 , 0 ,
            0 , fx , 0,
            -fy*cx , -fx*cy , fx*fy;
      Vector3d nn = pl_c.normal();
      Vector3d l = K * nn;

      Vector3d _spoint , _epoint;
      _spoint << spoint.x() , spoint.y() , 1.0;
      _epoint << epoint.x() , epoint.y() , 1.0;

      Vector2d obs(_measurement);
      obs << _spoint.transpose() * l , _epoint.transpose() * l;
      double len = sqrt(l.x()*l.x()+l.y()*l.y());
      obs = obs / len;
      _error = obs;

    }

    virtual void linearizeOplus();
    void setParams(const Vector2d& _spoint , const Vector2d& _epoint);
    void setPluckerLine(const PluckerLine& _plw);

    float fx ,fy , cx,  cy;
    Vector2d spoint , epoint;
    PluckerLine pl_w;

};

class EdgeRelativeSE3 : public BaseBinaryEdge<6, Vector6d , VertexSE3Expmap , VertexSE3Expmap>{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  public:
    EdgeRelativeSE3();

    bool read(std::istream& is);
    bool write(std::ostream& os) const;

    void computeError(){
      VertexSE3Expmap* vj = static_cast<VertexSE3Expmap*>(_vertices[1]);
      SE3Quat Tj(vj->estimate());
      VertexSE3Expmap* vi = static_cast<VertexSE3Expmap*>(_vertices[0]);
      SE3Quat _Ti(vi->estimate());
      SE3Quat _mTr(mTr.block<3,3>(0,0) , mTr.block<3,1>(0,3));

      SE3Quat dT = _Ti * Tj.inverse() * _mTr.inverse();
      //std::cout << "Ti : " << _Ti.to_homogeneous_matrix() << std::endl; 
      // std::cout <<" Tj : " << Tj.to_homogeneous_matrix() << std::endl<<"dT " << dT.to_homogeneous_matrix() << std::endl << "_mTr " << _mTr.to_homogeneous_matrix() << std::endl;
      Vector6d dr = dT.log();
      _error = dr;
      //std::cout << "relative pose error : " << _error.transpose() << std::endl;
      //std::cout << "information : " << information() << std::endl;
      //std::cout << "chi2(): " << error().dot(information()*error()) << std::endl;
    }

    virtual void linearizeOplus();

    void setParams(const Eigen::Matrix4d& _mTr);

    Eigen::Matrix4d mTr;
};

} // end namespace

#endif
