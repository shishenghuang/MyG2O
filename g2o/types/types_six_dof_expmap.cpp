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

#include "types_six_dof_expmap.h"

#include "../core/factory.h"
#include "../stuff/macros.h"

namespace g2o {

using namespace std;


Vector2d project2d(const Vector3d& v)  {
  Vector2d res;
  res(0) = v(0)/v(2);
  res(1) = v(1)/v(2);
  return res;
}

Vector3d unproject2d(const Vector2d& v)  {
  Vector3d res;
  res(0) = v(0);
  res(1) = v(1);
  res(2) = 1;
  return res;
}

VertexSE3Expmap::VertexSE3Expmap() : BaseVertex<6, SE3Quat>() {
}

VertexSE3Expmap::~VertexSE3Expmap(){
}


bool VertexSE3Expmap::read(std::istream& is) {
  Vector7d est;
  for (int i=0; i<7; i++)
    is  >> est[i];
  SE3Quat cam2world;
  cam2world.fromVector(est);
  setEstimate(cam2world.inverse());
  return true;
}

bool VertexSE3Expmap::write(std::ostream& os) const {
  SE3Quat cam2world(estimate().inverse());
  for (int i=0; i<7; i++)
    os << cam2world[i] << " ";
  return os.good();
}

VertexSO3Expmap::VertexSO3Expmap() : BaseVertex<3, SO3>(){

}

bool VertexSO3Expmap::read(std::istream& is){
  Vector3d est;
  for(int i = 0 ;i < 3 ; i ++){
    is >> est[i];

    SO3 plane(SO3::exp(est));
    setEstimate(plane);
    return true;
  }
}

bool VertexSO3Expmap::write(std::ostream& os) const{
  SO3 plane(estimate());
  Vector3d angle = plane.log();
  for(int i= 0; i < 3; i ++){
    os << angle[i] << " ";
  }
  return os.good();
}

EdgeSE3ProjectXYZ::EdgeSE3ProjectXYZ() : BaseBinaryEdge<2, Vector2d, VertexSBAPointXYZ, VertexSE3Expmap>() {
}

bool EdgeSE3ProjectXYZ::read(std::istream& is){
  for (int i=0; i<2; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeSE3ProjectXYZ::write(std::ostream& os) const {

  for (int i=0; i<2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}


void EdgeSE3ProjectXYZ::linearizeOplus() {
  VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
  SE3Quat T(vj->estimate());
  VertexSBAPointXYZ* vi = static_cast<VertexSBAPointXYZ*>(_vertices[0]);
  Vector3d xyz = vi->estimate();
  Vector3d xyz_trans = T.map(xyz);

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double z = xyz_trans[2];
  double z_2 = z*z;

  Matrix<double,2,3> tmp;
  tmp(0,0) = fx;
  tmp(0,1) = 0;
  tmp(0,2) = -x/z*fx;

  tmp(1,0) = 0;
  tmp(1,1) = fy;
  tmp(1,2) = -y/z*fy;

  _jacobianOplusXi =  -1./z * tmp * T.rotation().toRotationMatrix();

  _jacobianOplusXj(0,0) =  x*y/z_2 *fx;
  _jacobianOplusXj(0,1) = -(1+(x*x/z_2)) *fx;
  _jacobianOplusXj(0,2) = y/z *fx;
  _jacobianOplusXj(0,3) = -1./z *fx;
  _jacobianOplusXj(0,4) = 0;
  _jacobianOplusXj(0,5) = x/z_2 *fx;

  _jacobianOplusXj(1,0) = (1+y*y/z_2) *fy;
  _jacobianOplusXj(1,1) = -x*y/z_2 *fy;
  _jacobianOplusXj(1,2) = -x/z *fy;
  _jacobianOplusXj(1,3) = 0;
  _jacobianOplusXj(1,4) = -1./z *fy;
  _jacobianOplusXj(1,5) = y/z_2 *fy;
}

Vector2d EdgeSE3ProjectXYZ::cam_project(const Vector3d & trans_xyz) const{
  Vector2d proj = project2d(trans_xyz);
  Vector2d res;
  res[0] = proj[0]*fx + cx;
  res[1] = proj[1]*fy + cy;
  return res;
}


Vector3d EdgeStereoSE3ProjectXYZ::cam_project(const Vector3d & trans_xyz, const float &bf) const{
  const float invz = 1.0f/trans_xyz[2];
  Vector3d res;
  res[0] = trans_xyz[0]*invz*fx + cx;
  res[1] = trans_xyz[1]*invz*fy + cy;
  res[2] = res[0] - bf*invz;
  return res;
}

EdgeStereoSE3ProjectXYZ::EdgeStereoSE3ProjectXYZ() : BaseBinaryEdge<3, Vector3d, VertexSBAPointXYZ, VertexSE3Expmap>() {
}

bool EdgeStereoSE3ProjectXYZ::read(std::istream& is){
  for (int i=0; i<=3; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeStereoSE3ProjectXYZ::write(std::ostream& os) const {

  for (int i=0; i<=3; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

void EdgeStereoSE3ProjectXYZ::linearizeOplus() {
  VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
  SE3Quat T(vj->estimate());
  VertexSBAPointXYZ* vi = static_cast<VertexSBAPointXYZ*>(_vertices[0]);
  Vector3d xyz = vi->estimate();
  Vector3d xyz_trans = T.map(xyz);

  const Matrix3d R =  T.rotation().toRotationMatrix();

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double z = xyz_trans[2];
  double z_2 = z*z;

  _jacobianOplusXi(0,0) = -fx*R(0,0)/z+fx*x*R(2,0)/z_2;
  _jacobianOplusXi(0,1) = -fx*R(0,1)/z+fx*x*R(2,1)/z_2;
  _jacobianOplusXi(0,2) = -fx*R(0,2)/z+fx*x*R(2,2)/z_2;

  _jacobianOplusXi(1,0) = -fy*R(1,0)/z+fy*y*R(2,0)/z_2;
  _jacobianOplusXi(1,1) = -fy*R(1,1)/z+fy*y*R(2,1)/z_2;
  _jacobianOplusXi(1,2) = -fy*R(1,2)/z+fy*y*R(2,2)/z_2;

  _jacobianOplusXi(2,0) = _jacobianOplusXi(0,0)-bf*R(2,0)/z_2;
  _jacobianOplusXi(2,1) = _jacobianOplusXi(0,1)-bf*R(2,1)/z_2;
  _jacobianOplusXi(2,2) = _jacobianOplusXi(0,2)-bf*R(2,2)/z_2;

  _jacobianOplusXj(0,0) =  x*y/z_2 *fx;
  _jacobianOplusXj(0,1) = -(1+(x*x/z_2)) *fx;
  _jacobianOplusXj(0,2) = y/z *fx;
  _jacobianOplusXj(0,3) = -1./z *fx;
  _jacobianOplusXj(0,4) = 0;
  _jacobianOplusXj(0,5) = x/z_2 *fx;

  _jacobianOplusXj(1,0) = (1+y*y/z_2) *fy;
  _jacobianOplusXj(1,1) = -x*y/z_2 *fy;
  _jacobianOplusXj(1,2) = -x/z *fy;
  _jacobianOplusXj(1,3) = 0;
  _jacobianOplusXj(1,4) = -1./z *fy;
  _jacobianOplusXj(1,5) = y/z_2 *fy;

  _jacobianOplusXj(2,0) = _jacobianOplusXj(0,0)-bf*y/z_2;
  _jacobianOplusXj(2,1) = _jacobianOplusXj(0,1)+bf*x/z_2;
  _jacobianOplusXj(2,2) = _jacobianOplusXj(0,2);
  _jacobianOplusXj(2,3) = _jacobianOplusXj(0,3);
  _jacobianOplusXj(2,4) = 0;
  _jacobianOplusXj(2,5) = _jacobianOplusXj(0,5)-bf/z_2;

  //std::cout << "EdgeStereoSE3ProjectXYZ _jacobianOplusXi: " << std::endl << 
  //_jacobianOplusXi << std::endl << "EdgeStereoSE3ProjectXYZ _jacobianOplusXj:" << _jacobianOplusXj << std::endl;
}


//Only Pose
EdgeSE3ProjectXYZOnlyPose::EdgeSE3ProjectXYZOnlyPose() : BaseUnaryEdge<2, Vector2d, VertexSE3Expmap>(){

}

bool EdgeSE3ProjectXYZOnlyPose::read(std::istream& is){
  for (int i=0; i<2; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeSE3ProjectXYZOnlyPose::write(std::ostream& os) const {

  for (int i=0; i<2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}


void EdgeSE3ProjectXYZOnlyPose::linearizeOplus() {
  VertexSE3Expmap * vi = static_cast<VertexSE3Expmap *>(_vertices[0]);
  Vector3d xyz_trans = vi->estimate().map(Xw);

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double invz = 1.0/xyz_trans[2];
  double invz_2 = invz*invz;

  _jacobianOplusXi(0,0) =  x*y*invz_2 *fx;
  _jacobianOplusXi(0,1) = -(1+(x*x*invz_2)) *fx;
  _jacobianOplusXi(0,2) = y*invz *fx;
  _jacobianOplusXi(0,3) = -invz *fx;
  _jacobianOplusXi(0,4) = 0;
  _jacobianOplusXi(0,5) = x*invz_2 *fx;

  _jacobianOplusXi(1,0) = (1+y*y*invz_2) *fy;
  _jacobianOplusXi(1,1) = -x*y*invz_2 *fy;
  _jacobianOplusXi(1,2) = -x*invz *fy;
  _jacobianOplusXi(1,3) = 0;
  _jacobianOplusXi(1,4) = -invz *fy;
  _jacobianOplusXi(1,5) = y*invz_2 *fy;
}

Vector2d EdgeSE3ProjectXYZOnlyPose::cam_project(const Vector3d & trans_xyz) const{
  Vector2d proj = project2d(trans_xyz);
  Vector2d res;
  res[0] = proj[0]*fx + cx;
  res[1] = proj[1]*fy + cy;
  return res;
}

EdgeStereoSE3ProjectXYZOnlyPose::EdgeStereoSE3ProjectXYZOnlyPose():BaseUnaryEdge<3, Vector3d, VertexSE3Expmap>(){

}

Vector3d EdgeStereoSE3ProjectXYZOnlyPose::cam_project(const Vector3d & trans_xyz) const{
  const float invz = 1.0f/trans_xyz[2];
  Vector3d res;
  res[0] = trans_xyz[0]*invz*fx + cx;
  res[1] = trans_xyz[1]*invz*fy + cy;
  res[2] = res[0] - bf*invz;
  return res;
}


bool EdgeStereoSE3ProjectXYZOnlyPose::read(std::istream& is){
  for (int i=0; i<=3; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeStereoSE3ProjectXYZOnlyPose::write(std::ostream& os) const {

  for (int i=0; i<=3; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

void EdgeStereoSE3ProjectXYZOnlyPose::linearizeOplus() {
  VertexSE3Expmap * vi = static_cast<VertexSE3Expmap *>(_vertices[0]);
  Vector3d xyz_trans = vi->estimate().map(Xw);

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double invz = 1.0/xyz_trans[2];
  double invz_2 = invz*invz;

  _jacobianOplusXi(0,0) =  x*y*invz_2 *fx;
  _jacobianOplusXi(0,1) = -(1+(x*x*invz_2)) *fx;
  _jacobianOplusXi(0,2) = y*invz *fx;
  _jacobianOplusXi(0,3) = -invz *fx;
  _jacobianOplusXi(0,4) = 0;
  _jacobianOplusXi(0,5) = x*invz_2 *fx;

  _jacobianOplusXi(1,0) = (1+y*y*invz_2) *fy;
  _jacobianOplusXi(1,1) = -x*y*invz_2 *fy;
  _jacobianOplusXi(1,2) = -x*invz *fy;
  _jacobianOplusXi(1,3) = 0;
  _jacobianOplusXi(1,4) = -invz *fy;
  _jacobianOplusXi(1,5) = y*invz_2 *fy;

  _jacobianOplusXi(2,0) = _jacobianOplusXi(0,0)-bf*y*invz_2;
  _jacobianOplusXi(2,1) = _jacobianOplusXi(0,1)+bf*x*invz_2;
  _jacobianOplusXi(2,2) = _jacobianOplusXi(0,2);
  _jacobianOplusXi(2,3) = _jacobianOplusXi(0,3);
  _jacobianOplusXi(2,4) = 0;
  _jacobianOplusXi(2,5) = _jacobianOplusXi(0,5)-bf*invz_2;
}

EdgeSE3ProjectLineXYZ::EdgeSE3ProjectLineXYZ() : BaseBinaryEdge<2, Vector2d , VertexSBALineXYZ , VertexSE3Expmap >(){

}

bool EdgeSE3ProjectLineXYZ::read(std::istream& is){
  for(int i=0; i < 2; i ++){
    is >> _measurement[i];
  }
  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
  }

  return true;  
}

bool EdgeSE3ProjectLineXYZ::write(std::ostream& os) const{
    for (int i=0; i<2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

void EdgeSE3ProjectLineXYZ::setParams(const Vector2d& _pu , const Vector2d& _qu){
  pu = _pu;
  qu = _qu;  
}

void EdgeSE3ProjectLineXYZ::linearizeOplus(){
  VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
  SE3Quat T(vj->estimate());
  VertexSBALineXYZ *vi = static_cast<VertexSBALineXYZ *>(_vertices[0]);
  Vector6d line_p_q = vi->estimate();
  Vector3d line_p;//
  line_p << line_p_q[0] , line_p_q[1] , line_p_q[2];
  Vector3d line_q;
  line_q << line_p_q[3] , line_p_q[4] , line_p_q[5]; 

  //////////////////////////////////////////////////////////////////////////////////////////
  Vector3d h_pu; 
  h_pu << (pu[0] , pu[1] , 1.0);
  Vector3d h_qu;
  h_qu << qu[0] , qu[1] , 1.0;
  Vector3d l_pu_qu = h_pu.cross(h_qu);
  double l_pu_qu_n = sqrt(l_pu_qu[0]*l_pu_qu[0]+l_pu_qu[1]*l_pu_qu[1]);
  l_pu_qu = l_pu_qu / l_pu_qu_n;//l_pu_qu.norm();
  Vector2d l_pu_qu_v;
  l_pu_qu_v << l_pu_qu[0] , l_pu_qu[1];

  Vector3d line_p_trans = T.map(line_p);
  Vector3d line_q_trans = T.map(line_q);

  double px = line_p_trans[0];
  double py = line_p_trans[1];
  double pz = line_p_trans[2];
  Eigen::Matrix<double, 2,3> h_p = Eigen::Matrix<double , 2,3>::Zero();
  h_p << fx / pz , 0 , -fx * px / (pz * pz),
          0, fy / pz , -fy * py / (pz * pz);
          //0, 0, 1.0;
  
  double qx = line_q_trans[0];
  double qy = line_q_trans[1];
  double qz = line_q_trans[2];
  Eigen::Matrix<double, 2,3> h_q = Eigen::Matrix<double , 2,3>::Zero();
  h_q << fx / qz , 0 , -fx * qx / (qz * qz),
          0 , fy / qz , -fy * qy / (qz * qz);
          //0, 0, 1.0;
  
  Eigen::Matrix<double, 2,6> de_dpq = Eigen::Matrix<double, 2,6>::Zero();
  de_dpq.block<1,3>(0,0) = l_pu_qu_v.transpose() * h_p * T.rotation().toRotationMatrix();
  de_dpq.block<1,3>(1,3) = l_pu_qu_v.transpose() * h_q * T.rotation().toRotationMatrix();

  /////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////

  Eigen::Matrix3d p_skew = SO3::hat(line_p_trans) ; // skew_p 
  Eigen::Matrix<double, 3,6> p_skew_and_I = Eigen::Matrix<double, 3,6>::Zero();
  p_skew_and_I.block<3,3>(0,0) = - p_skew;
  p_skew_and_I.block<3,3>(0,3) = Eigen::Matrix3d::Identity();

  Eigen::Matrix3d q_skew = SO3::hat(line_q_trans); // skew_q
  Eigen::Matrix<double, 3,6> q_skew_and_I = Eigen::Matrix<double, 3,6>::Zero();
  q_skew_and_I.block<3,3>(0,0) = - q_skew;
  q_skew_and_I.block<3,3>(0,3) = Eigen::Matrix3d::Identity();

  Eigen::Matrix<double, 2,6> de_dT = Eigen::Matrix<double, 2,6>::Zero();
  de_dT.block<1,6>(0,0) = l_pu_qu_v.transpose() * h_p * p_skew_and_I;
  de_dT.block<1,6>(1,0) = l_pu_qu_v.transpose() * h_q * q_skew_and_I;
  /////////////////////////////////////////////////////////////////////////////////////////


  _jacobianOplusXi = de_dpq;
  _jacobianOplusXj = de_dT;
  //std::cout << "EdgeSE3ProjectLineXYZ _jacobianOplusXi : " << std::endl << _jacobianOplusXi << std::endl;
  //std::cout << "EdgeSE3ProjectLineXYZ _jacobianOplusXj : " << std::endl << _jacobianOplusXj << std::endl;


  // ok, finish

}

EdgeSE3ProjectLineXYZOnlyPose::EdgeSE3ProjectLineXYZOnlyPose() : BaseUnaryEdge<2, Vector2d , VertexSE3Expmap >(){

}

bool EdgeSE3ProjectLineXYZOnlyPose::read(std::istream& is){
  for(int i=0; i < 2; i ++){
    is >> _measurement[i];
  }
  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
  }

  return true;  
}

bool EdgeSE3ProjectLineXYZOnlyPose::write(std::ostream& os) const{
    for (int i=0; i<2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}


void EdgeSE3ProjectLineXYZOnlyPose::linearizeOplus(){
  VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[0]);
  SE3Quat T(vj->estimate());
  //VertexSBALineXYZ *vi = static_cast<VertexSBALineXYZ *>(_vertices[0]);
  //Vector6d line_p_q = vi->estimate();
  Vector3d line_p; //
  line_p << line_p_q[0] , line_p_q[1] , line_p_q[2];
  Vector3d line_q;
  line_q << line_p_q[3] , line_p_q[4] , line_p_q[5]; 

  //////////////////////////////////////////////////////////////////////////////////////////
  Vector3d h_pu;
  h_pu << pu[0] , pu[1] , 1.0;
  Vector3d h_qu;
  h_qu << qu[0] , qu[1] , 1.0;
  Vector3d l_pu_qu = h_pu.cross(h_qu);
  double l_pu_qu_n = sqrt(l_pu_qu[0]*l_pu_qu[0]+l_pu_qu[1]*l_pu_qu[1]);
  l_pu_qu = l_pu_qu / l_pu_qu_n;// l_pu_qu.norm();
  Vector2d l_pu_qu_v;
  l_pu_qu_v << l_pu_qu[0] , l_pu_qu[1];

  Vector3d line_p_trans = T.map(line_p);
  Vector3d line_q_trans = T.map(line_q);

  double px = line_p_trans[0];
  double py = line_p_trans[1];
  double pz = line_p_trans[2];
  Eigen::Matrix<double, 2,3> h_p = Eigen::Matrix<double , 2,3>::Zero();
  h_p << fx / pz , 0 , -fx * px / (pz * pz),
          0, fy / pz , -fy * py / (pz * pz);
          //0, 0, 1.0;
  
  double qx = line_q_trans[0];
  double qy = line_q_trans[1];
  double qz = line_q_trans[2];
  Eigen::Matrix<double, 2,3> h_q = Eigen::Matrix<double , 2,3>::Zero();
  h_q << fx / qz , 0 , -fx * qx / (qz * qz),
          0 , fy / qz , -fy * qy / (qz * qz);
          //0, 0, 1.0;
  
  Eigen::Matrix<double, 2,6> de_dpq = Eigen::Matrix<double, 2,6>::Zero();
  de_dpq.block<1,3>(0,0) = l_pu_qu_v.transpose() * h_p * T.rotation().toRotationMatrix();
  de_dpq.block<1,3>(1,3) = l_pu_qu_v.transpose() * h_q * T.rotation().toRotationMatrix();

  /////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////
  Eigen::Matrix3d p_skew = SO3::hat(line_p_trans) ; // skew_p 
  Eigen::Matrix<double, 3,6> p_skew_and_I = Eigen::Matrix<double, 3,6>::Zero();
  p_skew_and_I.block<3,3>(0,0) = - p_skew;
  p_skew_and_I.block<3,3>(0,3) = Eigen::Matrix3d::Identity();

  Eigen::Matrix3d q_skew = SO3::hat(line_q_trans); // skew_q
  Eigen::Matrix<double, 3,6> q_skew_and_I = Eigen::Matrix<double, 3,6>::Zero();
  q_skew_and_I.block<3,3>(0,0) = - q_skew;
  q_skew_and_I.block<3,3>(0,3) = Eigen::Matrix3d::Identity();

  Eigen::Matrix<double, 2,6> de_dT = Eigen::Matrix<double, 2,6>::Zero();
  de_dT.block<1,6>(0,0) = l_pu_qu_v.transpose() * h_p * p_skew_and_I;
  de_dT.block<1,6>(1,0) = l_pu_qu_v.transpose() * h_q * q_skew_and_I;
  /////////////////////////////////////////////////////////////////////////////////////////

  //_jacobianOplusXi = de_dpq;
  _jacobianOplusXi = de_dT;

  // ok, finish
}

void EdgeSE3ProjectLineXYZOnlyPose::setParams(const Vector2d& _pu , const Vector2d& _qu){
  pu = _pu;
  qu = _qu;
}

void EdgeSE3ProjectLineXYZOnlyPose::setParams(const Vector4f& puqu){
  pu << puqu[0] , puqu[1];
  qu << puqu[2] , puqu[3];
}

void EdgeSE3ProjectLineXYZOnlyPose::setMeasurement3DLine(const Vector6d& _line_p_q){
  line_p_q = _line_p_q;
}

EdgeSE3ProjectPlaneSO3::EdgeSE3ProjectPlaneSO3() : BaseBinaryEdge<3, Vector3d , VertexSO3Expmap , VertexSE3Expmap>(){

}

bool EdgeSE3ProjectPlaneSO3::read(std::istream& is){
  for(int i=0; i <= 2; i ++){
    is >> _measurement[i];
  }
  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }

  return true;  
}

bool EdgeSE3ProjectPlaneSO3::write(std::ostream& os) const{
    for (int i=0; i<=2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

void EdgeSE3ProjectPlaneSO3::linearizeOplus(){

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

  Vector3d q_y_v = n_q_y.vec();
  double q_y_v_norm = q_y_v.norm();
  double q_y_v_norm_3 = q_y_v_norm * q_y_v_norm * q_y_v_norm;
  double q4 = n_q_y.w();

  Eigen::Matrix3d J_l_u_inv = SO3::JacobianLInv(u);

  Eigen::Matrix<double ,3,4> M;
  M.block<3,3>(0,0) = 2 * acos(q4) * ( Eigen::Matrix3d::Identity() / q_y_v_norm - q_y_v * q_y_v.transpose()/q_y_v_norm_3);
  M.block<3,1>(0,3) = 2.0 * q_y_v / (q_y_v_norm * sqrt(1.0 - q4 * q4 + 0.0001));

  Eigen::Matrix<double, 4,4> Q;
  double y_norm = y.norm();
  double y_norm_3 = y_norm * y_norm * y_norm;
  Q = Eigen::Matrix4d::Identity() / y_norm - y * y.transpose() / y_norm_3;

  Vector3d pi_x_v;
  pi_x_v << pi_x.x() , pi_x.y() , pi_x.z();
  Eigen::Matrix<double, 4,6> H = Eigen::Matrix<double, 4,6>::Zero();
  H.block<3,3>(0,0) = T.rotation().toRotationMatrix().transpose() * SO3::hat(pi_x_v);
  H.block<1,3>(3,0) = T.translation().transpose() * SO3::hat(pi_x_v);
  H.block<1,3>(3,3) = pi_x_v.transpose();
  //////////////////////////////////////////////////////////////////////////////////////////

  // ok, finish
  _jacobianOplusXi = SO3::JacobianRInv(u); // Jr(u)^(-1)
  _jacobianOplusXj = J_l_u_inv * M * Q * H; // 3 x 6 matrix
  //////////////////////////////////////////////////////////////////////////////////////////

}
void EdgeSE3ProjectPlaneSO3::setParams(const SO3& plane_){
  plane_c = plane_;
}

EdgeSE3ProjectPlaneSO3OnlyPose::EdgeSE3ProjectPlaneSO3OnlyPose() : BaseUnaryEdge<3, Vector3d , VertexSE3Expmap>(){

}

bool EdgeSE3ProjectPlaneSO3OnlyPose::read(std::istream& is){
  for(int i=0; i <= 2; i ++){
    is >> _measurement[i];
  }
  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }

  return true;  
}

bool EdgeSE3ProjectPlaneSO3OnlyPose::write(std::ostream& os) const{
    for (int i=0; i<=2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

void EdgeSE3ProjectPlaneSO3OnlyPose::linearizeOplus(){

  VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[0]);
  SE3Quat T(vj->estimate());
  //VertexSO3Expmap *vi = static_cast<VertexSO3Expmap *>(_vertices[0]);
  SO3 plane = plane_g;//(vi->estimate());

  Quaterniond pi_x = plane_c.unit_quaternion();      // the observation
  Vector4d pi_x_4d;
  pi_x_4d << pi_x.x() , pi_x.y() , pi_x.z() , pi_x.w();
  
  Vector4d y = T.to_homogeneous_matrix().transpose() * pi_x_4d;
  Vector4d q_y = y / y.norm();
  Quaterniond n_q_y(q_y[3] , q_y[0] , q_y[1] , q_y[2]);
  SO3 q_T(n_q_y);
  SO3 q_u = q_T * plane.inverse();
  Vector3d u = q_u.log();

  Vector3d q_y_v = n_q_y.vec();
  double q_y_v_norm = q_y_v.norm();
  double q_y_v_norm_3 = q_y_v_norm * q_y_v_norm * q_y_v_norm;
  double q4 = n_q_y.w();

  Eigen::Matrix3d J_l_u_inv = SO3::JacobianLInv(u);

  Eigen::Matrix<double ,3,4> M;
  M.block<3,3>(0,0) = 2 * acos(q4) * ( Eigen::Matrix3d::Identity() / q_y_v_norm - q_y_v * q_y_v.transpose()/q_y_v_norm_3);
  M.block<3,1>(0,3) = 2.0 * q_y_v / (q_y_v_norm * sqrt(1.0 - q4 * q4 + 0.0001));

  Eigen::Matrix<double, 4,4> Q;
  double y_norm = y.norm();
  double y_norm_3 = y_norm * y_norm * y_norm;
  Q = Eigen::Matrix4d::Identity() / y_norm - y * y.transpose() / y_norm_3;

  Vector3d pi_x_v;
  pi_x_v << pi_x.x() , pi_x.y() , pi_x.z();
  Eigen::Matrix<double, 4,6> H = Eigen::Matrix<double, 4,6>::Zero();
  H.block<3,3>(0,0) = T.rotation().toRotationMatrix().transpose() * SO3::hat(pi_x_v);
  H.block<1,3>(3,0) = T.translation().transpose() * SO3::hat(pi_x_v);
  H.block<1,3>(3,3) = pi_x_v.transpose();
  //////////////////////////////////////////////////////////////////////////////////////////

  // ok, finish
  //_jacobianOplusXi = SO3::JacobianRInv(u); // Jr(u)^(-1)
  _jacobianOplusXi = J_l_u_inv * M * Q * H; // 3 x 6 matrix
  //////////////////////////////////////////////////////////////////////////////////////////

}
void EdgeSE3ProjectPlaneSO3OnlyPose::setParams(const SO3& plane_, const SO3& _plane_g){
  plane_c = plane_;
  plane_g = _plane_g;
}


EdgeSE3ProjectPluckerLine::EdgeSE3ProjectPluckerLine() : BaseBinaryEdge<2 , Vector2d , VertexSBAPluckerLine , VertexSE3Expmap>(){

}
void EdgeSE3ProjectPluckerLine::setParams(const Vector2d& _spoint , const Vector2d& _epoint){
  spoint = _spoint;
  epoint = _epoint;
}

void EdgeSE3ProjectPluckerLine::linearizeOplus(){
  VertexSE3Expmap* vj = static_cast<VertexSE3Expmap*>(_vertices[1]);
  SE3Quat Tcw(vj->estimate());
  VertexSBAPluckerLine *vi = static_cast<VertexSBAPluckerLine*>(_vertices[0]);
  PluckerLine pl_w = vi->estimate();

  Matrix3d Rcw = Tcw.rotation().toRotationMatrix();
  Vector3d tcw = Tcw.translation();

  PluckerLine pl_c = pl_w.transform(Tcw.to_homogeneous_matrix());
  Matrix3d K;
  K << fy , 0  , 0,
        0 , fx , 0,
        -fy*cx,-fx*cy,fx*fy;
  Vector3d l = K * pl_c.normal();
  Vector3d _spoint , _epoint;
  _spoint << spoint.x() , spoint.y() , 1.0;
  _epoint << epoint.x() , epoint.y() , 1.0;
  double es , ee;
  es = _spoint.transpose() * l;
  ee = _epoint.transpose() * l;

  double ln = sqrt(l.x()*l.x()+l.y()*l.y());

  Eigen::Matrix<double , 2, 3> del_dl;
  del_dl << _spoint.x() - l.x() * es / (ln * ln) , _spoint.y() - l.y() * es / (ln * ln) , 1.0,
            _epoint.x() - l.x() * ee / (ln * ln) , _epoint.y() - l.y() * ee / (ln * ln) , 1.0;
  del_dl = del_dl / ln;
  Eigen::Matrix<double , 3, 6> dl_dLc;
  dl_dLc.setZero();
  dl_dLc.block<3,3>(0,0) = K;

  Eigen::Matrix<double , 6, 6> dLc_dLw;
  dLc_dLw.setZero();
  dLc_dLw.block<3,3>(0,0) = Rcw;
  dLc_dLw.block<3,3>(3,3) = Rcw;
  dLc_dLw.block<3,3>(0,3) = skew(tcw) * Rcw;

  Eigen::Matrix<double , 6, 4> dLw_dthelta;
  Matrix3d U_w = pl_w.Umatrix();
  Matrix<double ,3, 2> W_w = pl_w.Wmatrix();
  double w1 = W_w(0,0);
  double w2 = W_w(1,1);

  dLw_dthelta.block<3,3>(0,0) = - skew(w1 * U_w.block<3,1>(0,0));
  dLw_dthelta.block<3,1>(0,3) = -w2 * U_w.block<3,1>(0,0);
  dLw_dthelta.block<3,3>(3,0) = - skew(w2 * U_w.block<3,1>(0,1));
  dLw_dthelta.block<3,1>(3,3) = w1 * U_w.block<3,1>(0,1);

  // jacobian of xi
  _jacobianOplusXi = del_dl * dl_dLc * dLc_dLw * dLw_dthelta;

  Matrix<double , 6, 6> dLc_dxi;
  dLc_dxi.setZero();
  Vector3d n_w = pl_w.normal();
  Vector3d v_w = pl_w.direction();
  dLc_dxi.block<3,3>(0,0) = - skew(Rcw * n_w) - skew(skew(tcw)*Rcw*v_w);
  dLc_dxi.block<3,3>(0,3) = - skew(Rcw * v_w);
  dLc_dxi.block<3,3>(3,0) = - skew(Rcw * v_w);

  // jacobian of xj
  _jacobianOplusXj = del_dl * dl_dLc * dLc_dxi;

  //std::cout << "_jacobianOplusXi : " << _jacobianOplusXi << std::endl;
  //std::cout << "_jacobianOplusXj : " << _jacobianOplusXj << std::endl;
}

void EdgeParallel::linearizeOplus() {

    VertexSBAPluckerLine *vi0 = static_cast<VertexSBAPluckerLine *>(_vertices[0]), *vi1 = static_cast<VertexSBAPluckerLine *>(_vertices[1]);
    PluckerLine pl_w0 = vi0->estimate(), pl_w1 = vi1->estimate();
    Eigen::Matrix<double , 6, 4> dLw_dthelta0;
    
    {
      Matrix3d U_w = pl_w0.Umatrix();
      Matrix<double ,3, 2> W_w = pl_w0.Wmatrix();
      double w1 = W_w(0,0);
      double w2 = W_w(1,1);

      dLw_dthelta0.block<3,3>(0,0) = - skew(w1 * U_w.block<3,1>(0,0));
      dLw_dthelta0.block<3,1>(0,3) = -w2 * U_w.block<3,1>(0,0);
      dLw_dthelta0.block<3,3>(3,0) = - skew(w2 * U_w.block<3,1>(0,1));
      dLw_dthelta0.block<3,1>(3,3) = w1 * U_w.block<3,1>(0,1);
      Eigen::Matrix<double , 3, 6> dE_dLw;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++)
          dE_dLw(i, j) = 0;
      }
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) if (j != i)
        {
          int k = 3 - i - j;
          if ((j - i + 3) % 3 == 1) dE_dLw(i, 3 + j) = pl_w1.direction()[k];
          else dE_dLw(i, 3 + j) = -pl_w1.direction()[k];
        }
      _jacobianOplusXi = dE_dLw * dLw_dthelta0;
    }

    
    Eigen::Matrix<double , 6, 4> dLw_dthelta1;
    {
      Matrix3d U_w = pl_w1.Umatrix();
      Matrix<double ,3, 2> W_w = pl_w1.Wmatrix();
      double w1 = W_w(0,0);
      double w2 = W_w(1,1);

      dLw_dthelta1.block<3,3>(0,0) = - skew(w1 * U_w.block<3,1>(0,0));
      dLw_dthelta1.block<3,1>(0,3) = -w2 * U_w.block<3,1>(0,0);
      dLw_dthelta1.block<3,3>(3,0) = - skew(w2 * U_w.block<3,1>(0,1));
      dLw_dthelta1.block<3,1>(3,3) = w1 * U_w.block<3,1>(0,1);
      Eigen::Matrix<double , 3, 6> dE_dLw;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++)
          dE_dLw(i, j) = 0;
      }
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) if (j != i)
        {
          int k = 3 - i - j;
          if ((j - i + 3) % 3 == 1) dE_dLw(i, 3 + j) = -pl_w0.direction()[k];
          else dE_dLw(i, 3 + j) = pl_w1.direction()[k];
        }
      _jacobianOplusXj = dE_dLw * dLw_dthelta1;
      // for (int i = 0; i < 3; i++) {
      //   for (int j = 0; j < 6; j++)
      //     cout << dE_dLw(i, j) << " ";
      //   cout << endl;
      // }
    }

}

bool EdgeSE3ProjectPluckerLine::write(std::ostream& os) const{
    for (int i=0; i<2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

bool EdgeSE3ProjectPluckerLine::read(std::istream& is){
  for(int i=0; i < 2; i ++){
    is >> _measurement[i];
  }
  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
  }

  return true;  
}

EdgeSE3ProjectPluckerLineOnlyPose::EdgeSE3ProjectPluckerLineOnlyPose() : BaseUnaryEdge<2 , Vector2d , VertexSE3Expmap>()
{
  
}

bool EdgeSE3ProjectPluckerLineOnlyPose::write(std::ostream& os) const{
    for (int i=0; i<2; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

bool EdgeSE3ProjectPluckerLineOnlyPose::read(std::istream& is){
  for(int i=0; i < 2; i ++){
    is >> _measurement[i];
  }
  for (int i=0; i<2; i++)
    for (int j=i; j<2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
  }

  return true;  
}

void EdgeSE3ProjectPluckerLineOnlyPose::linearizeOplus(){
  VertexSE3Expmap* vj = static_cast<VertexSE3Expmap*>(_vertices[0]);
  SE3Quat Tcw(vj->estimate());
  //VertexSBAPluckerLine *vi = static_cast<VertexSBAPluckerLine*>(_vertices[0]);
  //PluckerLine pl_w = ;//vi->estimate();

  Matrix3d Rcw = Tcw.rotation().toRotationMatrix();
  Vector3d tcw = Tcw.translation();

  PluckerLine pl_c = pl_w.transform(Tcw.to_homogeneous_matrix());
  Matrix3d K;
  K << fy , 0  , 0,
        0 , fx , 0,
        -fy*cx,-fx*cy,fx*fy;
  Vector3d l = K * pl_c.normal();
  Vector3d _spoint , _epoint;
  _spoint << spoint.x() , spoint.y() , 1.0;
  _epoint << epoint.x() , epoint.y() , 1.0;
  double es , ee;
  es = _spoint.transpose() * l;
  ee = _epoint.transpose() * l;

  double ln = sqrt(l.x()*l.x()+l.y()*l.y());

  Eigen::Matrix<double , 2, 3> del_dl;
  del_dl << _spoint.x() - l.x() * es / (ln * ln) , _spoint.y() - l.y() * es / (ln * ln) , 1.0,
            _epoint.x() - l.x() * ee / (ln * ln) , _epoint.y() - l.y() * ee / (ln * ln) , 1.0;
  del_dl = del_dl / ln;
  Eigen::Matrix<double , 3, 6> dl_dLc;
  dl_dLc.setZero();
  dl_dLc.block<3,3>(0,0) = K;

/*
  Eigen::Matrix<double , 6, 6> dLc_dLw;
  dLc_dLw.setZero();
  dLc_dLw.block<3,3>(0,0) = Rcw;
  dLc_dLw.block<3,3>(3,3) = Rcw;
  dLc_dLw.block<3,3>(0,3) = skew(tcw) * Rcw;

  Eigen::Matrix<double , 6, 4> dLw_dthelta;
  Matrix3d U_w = pl_w.Umatrix();
  Matrix<double ,3, 2> W_w = pl_w.Wmatrix();
  double w1 = W_w(0,0);
  double w2 = W_w(1,1);

  dLw_dthelta.block<3,3>(0,0) = - skew(w1 * U_w.block<3,1>(0,0));
  dLw_dthelta.block<3,1>(0,3) = -w2 * U_w.block<3,1>(0,0);
  dLw_dthelta.block<3,3>(3,0) = - skew(w2 * U_w.block<3,1>(0,1));
  dLw_dthelta.block<3,1>(3,3) = w1 * U_w.block<3,1>(0,1);

  // jacobian of xi
  _jacobianOplusXi = del_dl * dl_dLc * dLc_dLw * dLw_dthelta;
*/

  Matrix<double , 6, 6> dLc_dxi;
  dLc_dxi.setZero();
  Vector3d n_w = pl_w.normal();
  Vector3d v_w = pl_w.direction();
  dLc_dxi.block<3,3>(0,0) = - skew(Rcw * n_w) - skew(skew(tcw)*Rcw*v_w);
  dLc_dxi.block<3,3>(0,3) = - skew(Rcw * v_w);
  dLc_dxi.block<3,3>(3,0) = - skew(Rcw * v_w);

  // jacobian of xj
  _jacobianOplusXi = del_dl * dl_dLc * dLc_dxi;

}

void EdgeSE3ProjectPluckerLineOnlyPose::setParams(const Vector2d& _spoint , const Vector2d& _epoint){
  spoint = _spoint;
  epoint = _epoint;
}

void EdgeSE3ProjectPluckerLineOnlyPose::setPluckerLine(const PluckerLine& _plw){
  pl_w = _plw;
}

EdgeRelativeSE3::EdgeRelativeSE3() : BaseBinaryEdge<6, Vector6d , VertexSE3Expmap , VertexSE3Expmap>(){

}

bool EdgeRelativeSE3::write(std::ostream& os) const{
    for (int i=0; i<6; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<6; i++)
    for (int j=i; j<6; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

bool EdgeRelativeSE3::read(std::istream& is){
  for(int i=0; i < 6; i ++){
    is >> _measurement[i];
  }
  for (int i=0; i<6; i++)
    for (int j=i; j<6; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
  }

  return true;  
}

void EdgeRelativeSE3::linearizeOplus(){
    VertexSE3Expmap* vj = static_cast<VertexSE3Expmap*>(_vertices[1]);
    SE3Quat Tj(vj->estimate());
    VertexSE3Expmap* vi = static_cast<VertexSE3Expmap*>(_vertices[0]);
    SE3Quat _Ti(vi->estimate());
    SE3Quat _mTr(mTr.block<3,3>(0,0) , mTr.block<3,1>(0,3));

/*  2019-03-13改动
    SE3Quat dT = _Ti*Tj.inverse()*_mTr.inverse();
    Eigen::Quaterniond dr = dT.rotation();
    Eigen::Vector3d dt = dT.translation();
    SO3 dT_r(dr);
    Eigen::Matrix3d Jl_r_inv = SO3::JacobianLInv(dT_r.log());

    SE3Quat dP = _Ti*Tj.inverse();
    Eigen::Quaterniond dr_p = dP.rotation();
    Eigen::Matrix3d dR = dr_p.toRotationMatrix();
    Eigen::Vector3d dt_p = dP.translation();

    SE3Quat dM = _mTr.inverse();
    Eigen::Quaterniond dm_r = dM.rotation();
    Eigen::Vector3d dt_r = dM.translation();

    Eigen::Matrix<double , 6, 6> dT_dTi, dT_dTj;
    dT_dTi.setZero();
    dT_dTi.block<3,3>(0,0) = Jl_r_inv;
    dT_dTi.block<3,3>(3,0) = - Jl_r_inv * skew(dt);
    dT_dTi.block<3,3>(3,3) = Jl_r_inv;

    dT_dTj.setZero();
    dT_dTj.block<3,3>(0,0) = - Jl_r_inv * dR;
    dT_dTj.block<3,3>(3,0) = Jl_r_inv * dR * skew(dt_r);
    dT_dTj.block<3,3>(3,3) = - Jl_r_inv * dR;
    //////////////////////////////////////////////////
    // compute the Jacobian
    _jacobianOplusXi = dT_dTi;
    _jacobianOplusXj = dT_dTj;
    //std::cout << "Ji : " << std::endl << _jacobianOplusXi << std::endl << " Jj: " << std::endl << _jacobianOplusXj << std::endl;
*/

  SE3Quat Tij = _Ti*Tj.inverse()*_mTr.inverse();
  Eigen::Quaterniond rij = Tij.rotation();
  Eigen::Vector3d tij = Tij.translation();
  SO3 Rij(rij);
  Eigen::Matrix3d J_l_inv = SO3::JacobianLInv(Rij.log());
  Eigen::Matrix3d J_r_inv = SO3::JacobianRInv(Rij.log());

  Eigen::Matrix4d dTij = _mTr.to_homogeneous_matrix();
  Eigen::Matrix3d dRij = dTij.block<3,3>(0,0);
  Eigen::Vector3d dtij = dTij.block<3,1>(0,3);

  Eigen::Matrix<double , 6, 6> dT_dTi , dT_dTj;
  dT_dTi.setZero();
  dT_dTi.block<3,3>(0,0) = J_l_inv;
  dT_dTi.block<3,3>(3,0) = - J_l_inv * skew(tij);
  dT_dTi.block<3,3>(3,3) = J_l_inv;

  dT_dTj.setZero();
  dT_dTj.block<3,3>(0,0) = -J_r_inv * dRij;
  dT_dTj.block<3,3>(3,0) = -J_r_inv * skew(dtij) * dRij;
  dT_dTj.block<3,3>(3,3) = -J_r_inv * dRij;

  _jacobianOplusXi = dT_dTi;
  _jacobianOplusXj = dT_dTj; 

}

void EdgeRelativeSE3::setParams(const Eigen::Matrix4d& _mTr){
  mTr = _mTr;
}


} // end namespace
