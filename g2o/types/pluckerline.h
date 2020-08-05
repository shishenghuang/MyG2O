
/**
 * This is implementation of Plucker Line descriped in 
 * G. Zhang, J. H. Lee, J. Lim and I. H. Suh, "Building a 3-D Line-Based Map Using Stereo SLAM," in IEEE Transactions on Robotics, vol. 31, no. 6, pp. 1364-1377, Dec. 2015.
 * more details please refer to this paper
 * Author: Shi-Sheng Huang
 **/

#ifndef _PLUCKERLINE_H
#define _PLUCKERLINE_H

#include <iostream>
#include <vector>
#include "se3quat.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "se3_ops.h"

namespace g2o{

    using namespace Eigen;
    const double SMALL_EPS_LINE = 1e-5;

    class PluckerLine{
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        protected:
            Eigen::Vector3d _n;
            Eigen::Vector3d _v;
        public:
            PluckerLine(){
                _n.setZero();
                _v.setZero();
            }
            PluckerLine(const Eigen::Vector3d& _nn , const Eigen::Vector3d& _vv){
                _n = _nn;
                _v = _vv;
                //normalize();
            }
            ~PluckerLine(){
                
            }

            inline const Vector3d& normal() const { return _n;}
            inline const Vector3d& direction() const { return _v;}

            void normalize(){
                double len = _n.norm();
                _v = _v / len;
                _n = _n / len;
            }

            inline PluckerLine operator*(const PluckerLine& pl){
                // L = this * pl
                Vector3d nn1 = _n / _n.norm();
                Vector3d vv1 = _v / _v.norm();
                Vector3d nv1 = nn1.cross(vv1);//_v.cross(_n);
                //nv1 = nv1 / nv1.norm();

                Eigen::Matrix3d U1;
                U1.block<3,1>(0,0) = nn1;
                U1.block<3,1>(0,1) = vv1;
                U1.block<3,1>(0,2) = nv1;
                Eigen::Vector2d sigma1;
                sigma1 << _n.norm() , _v.norm();
                sigma1 = sigma1 / sigma1.norm();
                Eigen::Matrix2d W1;
                W1 << sigma1.x() , - sigma1.y(),
                        sigma1.y() , sigma1.x();

                Vector3d pln = pl.normal();
                Vector3d plv = pl.direction();
                Vector3d nn2 = pln / pln.norm();
                Vector3d vv2 = plv / plv.norm();
                Vector3d nv2 = vv2.cross(nn2);

                Eigen::Matrix3d U2;
                U2.block<3,1>(0,0) = vv2;
                U2.block<3,1>(0,1) = nn2;
                U2.block<3,1>(0,2) = nv2;
                Eigen::Vector2d sigma2;
                sigma2 << pln.norm() , plv.norm();
                sigma2 = sigma2 / sigma2.norm();
                Eigen::Matrix2d W2;
                W2 << sigma2.x() , -sigma2.y(),
                        sigma2.y() , sigma2.x();
                
                Eigen::Matrix3d U = U1 * U2;
                Eigen::Matrix2d W = W1 * W2;

                Eigen::Vector3d n_n = W(0,0) * U.block<3,1>(0,0);
                Eigen::Vector3d n_v = W(0,1) * U.block<3,1>(0,1);

                return PluckerLine(n_n , n_v); 

            }
            PluckerLine transform(const Eigen::Matrix4d& T){
                Eigen::Matrix3d R = T.block<3,3>(0,0);
                Eigen::Vector3d t = T.block<3,1>(0,3);
                Eigen::Matrix<double , 6,  6> H ;
                H.setZero();
                H.block<3,3>(0,0) = R;
                H.block<3,3>(3,3) = R;
                H.block<3,3>(0,3) = skew(t) * R;
                Vector6d L;
                L.segment<3>(0) = _n;
                L.segment<3>(3) = _v;

                Vector6d L1;
                L1 = H * L;
                //std::cout << "H : " << std::endl << H << std::endl << " L: " << L << std::endl << " L1 : " << L1 <<std::endl;
                return PluckerLine(L1.segment<3>(0) , L1.segment<3>(3));
            }
            inline double operator [](int i) const{
                assert(i < 6);
                if(i < 3){
                    return _n[i];
                }
                return _v[i-3];
            }

            inline Vector6d toVector() const{
                Vector6d v;
                v[0] = _n(0);
                v[1] = _n(1);
                v[2] = _n(2);
                v[3] = _v(0);
                v[4] = _v(1);
                v[5] = _v(2);
                return v;                                                
            }

            Vector4d log() const{
                //reference:
                //Building a 3-D Line-Based Map Using Stereo SLAM, TOR 2015
                Vector3d nn = _n / _n.norm();
                Vector3d vv = _v / _v.norm();
                Vector3d nv = _n.cross(_v);
                nv = nv / nv.norm();
                Eigen::Matrix3d R;
                R.block<3,1>(0,0) = nn;
                R.block<3,1>(0,1) = vv;
                R.block<3,1>(0,2) = nv;
                Eigen::Quaterniond qR(R);

                double n = qR.vec().norm();
                double w = qR.w();
                double squared_w = w*w;

                double two_atan_nbyw_by_n;
                if(n < SMALL_EPS_LINE){
                    assert(fabs(w) > SMALL_EPS_LINE);
                    two_atan_nbyw_by_n = 2. / w - 2. * (n*n)/(w*squared_w);
                }
                else{
                    if(fabs(w) < SMALL_EPS_LINE){
                        if(w > 0){
                            two_atan_nbyw_by_n = M_PI / n;
                        }
                        else{
                            two_atan_nbyw_by_n = -M_PI / n;
                        }

                    }
                    two_atan_nbyw_by_n = 2 * atan(n/w)/n;
                }

                Eigen::Vector3d thelta =  two_atan_nbyw_by_n*qR.vec();

                Eigen::Vector2d sigma;
                sigma << _n.norm() , _v.norm();
                sigma = sigma / sigma.norm();
                double angle = acos(sigma(0));

                Vector4d rec;
                rec << thelta(0) , thelta(1) , thelta(2) , angle;
                return rec;
            }

            static PluckerLine exp(const Vector4d& update){
                Eigen::Vector3d omega;// = update.segment<3>(0);
                for(int i = 0 ; i < 3 ;i ++){
                    omega[i] = update[i];
                }
                double angle = update[3];

                double theta = omega.norm();
                Matrix3d Omega = skew(omega);

                Matrix3d R;

                if(theta < 0.00001){
                    R = (Matrix3d::Identity() + Omega + Omega * Omega);
                }
                else{
                    R = (Matrix3d::Identity() + sin(theta)/theta * Omega
                        + (1.0 - cos(theta))/(theta*theta)*Omega*Omega);
                }
                double w1 = cos(angle);
                double w2 = sin(angle);

                Vector3d nn = w1 * R.block<3,1>(0,0);
                Vector3d vv = w2 * R.block<3,1>(0,1);

                return PluckerLine(nn , vv);


            }
            PluckerLine update(const Vector4d& update){
                Eigen::Vector3d omega;// = update.segment<3>(0);
                for(int i = 0 ; i < 3 ;i ++){
                    omega[i] = update[i];
                }
                double angle = update[3];

                double theta = omega.norm();
                Matrix3d Omega = skew(omega);

                Matrix3d R;

                if(theta < 0.00001){
                    R = (Matrix3d::Identity() + Omega + Omega * Omega);
                }
                else{
                    R = (Matrix3d::Identity() + sin(theta)/theta * Omega
                        + (1.0 - cos(theta))/(theta*theta)*Omega*Omega);
                }

                Matrix3d U = Umatrix();
                Vector2d W;
                W << _n.norm() , _v.norm();
                W = W / W.norm();
                Matrix2d W_angle;
                W_angle << cos(angle) , -sin(angle),
                            sin(angle) , cos(angle);
                
                W = W_angle * W;
                double w1 = W.x();
                double w2 = W.y();

                U = R * U;
                //U = U * R;

                return PluckerLine(w1 * U.block<3,1>(0,0) , w2 * U.block<3,1>(0,1));
            }
            Matrix3d Umatrix(){
                Matrix3d U;
                U.setZero();
                U.block<3,1>(0,0) = _n / _n.norm();
                U.block<3,1>(0,1) = _v / _v.norm();
                U.block<3,1>(0,2) = _n.cross(_v) / (_n.cross(_v).norm());
                return U;
            }
            Matrix<double , 3, 2> Wmatrix(){
                Matrix<double , 3, 2> W;
                W.setZero();
                W(0,0) = _n.norm();
                W(1,1) = _v.norm();
                double len = sqrt(W(0,0)*W(0,0)+W(1,1)*W(1,1));
                W = W / len;

                return W;
            }
    };

    inline std::ostream& operator << (std::ostream& out_str , const PluckerLine& pl){
        Eigen::Vector3d nn = pl.normal();
        Vector3d vv = pl.direction();
        out_str << nn.transpose() << " " << vv.transpose() << std::endl;
        return out_str;
    }
}

#endif