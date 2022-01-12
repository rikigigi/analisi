#ifndef TRICLINIC_H
#define TRICLINIC_H

#include "eigen_include.h"
template <class T>
class TriclinicLammpsCell {
public:
    using MatrixT = T;

    TriclinicLammpsCell (T * cell_) : cell{cell_} {
        if (cell(0,1)==0 &&
                cell(0,2)==0 &&
                cell(1,0)==0 &&
                cell(1,2)==0 &&
                cell(2,0)==0 &&
                cell(2,1)==0) {
            q = Eigen::Matrix<T,3,3>::Identity();
            r = cell;
            is_diagonal=true;
            return;
        }
        is_diagonal=false;
        auto qr = cell.transpose().householderQr();
        q = qr.householderQ();
        r = qr.matrixQR().template triangularView<Eigen::Upper>().transpose(); // lammps_tilt
        reflect();
    }
    void reflect(){
        for (int i=0;i<3;++i){
            if (r(i,i)<0) {
                Eigen::Matrix<T,3,3> m = Eigen::Matrix<T,3,3>::Identity();
                m(i,i)=-1;
                r = r*m;
                q = q*m;
            }
        }
    }
    void set_lammps_cell(T * cel,bool triclinic=true) const{
        //xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz
        cel[0]=0;
        cel[1]=r(0,0);
        cel[2]=0;
        cel[3]=r(1,1);
        cel[4]=0;
        cel[5]=r(2,2);
        if (triclinic){
            cel[6]=r(1,0);
            cel[7]=r(2,0);
            cel[8]=r(2,1);
        }
    }
    bool is_same_cell(const T* cel) const {
        for (int i =0; i<3;++i)
          for (int j =0; j<3;++j)
            if (cel[i*3+j] !=cell(i,j))
                return false;


        return true;
    }
    void rotate_vec(T* vec_) const{
        Eigen::Map<Eigen::Matrix<T,1,3>> vec(vec_);
        vec=vec*q;
    }
    T test() const {
        Eigen::Matrix<T,3,3> qr;
        qr = q*r.transpose();
        std::cerr << q << std::endl << r << std::endl<<cell << std::endl<<qr<<std::endl;
        return (qr-cell.transpose()).cwiseAbs2().sum();
    }
    bool isDiagonal() const {return is_diagonal;}
    const T * getCell() const {return cell.data();}
private:
    Eigen::Map<Eigen::Matrix<T,3,3>> cell;
    Eigen::Matrix<T,3,3> q,r;
    bool is_diagonal;
};

#endif // TRICLINIC_H
