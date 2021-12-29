#include "nevanlinna.h"
//LS: Generate file name
#include <stdio.h>


template <class T>
class Schur : precision_<T> {
private:
    using typename precision_<T>::nev_complex;
    using typename precision_<T>::nev_complex_vector;
    using typename precision_<T>::nev_complex_matrix;
    using typename precision_<T>::nev_complex_matrix_vector;
public:
    //check Nevanlinna/contractive interpolant existence condition
    Schur (std::string ifile, int imag_num, std::string ofile);
    //evaluation with 0 parametric function 
    void evaluation ();
private:
    int M; //number of Matsubara points
    imag_domain_data <T> imag; //theta values at Matsubara points (G -> NG -> theta)
    real_domain_data <T> real; //real frequency NG storage, at omega + i*eta
    nev_complex_vector phis; //phi_1 to phi_M
    nev_complex_matrix_vector abcds; //intermediate {a, b, c, d}s used to calculate phis
    //memoize intermediate abcds and calculate phis by iteration
    void core ();
};


template <class T>
Schur<T>::Schur (std::string ifile, int imag_num, std::string ofile) : imag(ifile, imag_num), real(ofile)  {
    M = imag_num;
    //fill the Pick matrix
    nev_complex_matrix Pick (M, M);
    nev_complex I {0., 1.};
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            nev_complex freq_i = (imag.freq()[i] - I) / (imag.freq()[i] + I);
            nev_complex freq_j = (imag.freq()[j] - I) / (imag.freq()[j] + I);
            nev_complex one {1., 0.};
            nev_complex nom = one - imag.val()[i] * std::conj(imag.val()[j]);
            nev_complex den = one - freq_i * std::conj(freq_j);
            Pick(i, j) = nom / den;
        }
    }
    //check the positive semi-definiteness of the Pick matrix using Cholesky decomposition
    Eigen::LLT<nev_complex_matrix> lltOfPick(Pick + nev_complex_matrix::Identity(M, M) * 1e-250);
    if(lltOfPick.info() == Eigen::NumericalIssue) 
        std::cerr << "Pick matrix is non positive semi-definite matrix in Schur method." << std::endl;
    else std::cerr << "Pick matrix is positive semi-definite." << std::endl;
}


template <class T>
void Schur<T>::core() {
    phis.resize(M);
    abcds.resize(M);
    
    //LS: creat data files
    int pc=64;
    char filename[128], sparamfile[128], abcdfile[128], hiwnfile[128];
    std::ofstream lsfile;
    lsfile.precision(pc);
    sprintf(sparamfile, "cpp_code/F%i/sparam_gaussian_1_N_%i.txt",pc, M);
    lsfile.open(sparamfile);
    lsfile.close();
    sprintf(abcdfile, "cpp_code/F%i/abcd_gaussian_1_N_%i.txt",pc, M);
    lsfile.open(abcdfile);
    lsfile.close();
    sprintf(hiwnfile, "cpp_code/F%i/hiwn_gaussian_1_N_%i_F64.txt", pc, M);
    lsfile.open(hiwnfile);
    lsfile.close();

    //LS: print wn and h(-g(iwn))
    for (int j = 0; j < M; j++){
        lsfile.open(hiwnfile, std::ios::app);
        lsfile << imag.freq()[j].imag()<< " " << \
            imag.val()[j].real()<< " " << imag.val()[j].imag()<<" "<< std::endl;
        lsfile.close();
    }
    
    phis[0] = imag.val()[0];
    for (int k = 0; k < M; k++) abcds[k] = nev_complex_matrix::Identity(2, 2);
    for (int j = 0; j < M - 1; j++) {
        //LS: creat theta file
        sprintf(filename, "cpp_code/F%i/theta_%i_gaussian_1_N_%i.txt",pc,j+1,M);
        lsfile.open(filename);
        lsfile.close();

        for (int k = j; k < M; k++) {
            nev_complex_matrix prod(2, 2);
            prod(0, 0) = (imag.freq()[k] - imag.freq()[j]) / (imag.freq()[k] - std::conj(imag.freq()[j]));
            prod(0, 1) = phis[j];
            prod(1, 0) = std::conj(phis[j])*
                        ((imag.freq()[k] - imag.freq()[j]) / (imag.freq()[k] - std::conj(imag.freq()[j])));
            prod(1, 1) = nev_complex{1., 0.};
            abcds[k] *= prod;

            //LS : print theta
            lsfile.open(filename, std::ios::app);
            lsfile << imag.freq()[k].imag()<< " " << \
                prod(0,0).real()<< " " << prod(0,0).imag()<<" "<< \
                prod(0,1).real()<< " " << prod(0,1).imag()<<" "<< \
                prod(1,0).real()<< " " << prod(1,0).imag()<<" "<<  \
                prod(1,1).real()<< " " << prod(1,1).imag()<<" "<< std::endl;
            lsfile.close();
        }
        phis[j + 1] = (- abcds[j + 1](1, 1) * imag.val()[j + 1] + abcds[j + 1](0, 1)) /
                        (abcds[j + 1](1, 0) * imag.val()[j + 1] - abcds[j + 1](0, 0)); 

        //LS: print phis
        lsfile.open(sparamfile, std::ios::app);
        lsfile << imag.freq()[j].imag()<<" "<<phis[j].real() << " " <<phis[j].imag()<< std::endl;
        if(j==M-2)
        {
        lsfile << imag.freq()[j+1].imag() <<" "<<phis[j+1].real() << " " <<phis[j+1].imag()<< std::endl;
        } 
        lsfile.close();

        //LS:print abcd
        lsfile.open(abcdfile, std::ios::app);
        lsfile << imag.freq()[j].imag()<< " " << \
            abcds[j+1](0,0).real()<< " " << abcds[j+1](0,0).imag()<<" "<< \
            abcds[j+1](0,1).real()<< " " << abcds[j+1](0,1).imag()<<" "<< \
            abcds[j+1](1,0).real()<< " " << abcds[j+1](1,0).imag()<<" "<<  \
            abcds[j+1](1,1).real()<< " " << abcds[j+1](1,1).imag()<<" "<< std::endl;
        lsfile.close();

        /* std::cout.precision(128); 
        std::cout << imag.freq()[j].imag()<<" "<<phis[j].real() << " " <<phis[j].imag()<< std::endl;
        if(j==M-2)
        {
        std::cout << imag.freq()[j+1].imag() <<" "<<phis[j+1].real() << " " <<phis[j+1].imag()<< std::endl;
        }  */ 
}
}


template <class T>
void Schur<T>::evaluation () {
    core();
    nev_complex I {0., 1.};
    nev_complex One {1., 0.};
    for (int i = 0; i < real.N_real(); i++) {
        nev_complex_matrix result = nev_complex_matrix::Identity(2, 2);
        nev_complex z = real.freq()[i];
        for (int j = 0; j < M; j++) {
            nev_complex_matrix prod(2, 2);
            prod(0, 0) = (z - imag.freq()[j]) / (z - std::conj(imag.freq()[j]));
            prod(0, 1) = phis[j];
            prod(1, 0) = std::conj(phis[j])*
                        ((z - imag.freq()[j]) / (z - std::conj(imag.freq()[j])));
            prod(1, 1) = nev_complex{1., 0.};
            result *= prod;
        }
        nev_complex param {0., 0.}; //theta_{M+1}, choose to be constant function 0 here
        nev_complex theta = (result(0, 0) * param + result(0, 1)) / (result(1, 0) * param + result(1, 1));
        //can output "real.freq(), a.real(), a.imag(), ..., d.imag()\n" into a file for optimization convenience
        real.val()[i] = I * (One + theta) / (One - theta); //inverse Mobius transform from theta to NG
    }
    real.write();
}
