#include <Eigen/Dense>
#include <string>

#pragma once

// Class declaration
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

class Householder {
    public:                 Householder         (int N);
        void                generateMatrix      ();
        void                sample              (std::string text);
        void                applyTransformation ();
        Eigen::MatrixXd     getMatrix           ();

    private:
        int                 N;
        Eigen::MatrixXd     M;

        Eigen::VectorXd     calculate_v(int i);
        Eigen::VectorXd     calculate_u(int i);
        Eigen::MatrixXd     calculate_S(int i);
        Eigen::MatrixXd     calculate_P(int i);
        void                calculate_A(int i);
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::