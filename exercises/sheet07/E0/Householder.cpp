#include <iostream>
#include <string>
#include "Householder.h"

// Method implementations
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Householder::Householder(int N) : N(N), M(Eigen::MatrixXd(N, N)) {} // Default constructor

// Public methods
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Householder::generateMatrix() {
    // Loop through rows
    for (int k = 0; k < N; k++) {
        // Loop through columns
        for (int l = 0; l < N; l++) {
            M(k, l) = k + l + (k == l? k : 0);
        }
    }
    sample("Sample of generated matrix:");

}

void Householder::sample(std::string text) {
    std::cout << text + "\n";
    std::cout << M.block(0,0,10,10) << "\n";  // Print the top-left 10x10 block of the matrix
}

void Householder::applyTransformation() {
    for (int i = 0; i < N - 2; i++) {
        calculate_A(i);
    }
    sample("\nSample of the tridiagonalized matrix:");
    std::cout << "\nRounded sample of tridiagonal matrix: " << M.block(0,0,10,10).array().round() << "\n";
}

Eigen::MatrixXd Householder::getMatrix() {
    return M;
}

// Private methods
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Eigen::VectorXd Householder::calculate_v(int i) {
    Eigen::VectorXd v = M.block(i+1, i, N-1-i, 1);
    return v;
}

Eigen::VectorXd Householder::calculate_u(int i) {
    Eigen::VectorXd v = calculate_v(i);
    double k = v.norm();
    // Match the sign of k with a21 
    k = abs(k);
    k = v(1) > 0? k : -k;

    Eigen::VectorXd u = v;
    u(0) = u(0) - k;
    u.normalize();
    return u;
}

Eigen::MatrixXd Householder::calculate_S(int i) {
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(N-i-1, N-i-1);
    Eigen::VectorXd u = calculate_u(i);
    Eigen::MatrixXd S = identity - 2 * (u * u.transpose());
    return S;
}

Eigen::MatrixXd Householder::calculate_P(int i) {
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(N, N);
    Eigen::MatrixXd S = calculate_S(i);
    P.block(0, 0, i+1, i+1) = Eigen::MatrixXd::Identity(i+1, i+1);
    P.block(i+1, i+1, N-1-i, N-1-i) = S;

    return P;
}

void Householder::calculate_A(int i) {
    Eigen::MatrixXd P = calculate_P(i);
    Eigen::MatrixXd A = P.transpose() * M * P;
    this->M = A;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::