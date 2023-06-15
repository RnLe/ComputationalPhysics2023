#include <Eigen/Dense>
#include <iostream>
#include <array>
#include <string>

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

        Eigen::VectorXd     calculate_v         (int i);
        Eigen::VectorXd     calculate_u         (int i);
        Eigen::MatrixXd     calculate_S         (int i);
        Eigen::MatrixXd     calculate_P         (int i);
        void                calculate_A         (int i);
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// Method implementations
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Constructors
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Householder::Householder(int N) : N(N), M(Eigen::MatrixXd(N, N)) {} // Default constructor

// Public methods
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Householder::generateMatrix() {
    // Loop through rows
    for (int k = 0; k < N; k++) {
        // Loop through columns
        for (int l = 0; l < N; l++) {
            M(k, l) = k + l + (k == l? k: 0);
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
    // M = M.array().round();
    M = M.unaryExpr([](double v) { return (abs(v) < 1e-8) ? 0.0 : v; });
    std::cout << "\nRounded sample of tridiagonal matrix:\n" << M.array().round() << "\n";
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
    k = v(0) >= 0? k : -k;

    Eigen::VectorXd u = -v;
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

    std::cout << "\n\nP " << i+1 << ":\n" << P;
    return P;
}

void Householder::calculate_A(int i) {
    Eigen::MatrixXd P = calculate_P(i);
    Eigen::MatrixXd A = P.transpose() * M * P;
    M = A;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

class QR {
    public:                 QR                  (Eigen::MatrixXd M);
        Eigen::MatrixXd     getMatrix           ();
        void                sample              (std::string text);
        void                applyTransformation (bool logging);
        void                checkIntegrity      ();

    private:
        Eigen::MatrixXd     M;
        bool                flipped = false;

        Eigen::MatrixXd     calculate_Q         ();
        Eigen::MatrixXd     calculate_P         (int i);
        void                flipMatrix          ();
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// Method implementations
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Constructors
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

QR::QR(Eigen::MatrixXd M) : M(M) {}

// Public methods
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Eigen::MatrixXd QR::getMatrix() {
    return M;
}

void QR::sample(std::string text) {
    std::cout << text + "\n";
    std::cout << M.block(0,0,10,10) << "\n";  // Print the top-left 10x10 block of the matrix
}

void QR::applyTransformation(bool logging) {
    Eigen::MatrixXd Q = calculate_Q();
    M = M * Q;
    if (logging) sample("\nTransformed matrix:\n");
}

void QR::checkIntegrity() {
    // Flags
    bool first = false, last = false, between = false;
    std::cout << "\n";

    // Check whether the matrix has zeroes on its diagonal
    // If it does, things get problematic. If only the first element is zero, the matrix can be flipped.
    // Effectively, only a zero at the last diagonal element can be handled.
    for (int i = 0; i < M.rows(); i++) {
        if (i == 0 && M(i, i) == 0) {
            std::cout << "First diagonal element is zero.\n";
            first = true;
        } else if(i == (M.rows()-1) && M(i, i) == 0) {
            std::cout << "Last diagonal element is zero.\n";
            last = true;
        } else if (!between && M(i, i) == 0) {
            std::cout << "Diagonal element inbetween is zero.\n";
            between = true;
        }
    }

    if (between) {
        std::cout << "QR algorithm not applicable.\n";
    } else if (first && !last) {
        std::cout << "QR algorithm applicable, but matrix has to be flipped.\n";
        flipMatrix();
        flipped = true;
    } else if (!first && last) {
        std::cout << "QR algorithm applicable.\n";
    } else if (first && last) {
        std::cout << "QR algorithm not applicable.\n";
    } else {
        std::cout << "QR algorithm applicable.\n";
    }
}

// Private methods
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Eigen::MatrixXd QR::calculate_Q() {
    // To prevent a recursive approach or the necessity to store all N P-matrices, a temporary matrix is introduced
    Eigen::MatrixXd P_final = calculate_P(0);
    Eigen::MatrixXd P_temp = P_final;
    for (int i = 1; i <= M.rows()-2; i++) {
        M = P_temp * M;
        P_temp = calculate_P(i);
        P_final *= P_temp;
        // std::cout << "\n" << P_temp << "\n";
    }
    return P_final;
}

Eigen::MatrixXd QR::calculate_P(int i) {
    // i = 0 equals P12
    double a = M(i, i);
    double k = M(i + 1, i);
    double t = k / a;
    double c = 1 / sqrt(t * t + 1);
    double s = t * c;

    Eigen::MatrixXd P = Eigen::MatrixXd::Identity(M.rows(), M.cols());
    Eigen::MatrixXd P_small(2, 2);
    P_small(0, 0) = c;
    P_small(0, 1) = s;
    P_small(1, 0) = -s;
    P_small(1, 1) = c;
    // std::cout << "\nMatrix P small i=" + std::to_string(i) + ":\n" << P_small << "\n";
    P.block(i, i, 2, 2) = P_small;

    std::cout << "\nMatrix P with i=" + std::to_string(i) + ":\n" << P << "\n";

    return P;
}

void QR::flipMatrix() {
    M.reverseInPlace();
    std::cout << "\nReversed matrix:\n" << M << "\n";
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::





// Main
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main(void) {

    // Initialization
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    std::array<int, 1> N = {10};    // Add more dimensions
    std::array<Eigen::MatrixXd, 3> Matrices;

    std::cout << "\nHouseholder transformation\n\n";

    for (int i = 0; i < N.size(); i++) {
        Householder M(N[i]);
        M.generateMatrix();
        M.applyTransformation();
        Matrices[i] = M.getMatrix();
    }

    for (int i = 0; i < N.size(); i++) {
        QR TriangleMatrix(Matrices[i]);
        TriangleMatrix.checkIntegrity();

        // Start QR-algorithm
        double threshold = 1e-20;   // Set the threshold here
        int max_iterations = 10;     // Set the iteration limit here
        Eigen::MatrixXd M_old;

        for (int iteration = 0; iteration < max_iterations; ++iteration) {
            M_old = TriangleMatrix.getMatrix();
            TriangleMatrix.applyTransformation(false);

            
            if ((M_old - TriangleMatrix.getMatrix()).norm() == threshold) {
                std::cout << "QR algorithm converged after " << iteration << " iterations.\n";
                break;
            }
            if (iteration == max_iterations-1) std::cout << "Iteration limit reached. " << iteration + 1 << " iterations.\n";
        }
        std::cout << "Triangle matrix\n" << TriangleMatrix.getMatrix(); 
    }
    
    return 0;
}