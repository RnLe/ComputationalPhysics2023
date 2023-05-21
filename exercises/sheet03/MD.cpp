#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <fstream>
#include <iomanip>
#include <streambuf>
#include <omp.h>
#include <chrono>

using namespace std;
using namespace Eigen;

class LoggingStreamBuffer : public streambuf {
public:
    LoggingStreamBuffer(streambuf* buf, bool enabled)
        : m_buffer(buf), m_enabled(enabled) {}

protected:
    virtual int_type overflow(int_type ch) {
        if (m_enabled) {
            return m_buffer->sputc(ch);
        }
        return ch;
    }

private:
    streambuf* m_buffer;
    bool m_enabled;
};

// ================================ Potential-class ================================================

// Virtual class from which concrete potentials can be inherited
// (Here only Lennard-Jones necessary, but so you can quickly implement other potentials).
class Potential
{
    public:
        virtual double      V               ( double r2 ) const = 0;  // Virtual function
        virtual Vector2d    F               ( Vector2d r ) const = 0; // Virtual function
};

class PotentialLJ: public Potential
{
    public:
        double              V               ( double r2 ) const;  // Overwrites virtual function
        Vector2d            F               ( Vector2d r ) const; // Overwrites virtual function
};

// For the potential, the square of the vector length is sufficient, which saves a root calculation.
double PotentialLJ::V ( double r2 ) const
{
    return 4. * (pow(r2, -6) - pow(r2, -3));
}

Vector2d PotentialLJ::F ( Vector2d r ) const
{
    return (24. * (2. * pow(r.squaredNorm(), -7) - pow(r.squaredNorm(), -4) )) * r;
}

// ------------------------------ End of Potential-class -------------------------------------------
// ================================ Thermostat class ===============================================

// Virtual class from which concrete thermostats can be inherited
class Thermostat
{
    public:
        virtual void        rescale         ( vector<Vector2d>& v, double T ) const = 0;
};

// No thermostat
class NoThermostat: public Thermostat
{
    public:
        void                rescale         ( vector<Vector2d>& v, double T ) const {} // does nothing
};

// Isokinetic thermostat for task d)
class IsokinThermostat: public Thermostat
{
    public:
        void                rescale         ( vector<Vector2d>& v, double T ) const;
};

void IsokinThermostat::rescale( vector<Vector2d>& v, double T ) const
{
    double T_current = 0.0;
    double N = v.size();

    // Calculate the current temperature
    for(const auto& vel : v) {
        T_current += 1. / 2. * vel.squaredNorm();
    }
    T_current /= (N - 1.);

    // Calculate the scaling factor
    double scale_factor = sqrt(T / T_current);

    // Rescale the velocities
    for(auto& vel : v) {
        vel *= scale_factor;
    }
}

// ------------------------------ End of Thermostat class --------------------------------------
// ================================ Data-Structs ===============================================

// Data set for time resolved data
// (structs basically the same as class, but all members are public by default)
struct Dataset
{
    double                  t, T, Ekin, Epot;
    Vector2d                vS;
    vector<Vector2d>        r;
};

// Return data of the MD simulation
// Data data(n) constructor; reserves memory and fills pair correlation function with 0s
struct Data
{
    vector<Dataset> datasets; // Time-resolved datasets.
    vector<double> rBin, g;   // Averaged pair correlation function
    vector<Vector2d> r;       // snapshot of the final position
                              // For task e) it may be useful to use r
                              // in the time-resolved datasets instead

                            Data            ( uint n, uint numBins, double binSize );
    void                    save            ( const string& filenameSets,
                                              const string& filenameG,
                                              const string& filenameR ) const;
};

Data::Data( uint n, uint numBins, double binSize ):
    datasets( n ),  // Initializer list, because it calls constructors of the members
    rBin( numBins ),
    g( numBins, 0. ),
    r( 0 )
{
}

void Data::save ( const string& filenameSets, const string& filenameG, const string& filenameR ) const
{
    clog << "Writing files\n";
    // Save time-resolved datasets
    ofstream fileSets(filenameSets);
    if (!fileSets.is_open()) {
        cerr << "Error: Unable to open file " << filenameSets << "\n";
        return;
    }

    clog << "Writing " << filenameSets << "..\n";
    fileSets << "t" << ", " << "T" << ", " << "Ekin" << ", " << "Epot" << ", " << "vS.x" << ", " << "vS.y" << "\n";
    for (const auto &dataset : datasets) {
        fileSets << dataset.t << ", " << dataset.T << ", " << dataset.Ekin << ", " << dataset.Epot << ", " << dataset.vS[0] << ", " << dataset.vS[1] << "\n";
    }
    fileSets.close();

    // Save averaged pair correlation function
    ofstream fileG(filenameG);
    if (!fileG.is_open()) {
        cerr << "Error: Unable to open file " << filenameG << "\n";
        return;
    }

    clog << "Writing " << filenameG << "..\n";
    fileG << "rBin" << ", " << "g" << "\n";
    for (size_t i = 0; i < rBin.size(); ++i) {
        fileG << rBin[i] << ", " << g[i] << "\n";
    }
    fileG.close();

    // Save the final positions
    ofstream fileR(filenameR);
    if (!fileR.is_open()) {
        cerr << "Error: Unable to open file " << filenameR << "\n";
        return;
    }

    clog << "Writing " << filenameR << "..\n";
    unsigned int particleCount = datasets[0].r.size();
        
    // Create headers
    for (unsigned int i = 0; i < particleCount; i++) {
        fileR << "x" << i+1 << ", " << "y" << i+1;
        if (i != particleCount - 1) {
            fileR << ", ";
        }
    }
    fileR << "\n";

    // Write data
    for (const auto &dataset : datasets) {
        for (unsigned int i = 0; i < particleCount; i++) {
            fileR << dataset.r[i][0] << ", " << dataset.r[i][1];
            if (i != particleCount - 1) {
                fileR << ", ";
            }
        }
        fileR << "\n";
    }
    fileR.close();
    clog << "\033[1;33mWriting files finished\n\033[0m";

}

// ------------------------------ End of Data-Structs --------------------------------------
// ================================ MD-Class ===============================================

class MD
{
    public:
                            MD              ( double L, uint N, uint particlesPerRow, double T,
                                            Potential& potential, Thermostat& thermostat,
                                            uint numBins = 1000 );

        void                equilibrate     ( const double dt, const unsigned int n );
        Data                measure         ( const double dt, const unsigned int n );

    private:
        vector<Vector2d>    r, v;
        double              L;
        uint                N;
        double              T;
        Potential&          potential;
        Thermostat&         thermostat;
        double              t = 0.;
        uint                particlesPerRow;

        uint                numBins;
        double              binSize;

        // Particles are moved in box [0,L]x[0,L].
        void                centerParticles ();

        // Calculations of important measured variables
        double              calcT           () const;
        double              calcEkin        () const;
        double              calcEpot        () const;
        Vector2d            calcvS          () const;
        Dataset             calcDataset     () const;

        // Calculation of the acceleration
        // To avoid redundant calculations, it may be useful to update the histogram
        // when calculating the accelerations, so it is passed here as a reference.
        vector<Vector2d>    calcAcc         ( vector<double>& hist ) const;

        // Calculation of the distance vector between particle r[i] and closest mirror particle of r[j].
        Vector2d            calcDistanceVec ( uint i, uint j ) const;
};

// Initialization of the system via constructor
MD::MD( double L, uint N, uint particlesPerRow, double T,
        Potential& potential, Thermostat& thermostat,
        uint numBins ):
    L(L),
    N(N),
    T(T),
    potential( potential ),
    thermostat( thermostat ),
    particlesPerRow( particlesPerRow ),
    numBins( numBins ),
    binSize( L / numBins )
{
    centerParticles();
    {
        IsokinThermostat tempThermo;
        tempThermo.rescale(v, T);
    }
}

// Integration without data acquisition for pure equilibration
void MD::equilibrate ( const double dt, const unsigned int n )
{
    clog << "Equilibration started..\n";
    // For verlet's algorithm, we need to define the first values for r once.
    vector<double> hist(numBins, 0.);
    vector<Vector2d> prevV(N);
    vector<Vector2d> prevA(N); 

    for (uint step = 0; step < n; step++) {
        // 1. Update positions and velocities using the Verlet algorithm
        // 2. Apply the thermostat to rescale velocities if needed

        // 1. Update positions and velocities
        vector<double> hist(numBins, 0.);
        vector<Vector2d> a = calcAcc(hist);
        prevA = a;
        for (uint i = 0; i < N; ++i) {
            r[i] = r[i] + v[i] * dt +  1./2. * prevA[i] * dt * dt;

            // Check if particle has moved out of bounds and correct if necessary
            for (int j = 0; j < 2; ++j) {  // loop over x and y coordinates
                if (r[i][j] < 0.) {
                    r[i][j] += L;
                } else if (r[i][j] >= L) {
                    r[i][j] -= L;
                }
            }
        }
        a = calcAcc(hist);
        for (uint i = 0; i < N; i++) { 
            v[i] = v[i] + 1./2. * (a[i] + prevA[i]) * dt;
        }
        prevA = a;
        prevV = v;

        // 2. Apply thermostat
        thermostat.rescale(v, T);

        t += dt;

        // Progress bar
        // Update progress
        if (step % (n / 100) == 0) { // Update every 1 percent
            clog << "\033[1;32mProgress: " << step / (n / 100) << "%\033[0m\r";
            clog.flush();
        }
    }

    clog << "\033[1;33mEquilibration finished\n\033[0m";
}

Data MD::measure ( const double dt, const unsigned int n )
{
    clog << "Measuring data\n";
    // #pragma omp parallel
    Data data(n, numBins, binSize);

    // For verlet's algorithm, we need to define the first values for r once.
    vector<Vector2d> prevV(N);
    vector<Vector2d> prevA(N);
    vector<double> hist(numBins, 0.);
    vector<Vector2d> a = calcAcc(hist);
    prevA = a;

    for (uint step = 0; step < n; step++) {
        // 1. Calculate and store the Dataset with the current state
        // 2. Update positions and velocities using the Verlet algorithm
        // 3. Apply the thermostat to rescale velocities if needed
        // 4. Update the pair correlation function histogram

        // 1. Calculate and store the Dataset
        data.datasets[step] = calcDataset();
        vector<double> hist(numBins, 0.);

        // 2. Update positions and velocities (Velocity Verlet Algorithm)
        for (uint i = 0; i < N; ++i) {
            r[i] = r[i] + v[i] * dt +  1./2. * prevA[i] * dt * dt;

            // Check if particle has moved out of bounds and correct if necessary
            for (int j = 0; j < 2; ++j) {  // loop over x and y coordinates
                if (r[i][j] < 0.) {
                    r[i][j] += L;
                } else if (r[i][j] >= L) {
                    r[i][j] -= L;
                }
            }
        }
        a = calcAcc(hist);
        for (uint i = 0; i < N; i++) { 
            v[i] = v[i] + 1./2. * (a[i] + prevA[i]) * dt;
        }
        prevA = a;
        prevV = v;

        // 3. Apply thermostat
        thermostat.rescale(v, T);

        // 4. Update the pair correlation function histogram
        for (uint bin = 0; bin < numBins; ++bin) {
            data.g[bin] += hist[bin];
        }
        t += dt;

        // Progress bar
        // Update progress
        if (step % (n / 100) == 0) { // Update every 1 percent
            clog << "\033[1;32mProgress: " << step / (n / 100) << "%\033[0m\r";
            clog.flush();
        }
    }

    // clog << "Normalize the pair correlation function histogram\n";
    // Normalize the pair correlation function histogram
    double norm = 1.0 / (M_PI * binSize * sqrt(N) * n * 2.5);
    for (uint bin = 0; bin < numBins; bin++) {
        double rInner = bin * binSize;
        double rOuter = (bin + 1) * binSize;
        double area = M_PI * (rOuter * rOuter - rInner * rInner);
        data.g[bin] *= norm / area;
        data.rBin[bin] = (rOuter + rInner) * 0.5;
    }

    // Save the final positions
    data.r = r;

    clog << "\033[1;33mMeasuring finished\n\033[0m";
    return data;
}

void MD::centerParticles()
{
    clog << "Centering particles\n";
    // The positions are supposed to be equidistant, while the velocities are random.
    // Initialize the positions on a square grid

    // Positions
    double equidistance = L / static_cast<double>(particlesPerRow);         // This is just 2..
    // clog << "Equidistance: " << equidistance << "\n";
    for (uint i = 0; i < particlesPerRow; i++) {
        for (uint k = 0; k < particlesPerRow; k++) {
            r.push_back(Vector2d(i * equidistance + 1., k * equidistance + 1.));
        } 
    }
    
    // Velocities
    // Keep track of the total velocity
    Vector2d totalVelocity = Vector2d::Zero();

    mt19937 rnd;
    uniform_real_distribution<double> dist(-1, 1);

    for (uint i = 0; i < N; i++) {
        v.push_back(Vector2d(dist(rnd), dist(rnd)));
        totalVelocity += v[i];
    }
    // clog << "Total velocity before normalization: " << totalVelocity << "\n";
    // Average velocity of all particles
    // We keep track of this so we can subtract this value from each particles velocity
    totalVelocity /= static_cast<double>(N);
    // clog << "Average velocity (after normalization): " << totalVelocity << "\n";

    for (uint i = 0; i < N; i++) {
        v[i] -= totalVelocity;
    }
}

double MD::calcT() const
{
    double Ekin = calcEkin();
    return Ekin / (N - 1);    // Derived from virial theorem. 2D system has only two degrees of freedom.
}

double MD::calcEkin() const
{   
    // We remember that m = 1
    double totalEkin = 0.0;
    for (uint i = 0; i < N; i++) {
        totalEkin += 0.5 * v[i].squaredNorm();     // Dot product squared
    }
    return totalEkin;
}

double MD::calcEpot() const
{
    double totalEpot = 0.0;
    for (uint i = 0; i < N; i++) {
        for (uint k = i + 1; k < N; k++) {
            totalEpot += potential.V(calcDistanceVec(i,k).squaredNorm()) - potential.V(pow(0.5 * L, 2));       // Squared distance, which is respected in the potential function.
        }
    }
    return totalEpot;
}

Vector2d MD::calcvS() const
{
    Vector2d totalVelocities = Vector2d::Zero();
    for (uint i = 0; i < N; i++) {
        totalVelocities += v[i];
    }
    return totalVelocities / static_cast<double>(N);
}

Dataset MD::calcDataset() const
{
    Dataset dataset;
    dataset.t = t;
    dataset.T = calcT();
    dataset.Ekin = calcEkin();
    dataset.Epot = calcEpot();
    dataset.vS = calcvS();
    dataset.r = r;
    return dataset;
}

Vector2d MD::calcDistanceVec( uint i, uint k ) const
{
    // Minimum-image convention (PCB particle bookkeeping)
    // Vector pointing to r[k]
    Vector2d r_ik = r[k] - r[i];
    r_ik(0) -= L * floor((r_ik(0) / L) + 0.5);
    r_ik(1) -= L * floor((r_ik(1) / L) + 0.5);

    return r_ik;
}


vector<Vector2d> MD::calcAcc( vector<double>& hist ) const
{
    vector<Vector2d> acc(N, Vector2d::Zero());
    
    for (uint i = 0; i < N - 1; i++) {
        for (uint k = i + 1; k < N; k++) {
            Vector2d distanceVec = calcDistanceVec(i, k);
            double r2 = distanceVec.squaredNorm();

            // Check if the distance is within the cutoff radius
            if (r2 < (0.5 * L) * (0.5 * L)) {
                Vector2d force = potential.F(distanceVec);
                
                acc[i] -= force;
                acc[k] += force;

                // Update histogram for the pair correlation function g(r)
                uint bin = static_cast<uint>(sqrt(r2) / binSize);
                if (bin < numBins) {
                    hist[bin]++;
                }
            }
        }
    }

    return acc;
}

// ------------------------------ End of MD-class ------------------------------------------


int main(void)
{
    // Flags
    bool b = true;
    bool c = true;
    bool d = true;
    bool logging_enabled = true;

    // Timer
    auto start_time = chrono::high_resolution_clock::now();
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();

    // Create the custom stream buffer and connect it to clog
    LoggingStreamBuffer log_buffer(clog.rdbuf(), logging_enabled);
    clog.rdbuf(&log_buffer);
    clog << "Logging enabled\n";

    PotentialLJ      LJ;
    NoThermostat     noThermo;
    IsokinThermostat isoThermo;

    string n[3] = { "4", "8", "16" };
    uint partPerRow;                // Particles per row. I assume that the whole x-axis is meant by that.
    uint N;
    double L;
    int numBins                     = 512;

    if (b) {
        clog << "\033[1;34mb) Equilibration test\n\033[0m";
        // b) Equilibration test
        for (int i = 0; i < 3; i++) {
            start_time = chrono::high_resolution_clock::now();

            const double T          = 1.0;
            const double dt         = 0.01;
            const uint steps        = 10000;

            partPerRow                 = stoi(n[i]);
            N                          = stoi(n[i]) * stoi(n[i]);
            L                          = 2. * stoi(n[i]);

            MD md( L, N, partPerRow, T, LJ, noThermo, numBins );
            md.measure( dt, steps ).save( "b)set_n" + n[i] + ".dat", "b)g_n" + n[i] + ".dat", "b)r_n" + n[i] + ".dat" );
            end_time = chrono::high_resolution_clock::now();
            duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
            clog << "\033[1;34mTime needed for b), " << n[i] << " particles: " << duration / 1000.0 << " Seconds\n\n\033[0m";
        }
        clog << "\033[1;33mb) Finished\n\n\033[0m";
    }

    partPerRow                 = stoi(n[1]);
    N                          = stoi(n[1]) * stoi(n[1]);
    L                          = 2. * stoi(n[1]);
    string TstringVec[3]       = { "0.01", "1", "100" };

    if (c) {
        clog << "\033[1;34mc) Pair correlation function\n\033[0m";
        // c) Pair correlation function
        
        for ( auto& Tstring: TstringVec ) {
            start_time = chrono::high_resolution_clock::now();

            const double T          = stod(Tstring);
            const double dt         = (Tstring == "0.01") ? 0.001 : 0.001;     // For 0.01 K, use a smaller step size
            const uint equiSteps    = 100000;
            const uint steps        = 10000;

            MD md( L, N, partPerRow, T, LJ, noThermo, numBins );
            md.equilibrate( dt, equiSteps );
            md.measure( dt, steps ).save( "c)set" + Tstring + ".dat", "c)g" + Tstring + ".dat", "c)r" + Tstring + ".dat" );

            end_time = chrono::high_resolution_clock::now();
            duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
            clog << "\033[1;34mTime needed for c), " << Tstring << " kelvin: " << duration / 1000.0 << " Seconds\n\n\033[0m";
        }
        clog << "\033[1;33mc) Finished\n\n\033[0m";
    }    

    if (d) {
        clog << "\033[1;34md) Isokinetic thermostat\n\033[0m";
        // d) Thermostat
        for ( auto& Tstring: TstringVec ) {
            start_time = chrono::high_resolution_clock::now();

            const double T          = stod(Tstring);
            const double dt         = (Tstring == "0.01") ? 0.001 : 0.001;     // For 0.01 K, use a smaller step size
            const uint equiSteps    = 100000;
            const uint steps        = 10000;

            MD md( L, N, partPerRow, T, LJ, isoThermo, numBins );
            md.equilibrate( dt, equiSteps );
            md.measure( dt, steps ).save( "d)set" + Tstring + ".dat", "d)g" + Tstring + ".dat", "d)r" + Tstring + ".dat" );

            end_time = chrono::high_resolution_clock::now();
            duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
            clog << "\033[1;34mTime needed for d), " << Tstring << " kelvin: " << duration / 1000.0 << " Seconds\n\n\033[0m";
        }
    }

    return 0;
}
