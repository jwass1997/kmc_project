#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

class LiteLaplaceCircle {
public:
    // Construct a polar grid on a disk of radius R (same units you’ll use later).
    // Nr: radial nodes (>=2), Nt: angular nodes (>=8).
    // useSOR/omega: relaxation (1 < omega < 2).
    explicit LiteLaplaceCircle(double R, int Nr, int Nt,
                               bool useSOR = true, double omega = 1.85);

    // Reset the field to a constant.
    void clearSolution(double value = 0.0);

    // Define an electrode arc on the outer boundary by angle [radians].
    // begin/end may be in any order; they’re wrapped to [0, 2π).
    void setElectrode(double voltage, double begin, double end);

    // Convenience: define electrode by center [deg] & linear width (same units as R).
    // Example: R=150 nm, width=60 nm -> width_rad = 60/150.
    void setElectrodeLinearWidth(double voltage, double center_deg, double width_linear);

    // Convenience: define electrode by center [deg] & angular width [deg].
    void setElectrodeDeg(double voltage, double center_deg, double width_deg);

    // Update voltage of an existing electrode (by creation order).
    void updateElectrodeVoltage(int electrodeIndex, double voltage);

    // Solve Laplace with SOR/Jacobi until max change < tol or max_iters reached.
    int run(int max_iters = 10000, double tol = 1e-9);

    // Interpolate φ at Cartesian (x,y).
    double getPotential(double x, double y) const;

    // Accessors.
    int Nr() const { return Nr_; }
    int Ntheta() const { return Nt_; }
    double radius() const { return R_; }
    double dr() const { return dr_; }
    double dtheta() const { return dth_; }
    const std::vector<double>& rawGrid() const { return phi_; }

    std::vector<double> gridCopy() const { return phi_; }

private:
    // Geometry
    double R_;
    int Nr_, Nt_;
    double dr_, dth_;

    // Iteration
    bool useSOR_;
    double omega_;

    // Field (stored on [i=0..Nr-1, j=0..Nt-1])
    std::vector<double> phi_;

    // Outer boundary bookkeeping
    std::vector<bool>   outer_dirichlet_mask_;   // which theta cells are Dirichlet
    std::vector<double> outer_dirichlet_value_;  // their voltages
    std::vector<std::vector<int>> electrodes_;   // per-electrode θ indices (for updates)

    // --- helpers ---
    inline int idx(int i, int j) const { return i * Nt_ + ((j % Nt_) + Nt_) % Nt_; }
    inline int jp(int j) const { return (j + 1) % Nt_; }
    inline int jm(int j) const { return (j - 1 + Nt_) % Nt_; }
    inline double& at(int i, int j) { return phi_[idx(i, j)]; }
    inline const double& atc(int i, int j) const { return phi_[idx(i, j)]; }

    static double wrapAngle(double a);
    int thetaToIndex(double theta) const;

    // Apply BCs on the outer ring:
    // Dirichlet on electrode θ-cells; Neumann (∂φ/∂r = 0) elsewhere.
    void applyOuterBoundaryToSolution();
    void enforceOuterBC();

    // Internal “set electrode” that records only the indices created by this call.
    void setElectrodeRecordIndices(double voltage, double begin_wrapped, double end_wrapped);
};
