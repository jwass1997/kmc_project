#include "LiteLaplaceCircle.h"

LiteLaplaceCircle::LiteLaplaceCircle(double R, int Nr, int Nt,
                                     bool useSOR, double omega)
    : R_(R), Nr_(Nr), Nt_(Nt), useSOR_(useSOR), omega_(omega)
{
    if (R_ <= 0) throw std::invalid_argument("R must be > 0");
    if (Nr_ < 2) throw std::invalid_argument("Nr must be >= 2");
    if (Nt_ < 8) throw std::invalid_argument("Nt must be >= 8");
    if (useSOR_ && (omega_ <= 1.0 || omega_ >= 2.0))
        throw std::invalid_argument("omega must be in (1,2)");

    dr_  = R_ / (Nr_ - 1);
    dth_ = 2.0 * M_PI / Nt_;

    phi_.assign(Nr_ * Nt_, 0.0);
    outer_dirichlet_mask_.assign(Nt_, false);
    outer_dirichlet_value_.assign(Nt_, 0.0);
}

void LiteLaplaceCircle::clearSolution(double value) {
    std::fill(phi_.begin(), phi_.end(), value);
}

double LiteLaplaceCircle::wrapAngle(double a) {
    double t = std::fmod(a, 2.0 * M_PI);
    return (t < 0.0) ? t + 2.0 * M_PI : t;
}

int LiteLaplaceCircle::thetaToIndex(double theta) const {
    double t = wrapAngle(theta);
    int j = static_cast<int>(std::round(t / dth_)) % Nt_;
    if (j < 0) j += Nt_;
    return j;
}

void LiteLaplaceCircle::setElectrode(double voltage, double begin, double end) {
    begin = wrapAngle(begin);
    end   = wrapAngle(end);
    setElectrodeRecordIndices(voltage, begin, end);
    applyOuterBoundaryToSolution();
}

void LiteLaplaceCircle::setElectrodeLinearWidth(double voltage, double center_deg, double width_linear) {
    const double center_rad = center_deg * M_PI / 180.0;
    const double half_width_rad = 0.5 * (width_linear / R_); // arc-length -> angle
    setElectrode(voltage, center_rad - half_width_rad, center_rad + half_width_rad);
}

void LiteLaplaceCircle::setElectrodeDeg(double voltage, double center_deg, double width_deg) {
    const double center_rad = center_deg * M_PI / 180.0;
    const double half_width_rad = 0.5 * (width_deg * M_PI / 180.0);
    setElectrode(voltage, center_rad - half_width_rad, center_rad + half_width_rad);
}

void LiteLaplaceCircle::setElectrodeRecordIndices(double voltage, double begin_wrapped, double end_wrapped) {
    std::vector<int> this_arc;
    auto mark_index = [&](int j) {
        outer_dirichlet_mask_[j]  = true;
        outer_dirichlet_value_[j] = voltage;
        this_arc.push_back(j);
    };

    if (begin_wrapped <= end_wrapped) {
        int j0 = thetaToIndex(begin_wrapped);
        int j1 = thetaToIndex(end_wrapped);
        for (int j = j0; j <= j1; ++j) mark_index(j % Nt_);
    } else {
        int j0 = thetaToIndex(begin_wrapped);
        int j1 = thetaToIndex(end_wrapped);
        for (int j = j0; j < j0 + Nt_; ++j) {
            int jj = j % Nt_;
            mark_index(jj);
            if (jj == j1) break;
        }
    }
    electrodes_.push_back(std::move(this_arc));
}

void LiteLaplaceCircle::updateElectrodeVoltage(int electrodeIndex, double voltage) {
    if (electrodeIndex < 0 || electrodeIndex >= static_cast<int>(electrodes_.size()))
        throw std::out_of_range("Bad electrode index");
    for (int j : electrodes_[electrodeIndex]) {
        outer_dirichlet_mask_[j]  = true;
        outer_dirichlet_value_[j] = voltage;
    }
    applyOuterBoundaryToSolution();
}

void LiteLaplaceCircle::applyOuterBoundaryToSolution() {
    // Dirichlet on electrode cells, Neumann (∂φ/∂r=0) elsewhere.
    for (int j = 0; j < Nt_; ++j) {
        if (outer_dirichlet_mask_[j]) {
            at(Nr_ - 1, j) = outer_dirichlet_value_[j];
        } else {
            at(Nr_ - 1, j) = at(Nr_ - 2, j); // copy from interior ring
        }
    }
}
void LiteLaplaceCircle::enforceOuterBC() { applyOuterBoundaryToSolution(); }

/* int LiteLaplaceCircle::run(int max_iters, double tol) {
    for (int iter = 0; iter < max_iters; ++iter) {
        double max_change = 0.0;

        // Outer BC each sweep
        enforceOuterBC();

        // Center regularity: φ(0,θ) = mean(φ(1,θ))
        double mean_first_ring = 0.0;
        for (int j = 0; j < Nt_; ++j) mean_first_ring += at(1, j);
        mean_first_ring /= Nt_;
        for (int j = 0; j < Nt_; ++j) {
            max_change = std::max(max_change, std::abs(at(0, j) - mean_first_ring));
            at(0, j) = mean_first_ring;
        }

        // Interior update (Jacobi-like with optional SOR relaxation)
        for (int i = 1; i <= Nr_ - 2; ++i) {
            const double r   = i * dr_;
            const double rpp = r + 0.5 * dr_;
            const double rmm = r - 0.5 * dr_;
            const double cth = 1.0 / (r * r * dth_ * dth_);

            for (int j = 0; j < Nt_; ++j) {
                const double phi_im = at(i - 1, j);
                const double phi_ip = at(i + 1, j);
                const double phi_jm = at(i, jm(j));
                const double phi_jp = at(i, jp(j));

                // Discrete (1/r)∂r(r∂rφ) + (1/r^2)∂θθφ = 0
                const double A   = (rpp + rmm) / (r * dr_ * dr_) + 2.0 * cth;
                const double rhs = (rpp * phi_ip + rmm * phi_im) / (r * dr_ * dr_)
                                 + cth * (phi_jp + phi_jm);

                double phi_new = rhs / A;
                if (useSOR_) {
                    const double phi_old = at(i, j);
                    phi_new = phi_old + omega_ * (phi_new - phi_old);
                }
                max_change = std::max(max_change, std::abs(phi_new - at(i, j)));
                at(i, j) = phi_new;
            }
        }

        // Re-enforce BC
        enforceOuterBC();

        if (max_change < tol) return iter + 1;
    }
    return max_iters;
} */

int LiteLaplaceCircle::run(int max_iters, double tol) {
    // Precompute per-ring radial coefficients (independent of theta)
    struct Coeff {
        double ar, br, cth, invA;
    };
    std::vector<Coeff> cf(Nr_);
    for (int i = 1; i <= Nr_ - 2; ++i) {
        const double r   = i * dr_;
        const double rpp = r + 0.5 * dr_;
        const double rmm = r - 0.5 * dr_;

        const double ar   = rpp / (r * dr_ * dr_);
        const double br   = rmm / (r * dr_ * dr_);
        const double cth  = 1.0 / (r * r * dth_ * dth_);
        const double A    = (rpp + rmm) / (r * dr_ * dr_) + 2.0 * cth; // = ar+br+2*cth
        const double invA = 1.0 / A;

        cf[i] = { ar, br, cth, invA };
    }

    for (int iter = 0; iter < max_iters; ++iter) {
        double max_change = 0.0;

        // Enforce outer boundary (Dirichlet on electrodes; Neumann via copy) each sweep.
        enforceOuterBC();

        // Center regularity: φ(0,θ) = mean(φ(1,θ))
        double mean_first_ring = 0.0;
        for (int j = 0; j < Nt_; ++j) mean_first_ring += at(1, j);
        mean_first_ring /= Nt_;
        for (int j = 0; j < Nt_; ++j) {
            const double old0 = at(0, j);
            const double diff = std::abs(old0 - mean_first_ring);
            if (diff > max_change) max_change = diff;
            at(0, j) = mean_first_ring;
        }

        // Red–black Gauss–Seidel (in-place), optional SOR
        for (int color = 0; color < 2; ++color) {
            for (int i = 1; i <= Nr_ - 2; ++i) {
                //const auto [ar, br, cth, invA] = cf[i];
                const double ar  = cf[i].ar;
                const double br  = cf[i].br;
                const double cth = cf[i].cth;
                const double invA= cf[i].invA;

                for (int j = 0; j < Nt_; ++j) {
                    if (((i + j) & 1) != color) continue;

                    const double phi_im = at(i - 1, j);
                    const double phi_ip = at(i + 1, j);     // at(i+1, j) is valid; for i=Nr_-2 this is the boundary ring already set by enforceOuterBC()
                    const double phi_jm = at(i, jm(j));
                    const double phi_jp = at(i, jp(j));

                    const double rhs     = ar * phi_ip + br * phi_im + cth * (phi_jp + phi_jm);
                    const double phi_old = at(i, j);
                    double phi_new       = rhs * invA;

                    if (useSOR_) {
                        phi_new = phi_old + omega_ * (phi_new - phi_old);
                    }

                    const double delta = std::abs(phi_new - phi_old);
                    if (delta > max_change) max_change = delta;

                    at(i, j) = phi_new;
                }
            }

            // Keep the boundary consistent between colors
            enforceOuterBC();
        }

        // Final boundary touch this iteration
        enforceOuterBC();

        if (max_change < tol) return iter + 1;
    }
    return max_iters;
}

double LiteLaplaceCircle::getPotential(double x, double y) const {
    const double r = std::hypot(x, y);
    if (r >= R_) {
        int j = thetaToIndex(std::atan2(y, x));
        return atc(Nr_ - 1, j);
    }
    const double fr = r / dr_;
    int i0 = std::clamp(static_cast<int>(std::floor(fr)), 0, Nr_ - 2);
    const double tr = fr - i0;

    const double th = wrapAngle(std::atan2(y, x));
    const double ft = th / dth_;
    int j0 = static_cast<int>(std::floor(ft)) % Nt_; if (j0 < 0) j0 += Nt_;
    const double tt = ft - std::floor(ft);
    const int j1 = (j0 + 1) % Nt_;

    const double p00 = atc(i0,     j0);
    const double p10 = atc(i0 + 1, j0);
    const double p01 = atc(i0,     j1);
    const double p11 = atc(i0 + 1, j1);

    const double p0 = p00 * (1 - tr) + p10 * tr;
    const double p1 = p01 * (1 - tr) + p11 * tr;
    return p0 * (1 - tt) + p1 * tt;
}
