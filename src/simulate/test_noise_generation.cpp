//------------------------------------------------------------------------------
//
//   Standalone test for Quantum Noise Generation
//   Based on src/simulate/llg_quantum.cpp
//
//   Compile with: g++ -o test_noise test_noise_generation.cpp -lfftw3 -lm
//
//------------------------------------------------------------------------------

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <random>
#include <complex>
#include <iomanip>
#include <fstream>
#include <fftw3.h>

// Configuration parameters
struct SimulationParams {
    double A;
    double Gamma;
    double omega0;
    double S0;
    double T;
    int noise_type; // 0: Classical, 1: Quantum, 2: Semiquantum
    double dt;
    int N_steps;
};

// Global params for the functions to access
SimulationParams params;

// Function prototypes
double PSD(const double& omega, const double& T);
double estimate_cutoff_omega_cdf(double T, double target_frac);
void calculate_random_fields(int realizations, int n_fine, double dt_fine, int M, double T, int n_coarse, std::vector<double>& noise_field);

// Power Spectral Density function
double PSD(const double& omega, const double& T) {
    const double A = params.A;
    const double Gamma = params.Gamma;
    const double omega0 = params.omega0;

    double x = (T > 1e-12) ? omega / (2 * T) : omega;
    double lorentzian_denom = (omega0 * omega0 - omega * omega) * (omega0 * omega0 - omega * omega) + Gamma * Gamma * omega * omega;
    if (lorentzian_denom < 1e-12) lorentzian_denom = 1e-12;
    double lorentzian = A * Gamma * omega / lorentzian_denom;
    double coth = (x < 1e-10) ? 1.0 / x : 1.0 / tanh(x);

    switch (params.noise_type) {
        case 0: // Classical Noise
            return 2 * T * A * Gamma / lorentzian_denom;
        case 1: // Quantum Noise
            if (omega == 0) return 2 * T * A * Gamma / (omega0 * omega0 * omega0 * omega0);
            else return coth * lorentzian;
        case 2: // Semiquantum Noise
            if (omega == 0) return 2 * T * A * Gamma / (omega0 * omega0 * omega0 * omega0);
            else return (coth - 1) * lorentzian;
        default:
            return 1.0 / x * lorentzian;
    }
}

// Estimate cutoff frequency
double estimate_cutoff_omega_cdf(double T, double target_frac) {
    const double omega0 = params.omega0;
    if (omega0 <= 0) return 1.0;
    double omega_max = 10.0 * omega0;
    int steps = 50000;
    std::vector<double> psd_vals(steps + 1);
    double domega = omega_max / steps;
    for (int i = 0; i <= steps; ++i) {
        double omega = i * domega;
        psd_vals[i] = PSD(omega, T);
    }
    double total_area = 0.0;
    for (int i = 0; i < steps; ++i) {
        total_area += 0.5 * (psd_vals[i] + psd_vals[i + 1]) * domega;
    }
    if (total_area <= 1e-12) return omega_max;
    double cum_area = 0.0;
    for (int i = 0; i <= steps; ++i) {
        if (i > 0) cum_area += 0.5 * (psd_vals[i - 1] + psd_vals[i]) * domega;
        if (cum_area >= target_frac * total_area) {
            return i * domega;
        }
    }
    return omega_max;
}

// Calculate random fields using FFTW
void calculate_random_fields(int realizations, int n_fine, double dt_fine, int M, double T, int n_coarse, std::vector<double>& noise_field) {
    if (n_coarse < 0) {
        n_coarse = (n_fine > 0) ? ((n_fine - 1) / M + 1) : 0;
    }

    const double S0 = params.S0;
    const double inv_sqrt_S0 = (S0 > 0) ? 1.0 / std::sqrt(S0) : 1.0;

    const double dt_coarse = dt_fine * M;

    double* in = (double*)fftw_malloc(sizeof(double) * n_coarse);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (n_coarse / 2 + 1));
    double* result = (double*)fftw_malloc(sizeof(double) * n_coarse);

    fftw_plan forward = fftw_plan_dft_r2c_1d(n_coarse, in, out, FFTW_MEASURE);
    fftw_plan backward = fftw_plan_dft_c2r_1d(n_coarse, out, result, FFTW_MEASURE);

    const double norm_factor = 1.0 / n_coarse;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0, 1.0 / std::sqrt(dt_fine));

    noise_field.resize(realizations * n_coarse);

    std::vector<double> sqrt_PSD_coarse(n_coarse / 2 + 1);
    double df_coarse = 1.0 / (n_coarse * dt_coarse);
    for (int i = 0; i <= n_coarse / 2; ++i) {
        double omega = 2.0 * M_PI * i * df_coarse;
        sqrt_PSD_coarse[i] = std::sqrt(PSD(omega, T));
    }

    for (int r = 0; r < realizations; ++r) {
        for (int i = 0; i < n_coarse; ++i) {
            in[i] = dist(gen);
        }

        fftw_execute(forward);

        for (int i = 0; i <= n_coarse / 2; i++) {
            const double magnitude = sqrt_PSD_coarse[i];
            out[i][0] *= magnitude;
            out[i][1] *= magnitude;
        }

        fftw_execute(backward);

        const double scale = std::sqrt(dt_fine / dt_coarse);
        for (int j = 0; j < n_coarse; ++j) {
            size_t index = (size_t)j + (size_t)r * (size_t)n_coarse;
            noise_field[index] = result[j] * norm_factor * inv_sqrt_S0 * scale;
        }
    }

    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    fftw_free(in);
    fftw_free(out);
    fftw_free(result);
}

int main(int argc, char* argv[]) {
    // Default parameters
    params.A = 1.0;
    params.Gamma = 0.1;
    params.omega0 = 1.0;
    params.S0 = 1.0;
    params.T = 300.0;
    params.noise_type = 1; // 1: Quantum
    params.dt = 1e-15;
    params.N_steps = 10000;

    // Simple command line argument parsing (optional override)
    if (argc > 1) params.N_steps = std::atoi(argv[1]);
    if (argc > 2) params.dt = std::atof(argv[2]);
    if (argc > 3) params.A = std::atof(argv[3]);
    if (argc > 4) params.Gamma = std::atof(argv[4]);
    if (argc > 5) params.omega0 = std::atof(argv[5]);
    if (argc > 6) params.T = std::atof(argv[6]);
    if (argc > 7) params.noise_type = std::atoi(argv[7]);

    // Calculate derived parameters
    double omega_cutoff = estimate_cutoff_omega_cdf(params.T, 0.99999);
    int M_decimation = static_cast<int>(std::ceil((M_PI / omega_cutoff) / params.dt));
    int n_coarse = (params.N_steps > 0) ? ((params.N_steps - 1) / M_decimation + 1) : 0;
    
    std::cout << "Simulation Parameters:" << std::endl;
    std::cout << "  N_steps: " << params.N_steps << std::endl;
    std::cout << "  dt: " << params.dt << std::endl;
    std::cout << "  A: " << params.A << std::endl;
    std::cout << "  Gamma: " << params.Gamma << std::endl;
    std::cout << "  omega0: " << params.omega0 << std::endl;
    std::cout << "  T: " << params.T << std::endl;
    std::cout << "  Noise Type: " << params.noise_type << " (0:Classical, 1:Quantum, 2:Semiquantum)" << std::endl;
    std::cout << "Derived Parameters:" << std::endl;
    std::cout << "  M_decimation: " << M_decimation << std::endl;
    std::cout << "  n_coarse: " << n_coarse << std::endl;

    std::vector<double> noise_field;
    // Generate 1 realization
    calculate_random_fields(1, params.N_steps, params.dt, M_decimation, params.T, n_coarse, noise_field);

    // Output to file
    std::string filename = "noise_output.dat";
    std::ofstream outfile(filename);
    if (outfile.is_open()) {
        outfile << "# Time NoiseValue\n";
        for (size_t i = 0; i < noise_field.size(); ++i) {
            outfile << i * params.dt * M_decimation << " " << noise_field[i] << "\n";
        }
        outfile.close();
        std::cout << "Noise data written to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open output file " << filename << std::endl;
    }

    return 0;
}
