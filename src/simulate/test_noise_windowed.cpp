//------------------------------------------------------------------------------
//
//   Windowed Noise Generation Test
//   Based on src/simulate/test_noise_generation.cpp
//
//   Compile with: g++ -o test_noise_windowed test_noise_windowed.cpp -lfftw3 -lm
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
#include <cstring>

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
    int window_size; // Size of the small window (must be divisible by 6)
};

SimulationParams params;

// Function prototypes
double PSD(const double& omega, const double& T);
double estimate_cutoff_omega_cdf(double T, double target_frac);
void calculate_random_fields_windowed(int realizations, int n_fine, double dt_fine, int M, double T, int n_coarse_total, std::vector<double>& noise_field);

// Power Spectral Density function (Same as before)
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

// Estimate cutoff frequency (Same as before)
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

// Windowed Random Field Calculation
void calculate_random_fields_windowed(int realizations, int n_fine, double dt_fine, int M, double T, int n_coarse_total, std::vector<double>& noise_field) {
    
    int window_size = params.window_size;
    if (window_size % 6 != 0) {
        std::cerr << "Error: Window size must be divisible by 6." << std::endl;
        return;
    }
    int segment_size = window_size / 6;

    const double S0 = params.S0;
    const double inv_sqrt_S0 = (S0 > 0) ? 1.0 / std::sqrt(S0) : 1.0;
    const double dt_coarse = dt_fine * M;
    const double norm_factor = 1.0 / window_size;
    const double scale = std::sqrt(dt_fine / dt_coarse);

    // FFTW Setup for the window size
    double* in = (double*)fftw_malloc(sizeof(double) * window_size);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (window_size / 2 + 1));
    double* result = (double*)fftw_malloc(sizeof(double) * window_size);

    fftw_plan forward = fftw_plan_dft_r2c_1d(window_size, in, out, FFTW_MEASURE);
    fftw_plan backward = fftw_plan_dft_c2r_1d(window_size, out, result, FFTW_MEASURE);

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0, 1.0 / std::sqrt(dt_fine));

    // Precompute PSD for the window
    std::vector<double> sqrt_PSD_window(window_size / 2 + 1);
    double df_window = 1.0 / (window_size * dt_coarse);
    for (int i = 0; i <= window_size / 2; ++i) {
        double omega = 2.0 * M_PI * i * df_window;
        sqrt_PSD_window[i] = std::sqrt(PSD(omega, T));
    }

    // Buffer to hold white noise
    std::vector<double> white_noise_buffer(window_size);

    // Resize output vector
    noise_field.clear();
    noise_field.reserve(n_coarse_total * realizations);

    for (int r = 0; r < realizations; ++r) {
        int generated_samples = 0;
        bool first_run = true;

        while (generated_samples < n_coarse_total) {
            
            if (first_run) {
                // Fill entire buffer with new random numbers
                for (int i = 0; i < window_size; ++i) {
                    white_noise_buffer[i] = dist(gen);
                }
            } else {
                // Shift buffer: Move last 2 segments (5 and 6) to front (1 and 2)
                // Indices: 4*segment_size to 6*segment_size-1  --> 0 to 2*segment_size-1
                std::memmove(white_noise_buffer.data(), 
                             white_noise_buffer.data() + 4 * segment_size, 
                             2 * segment_size * sizeof(double));
                
                // Fill the rest (segments 3, 4, 5, 6) with new random numbers
                for (int i = 2 * segment_size; i < window_size; ++i) {
                    white_noise_buffer[i] = dist(gen);
                }
            }

            // Copy to FFT input
            for(int i=0; i<window_size; ++i) in[i] = white_noise_buffer[i];

            // FFT
            fftw_execute(forward);

            // Apply PSD
            for (int i = 0; i <= window_size / 2; i++) {
                const double magnitude = sqrt_PSD_window[i];
                out[i][0] *= magnitude;
                out[i][1] *= magnitude;
            }

            // Inverse FFT
            fftw_execute(backward);

            // Process Output
            if (first_run) {
                // Take segments 1-5
                int count = 5 * segment_size;
                for (int j = 0; j < count; ++j) {
                    if (noise_field.size() < (size_t)(n_coarse_total * (r + 1))) {
                         noise_field.push_back(result[j] * norm_factor * inv_sqrt_S0 * scale);
                    }
                }
                generated_samples += count;
                first_run = false;
            } else {
                // Take segments 2-5 (indices 1*segment_size to 5*segment_size - 1)
                // This corresponds to the "2-4th" request interpreted as 2,3,4,5 to match the shift of 4 segments.
                // Segment 1 (index 0..S-1) is discarded (overlap with previous Segment 5).
                // Segment 2 (index S..2S-1) is the recalculated previous Segment 6.
                // Segments 3,4,5,6 are new. We output 2,3,4,5. Segment 6 is kept for next overlap.
                
                int start_idx = 1 * segment_size; // Start of segment 2
                int end_idx = 5 * segment_size;   // End of segment 5 (exclusive)
                
                for (int j = start_idx; j < end_idx; ++j) {
                    if (noise_field.size() < (size_t)(n_coarse_total * (r + 1))) {
                        noise_field.push_back(result[j] * norm_factor * inv_sqrt_S0 * scale);
                    }
                }
                generated_samples += (end_idx - start_idx);
            }
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
    params.window_size = 6000; // Default window size

    // Simple command line argument parsing
    if (argc > 1) params.N_steps = std::atoi(argv[1]);
    if (argc > 2) params.dt = std::atof(argv[2]);
    if (argc > 3) params.A = std::atof(argv[3]);
    if (argc > 4) params.Gamma = std::atof(argv[4]);
    if (argc > 5) params.omega0 = std::atof(argv[5]);
    if (argc > 6) params.T = std::atof(argv[6]);
    if (argc > 7) params.noise_type = std::atoi(argv[7]);
    if (argc > 8) params.window_size = std::atoi(argv[8]);

    // Calculate derived parameters
    double omega_cutoff = estimate_cutoff_omega_cdf(params.T, 0.99999);
    int M_decimation = static_cast<int>(std::ceil((M_PI / omega_cutoff) / params.dt));
    int n_coarse = (params.N_steps > 0) ? ((params.N_steps - 1) / M_decimation + 1) : 0;
    
    std::cout << "Windowed Simulation Parameters:" << std::endl;
    std::cout << "  N_steps: " << params.N_steps << std::endl;
    std::cout << "  dt: " << params.dt << std::endl;
    std::cout << "  Window Size: " << params.window_size << std::endl;
    std::cout << "  M_decimation: " << M_decimation << std::endl;
    std::cout << "  n_coarse: " << n_coarse << std::endl;

    std::vector<double> noise_field;
    calculate_random_fields_windowed(1, params.N_steps, params.dt, M_decimation, params.T, n_coarse, noise_field);

    // Output to file
    std::string filename = "noise_windowed.dat";
    std::ofstream outfile(filename);
    if (outfile.is_open()) {
        outfile << "# Time NoiseValue\n";
        for (size_t i = 0; i < noise_field.size(); ++i) {
            outfile << i * params.dt * M_decimation << " " << noise_field[i] << "\n";
        }
        outfile.close();
        std::cout << "Windowed noise data written to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open output file " << filename << std::endl;
    }

    return 0;
}
