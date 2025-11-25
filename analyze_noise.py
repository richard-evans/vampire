import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def load_data(filename):
    try:
        # Assuming space-separated values with comments starting with #
        data = pd.read_csv(filename, sep='\s+', comment='#', header=None, names=['Time', 'Noise'])
        return data
    except FileNotFoundError:
        print(f"Error: File {filename} not found.")
        return None

def calculate_autocorrelation(x):
    """
    Compute the autocorrelation of the signal, based on the convolution theorem.
    """
    xp = x - np.mean(x)
    f = np.fft.fft(xp)
    p = np.absolute(f)**2
    pi = np.fft.ifft(p)
    return np.real(pi)[:len(x)//2] / np.sum(xp**2)

def main():
    file_monolithic = 'noise_output.dat'
    file_windowed = 'noise_windowed.dat'

    print(f"Loading {file_monolithic}...")
    data_mono = load_data(file_monolithic)
    
    print(f"Loading {file_windowed}...")
    data_window = load_data(file_windowed)

    if data_mono is None or data_window is None:
        return

    # --- Histograms ---
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.hist(data_mono['Noise'], bins=100, density=True, alpha=0.6, label='Monolithic', color='blue')
    plt.hist(data_window['Noise'], bins=100, density=True, alpha=0.6, label='Windowed', color='orange')
    plt.title('Noise Histogram')
    plt.xlabel('Noise Value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, alpha=0.3)

    # --- Autocorrelation ---
    print("Calculating autocorrelation for Monolithic data...")
    acf_mono = calculate_autocorrelation(data_mono['Noise'].values)
    
    print("Calculating autocorrelation for Windowed data...")
    acf_window = calculate_autocorrelation(data_window['Noise'].values)

    # Time axis for ACF (assuming constant dt)
    dt = data_mono['Time'].iloc[1] - data_mono['Time'].iloc[0]
    lags = np.arange(len(acf_mono)) * dt

    plt.subplot(1, 2, 2)
    plt.plot(lags, acf_mono, label='Monolithic', color='blue', alpha=0.8)
    plt.plot(lags, acf_window, label='Windowed', color='orange', alpha=0.8, linestyle='--')
    plt.title('Autocorrelation Function')
    plt.xlabel('Time Lag')
    plt.ylabel('Autocorrelation')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 0.5) # Zoom in on the beginning
    


    plt.tight_layout()
    plt.savefig('noise_analysis.png')
    print("Plot saved to noise_analysis.png")
    plt.show()

if __name__ == "__main__":
    main()
