## Global Performace Settings

fft_configs = Dict(
    # threads for fft
    "FFT_NUM_THREADS" => nthreads(),

    # fft measurement flag
    "FFT_FLAG" => FFTW.MEASURE,
)
