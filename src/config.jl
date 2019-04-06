## Global Performace Settings

FFT_CONFIGS = Dict(
    # threads for fft
    "FFT_NUM_THREADS" => Threads.nthreads(),

    # fft measurement flag
    "FFT_FLAG" => FFTW.MEASURE,
)
