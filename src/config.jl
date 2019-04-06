## Global Performace Settings

parameters = Dict(
    "FFT" => Dict(
        # threads for fft
        "NUM_THREADS" => Threads.nthreads(),
        # fft measurement flag
        "FLAG" => FFTW.MEASURE,
    ),
)
