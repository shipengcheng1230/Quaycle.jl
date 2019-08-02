## Global Performace Settings

const parameters = Dict(
    "FFT" => Dict(
        # FFTW threads
        "NUM_THREADS" => Threads.nthreads(),
        # FFTW measurement flag
        "FLAG" => FFTW.MEASURE,
    ),
    "BLAS" => Dict(
        # BLAS threads
        "NUM_THREADS" => Threads.nthreads(),
    )
)

BLAS.set_num_threads(parameters["BLAS"]["NUM_THREADS"])
FFTW.set_num_threads(parameters["FFT"]["NUM_THREADS"])
