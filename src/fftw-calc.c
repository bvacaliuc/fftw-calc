// fftw-calc.c - A simple Tool for estimating computational complexity of streaming FFTs
//
// Acknowledgements:
//   Matteo Frigo and Steven G. Johnson.
//     http://www.fftw.org/fftw3_doc/Using-Plans.html#Using-Plans
//   Pioneer Hill Software, LLC
//     https://www.spectraplus.com/DT_help/overlap_percentage.htm
//   jonolafur (http://sourceforge.net/projects/cli-fft/)
//   Pau Garcia i Quiles (http://www.elpauer.org/stuff/learning_cmake.pdf)
//   Alexander J. Yee
//     https://github.com/Mysticial/Flops
//     http://stackoverflow.com/questions/8389648/how-to-achieve-4-flops-per-cycle

#include <stdio.h>
#include <stdlib.h>
//#include <stdbool.h>  // NB: dealt with in cpuid.h, may need to refactor/rethink 'bool'
#include <ctype.h>  // for tolower()
#include <errno.h>
#include <time.h>   // for clock()
#include <fftw3.h>
#include "cpuid.h"

// http://stackoverflow.com/questions/1546789/clean-code-to-printf-size-t-in-c-or-nearest-equivalent-of-c99s-z-in-c
#if defined(_MSC_VER) || defined(__MINGW32__) //__MINGW32__ should go before __GNUC__
#define __SIZE_T_SPECIFIER    "%Iu"
#define __SSIZE_T_SPECIFIER   "%Id"
#define __PTRDIFF_T_SPECIFIER "%Id"
#else //elif defined(__GNUC__)
#define __SIZE_T_SPECIFIER    "%zu"
#define __SSIZE_T_SPECIFIER   "%zd"
#define __PTRDIFF_T_SPECIFIER "%zd"
#endif

void usage( const char *argv0 )
{
	printf("Usage:"
	"\n%s [N [FS [WIN [EST]]]]"
	"\nWhere:"
	"\n\tN:\tFFT length in points (default 1024)"
	"\n\tFS:\tSampling frequency in Hz (default 1e6/1MHz)"
	"\n\tWIN:\tOverlap Window in %% (default 0%%/no overlap)"
	"\n\tEST:\tEstimate(1) or Measure(0) (default 0/measure)"
	"\n\nExamples:"
	"\n\n%s 32768 28e6 0.1 1"
	"\nCalculates the FFT plan and estimates operations for a spectrum analyzer"
	"\noperating at 28MHz (complex samples) with an FFT size of 32K and 10%% overlap."
	"\n\n%s 32K 28MHz 0.1 0"
	"\nSimilar to above, but illustrates the use of the suffixes, and forces measurement"
	"\n\n", argv0, argv0, argv0);
	exit(1);
}

void main(int argc, char *argv[], char *envp[])
{
	fftw_complex *in, *out;
	fftw_plan p;
	size_t N = 1024;
	double _FS, FS = 1e6;	// 1MHz sample rate
	double WIN = 0.0;	// 0% window overlap
	double nadd,nmul,nfma,nflops,ntotal,factor;
	char *units, *fs_units;
	bool EST = false;
	char *endp = NULL;
    cpuid_info_t cpu;
    const char *_mp = "";

	if(argc>1 && ( *argv[1]=='?' || *argv[1]=='-' ) ) usage(argv[0]);

	errno = 0;
	if(argc>1) N   = strtoul(argv[1],&endp,0); if(errno) perror("N"); errno = 0;
	     if(endp && tolower(*endp)=='k') N *= 1024;
	else if(endp && tolower(*endp)=='m') N *= 1024*1024;
	if(argc>2) FS  = strtod(argv[2],&endp); if(errno) perror("FS(Hz)"); errno = 0;
	     if(endp && *endp=='k') FS *= 1000.0;	// common use would be to qualify with a 2nd character
	else if(endp && *endp=='M') FS *= 1000.0*1000.0;
	else if(endp && *endp=='G') FS *= 1000.0*1000.0*1000.0;
	else if(endp && *endp=='T') FS *= 1000.0*1000.0*1000.0*1000.0;  // I'm dreaming of the day...
	else if(endp && *endp=='m') FS /= 1000.0;	// ok, kind of silly
	if(argc>3) WIN = strtod(argv[3],NULL); if(errno) perror("WINDOW(%)"); errno = 0;
	if(argc>4) EST = atoi(argv[4])?true:false;

    cpuid_get_info( &cpu );
    _FS = FS;

    // http://www.fftw.org/fftw3_doc/Usage-of-Multi_002dthreaded-FFTW.html#Usage-of-Multi_002dthreaded-FFTW
#if defined(_OPENMP) // || defined(_POSIX_THREADS)
    if( cpu.threads > 1 )
    {
        fftw_init_threads();
        fftw_plan_with_nthreads(cpu.threads);
        _mp = "-omp";
    }
#endif

	printf("FFT("__SIZE_T_SPECIFIER", %s) :\n",N,EST?"estimated":"measured");
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, EST?FFTW_ESTIMATE:FFTW_MEASURE);
	fftw_flops(p,&nadd,&nmul,&nfma);
	fftw_print_plan(p);
	nflops = ntotal = nadd+nmul+((cpu.exts.XOP||cpu.exts.FMA3||cpu.exts.FMA4)?nfma:2*nfma);
	
	printf("\nFLOPS: add=%.0f mul=%.0f fma=%.0f total=%.0f Flops/frame\n",
		nadd, nmul, nfma, ntotal);
	factor  = (1.0 + ( 1.0 / ( 1.0 - WIN ))) / 2.0 ;	// additional FFTs required due to overlap
	ntotal *= (FS * factor) / N;						// FFTs = FS/N times the factor due to overlap
	     if( ntotal > 5e17) { ntotal /= 1e18;units = "ExaFlops"; }
	else if( ntotal > 5e14) { ntotal /= 1e15;units = "PFlops"; }
	else if( ntotal > 5e11) { ntotal /= 1e12;units = "TFlops"; }
	else if( ntotal > 5e8 ) { ntotal /= 1e9; units = "GFlops"; }
	else if( ntotal > 5e5 ) { ntotal /= 1e6; units = "MFlops"; }
	else if( ntotal > 5e2 ) { ntotal /= 1e3; units = "KFlops"; }
    else                    {                units = "Flops";  }
	     if( FS > 5e11) { FS /= 1e12;fs_units = "THz"; }
	else if( FS > 5e8 ) { FS /= 1e9; fs_units = "GHz"; }
	else if( FS > 5e5 ) { FS /= 1e6; fs_units = "MHz"; }
	else if( FS > 5e2 ) { FS /= 1e3; fs_units = "KHz"; }
	else                {            fs_units = "Hz";  }
	printf("FS=%.2f%s, %.2f%% overlap, %.2f %s (%s method)\n",
		FS, fs_units, WIN*100.0, ntotal, units, EST?"by estimate":"by measure");

#ifndef FFTW_DLL
    // TODO: doesn't work for MSVC build
    printf("FFTw Version = %s%s\n", fftw_version, _mp);
#else
    printf("FFTw Version = TBD%s\n", _mp);
#endif
    printf("Current CPU = %s\n", cpu.name.str);
    printf("CPU Threads = %d\n", cpu.threads);

	// TODO: actually compute some representative FFTs, timing them and
	//       extrapolate the performance on *THIS* machine as configured
	if( !EST )
	{
		clock_t start, stop;
		double  elapsed, fps;
		int ii, ffts = 1 * (int)((_FS * factor) / N);

		//printf("_FS=%g, factor=%g, N=%zu, ffts=%d\n", _FS, factor, N, ffts);

		// the total amount of work necessary to go through 1sec of input data...
		start = clock();
		for(ii=0;ii<ffts;ii++) fftw_execute(p); /* repeat as needed */
		stop = clock();

		elapsed = ((double)stop - (double)start) / (double)CLOCKS_PER_SEC;

		//printf("stop=%g, start=%g, elapsed=%.3f sec\n", (double)stop, (double)start, elapsed);

		fps = nflops * (double)ffts/elapsed;
		     if( fps > 5e17) { fps /= 1e18;units = "ExaFlops"; }
		else if( fps > 5e14) { fps /= 1e15;units = "PFlops"; }
		else if( fps > 5e11) { fps /= 1e12;units = "TFlops"; }
		else if( fps > 5e8 ) { fps /= 1e9; units = "GFlops"; }
		else if( fps > 5e5 ) { fps /= 1e6; units = "MFlops"; }
		else if( fps > 5e2 ) { fps /= 1e3; units = "KFlops"; }
		else                 {             units = "Flops";  }

		printf("%d FFTs in %.3f sec (%.2f %s)\n", ffts, elapsed, fps, units);
		if( elapsed > 1.1 ) printf("*** this CPU/configuration will not meet your specification ***\n");
		else if( elapsed > 0.9 ) printf("*** this configuration is close to full utilization on this CPU ***\n");
	}

	fftw_destroy_plan(p);
	fftw_free(in); fftw_free(out);
#if defined(_OPENMP) // || defined(_POSIX_THREADS)
	if( cpu.threads > 1 ) fftw_cleanup_threads();
#endif
}

