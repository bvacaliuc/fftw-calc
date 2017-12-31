// fftw-io.c - A Tool for demonstrating the input/output loading/unloading of FFTs
//
// Acknowledgements:
//   MIT Haystack Software Developers
//     https://www.haystack.mit.edu/edu/undergrad/srt/newsrtsource_ver7.tar.gz
//     https://www.haystack.mit.edu/edu/undergrad/srt/
//   Michael Bader
//     https://www5.in.tum.de/lehre/vorlesungen/asc/ss12/rdft.pdf
//   J. Fessler
//     https://web.eecs.umich.edu/~fessler/course/451/l/pdf/c6.pdf
//   The Scientist and Engineer's Guide to Digital Signal Processing (chapter 12)
//     http://www.dspguide.com/CH12.PDF
//   Robert Matusiak, Texas Instruments, SPRA291
//     http://www.ti.com/lit/an/spra291/spra291.pdf
//   Enrico M. Staderini
//     http://www.staderini.eu/wordpress/?page_id=8
//   J. Shima
//     http://www.hyperdynelabs.com/dspdude/papers/COMPUTING%20THE%20FFT%20OF%20TWO%20REAL%20SIGNALS%20USING%20A%20SINGLE%20FFT.pdf

#include <stdio.h>
#include <stdlib.h> // for exit()
#include <ctype.h>  // for tolower()
#include <errno.h>
#include <fftw3.h>

// How many lines of context to produce around the edge conditions
#define CONTEXT 4

// from vspectra_fftw.c (101-103) loading the fftw_complex input array
// with complex samples in support of the RTL-SDR data stream.
void input_complex(fftw_complex *in, size_t N, const char *tag)
{
    size_t i, blsiz2=N, blsiz=N*2;              // NB: differs from real
    uint32_t bufferRead[blsiz];                 // NB: vspectra_fftw.c usage

    fprintf(stdout,"\n# %s\n", tag);
    fprintf(stdout,"# re{sample #},im{sample #}\n");

    for (i = 0; i < blsiz; i++)
        bufferRead[i] = (uint32_t) (i);         // NB: value is sample #
    for (i = 0; i < blsiz2; i++) {
        in[i][0] = (double) (bufferRead[2 * i]);
        in[i][1] = (double) (bufferRead[2 * i + 1]);

        if(i<CONTEXT || i>(blsiz2-CONTEXT))
            fprintf(stdout,"%.0lf,%.0lf\n", in[i][0], in[i][1]);
        else if(i==CONTEXT)
            fprintf(stdout,"...\n");
    }
}

// from vspectra_pci_fftw.c (128-133) loading the fftw_complex input array
// with real samples in support of sampled data from a PCI-2040-12 ADC card.
void input_real(fftw_complex *in, size_t N, const char *tag)
{
    size_t i, j, kk, kkk, blsiz2=N/2, blsiz=N;  // NB: differs from complex
    uint32_t buffer1[blsiz*2];                  // NB: vspectra_pci_fftw.c usage

    fprintf(stdout,"\n# %s\n", tag);
    fprintf(stdout,"# re{sample #},im{sample #}\n");

    for (i = 0; i < blsiz*2; i++)
        buffer1[i] = (uint32_t) (i);            // NB: value is sample #
    for (kk = 0; kk < 2; kk++) {
        kkk = kk * blsiz;
        for (j = 0; j < blsiz; j++) {
            // NB: fft is only done when (kk%2)==1
            if (kk % 2 == 0)
                in[j][0] = (double) (buffer1[j + kkk]);
            else
                in[j][1] = (double) (buffer1[j + kkk]);
        }
    }
    for (i = 0; i < blsiz; i++) {
        if(i<CONTEXT || i>(blsiz-CONTEXT))
            fprintf(stdout,"%.0lf,%.0lf\n", in[i][0], in[i][1]);
        else if(i==CONTEXT)
            fprintf(stdout,"...\n");
    }
}

void dump(fftw_complex *x, size_t N, const char *tag)
{
    size_t ii;

    fprintf(stderr,"# %s (complete list of FFT re/im data)\n", tag);
    fprintf(stderr,"# re{n},im{n}\n");

    for(ii=0;ii<N;ii++) fprintf(stderr,"%.6lf,%.6lf\n", x[ii][0], x[ii][1]);
}

// from vspectra_fftw.c (112-124) unloading the fftw_complex output array
// with power spectrum estimate in support of the RTL-SDR data stream.
void output_complex(fftw_complex *out, size_t N, const char *tag, int prec)
{
    size_t i, j, blsiz2=N, blsiz=N*2;           // NB: differs from real
    double vspec[blsiz2], scale;
    int numm = 0;

    fprintf(stdout,"# %s\n", tag);
    if(prec<0)
        fprintf(stdout,"# formula to produce vspec{n}\n");
    else
        fprintf(stdout,"# re{n},im{n},vspec{n}\n");

    for (i = 0; i < blsiz2; i++) {
        vspec[i] = 0.0;
    }
    for (i = 0; i < blsiz2; i++) {
        if (i < blsiz2 / 2)
            j = i + blsiz2 / 2;
        else
            j = i - blsiz2 / 2;
        vspec[j] += out[i][0] * out[i][0] + out[i][1] * out[i][1];
        if(prec<0 && (i<CONTEXT || i>(N-CONTEXT) || j<CONTEXT || j>(N-CONTEXT)))
            fprintf(stdout,"vspec[%04d]=sqr(re[%d])+sqr(im[%d])\n", j, i, i);
        else if(prec<0 && (i==CONTEXT || j==CONTEXT))
            fprintf(stdout,"...\n");
        numm++;
    }
    if(prec<0) return;
    scale = 1.0/(double)numm;
    for (i = 0; i < blsiz2; i++) {
        fprintf(stdout,"%.*lf,%.*lf,%.6lf\n", prec, out[i][0], prec, out[i][1], vspec[i]*scale);
    }
}

// from vspectra_pci_fftw.c (157-159) unloading the fftw_complex output array
// with power spectrum estimate in support of the PCI-2040-12 ADC data stream.
void output_real(fftw_complex *out, size_t N, const char *tag, int prec)
{
    size_t i, kk, blsiz2=N/2, blsiz=N;          // NB: differs from complex
    double rre, aam, rre2, aam2;
    double vspec[blsiz2], scale;
    int numm = 0;

    fprintf(stdout,"# %s\n", tag);
    if(prec<0)
        fprintf(stdout,"# formula to produce vspec{n}\n");
    else
        fprintf(stdout,"# vspec{n}\n");

    for (i = 0; i < blsiz2; i++) {
        vspec[i] = 0.0;
    }
    for (kk = 0; kk < 2; kk++) {
        // NB: fft is only done when (kk%2)==1
        if (kk % 2 == 1) {
            for (i = 0; i < blsiz2; i++) {
                if (i >= 1) {
                    rre = out[i][0] + out[(blsiz - i)][0];
                    aam = out[i][1] - out[(blsiz - i)][1];
                    aam2 = -out[i][0] + out[(blsiz - i)][0];
                    rre2 = out[i][1] + out[(blsiz - i)][1];
                    if(prec<0 && (i<CONTEXT || i>(blsiz2-CONTEXT)))
                        fprintf(stdout,"vspec[%04d]=sqr(%s[%d]%s[%d])+sqr(%s[%d]%s[%d])+sqr(%s[%d]%s[%d])+sqr(%s[%d]%s[%d])\n",
                            i,
                            "re", i, "+re", (blsiz - i),
                            "im", i, "-im", (blsiz - i),
                            "-re", i,"+re", (blsiz - i),
                            "im", i, "+im", (blsiz - i));
                    else if(prec<0 && i==CONTEXT)
                        fprintf(stdout,"...\n");
                } else {
                    rre = out[i][0] + out[0][0];
                    aam = out[i][1] - out[0][1];
                    aam2 = -out[i][0] + out[0][0];
                    rre2 = out[i][1] + out[0][1];
                    if(prec<0)
                        fprintf(stdout,"vspec[%04d]=sqr(%s[%d]%s[%d])+sqr(%s[%d]%s[%d])+sqr(%s[%d]%s[%d])+sqr(%s[%d]%s[%d])\n",
                            i,
                            "re", i, "+re", 0,
                            "im", i, "-im", 0,
                            "-re", i,"+re", 0,
                            "im", i, "+im", 0);
                }
                vspec[i] += rre * rre + aam * aam + rre2 * rre2 + aam2 * aam2;
            }
        }
    }
    if(prec<0) return;
    // TODO: strange, the vspectra_pci_fftw.c code only increments the scaling factor
    // once outside of the kk loop (i.e. it divides by two every 1M data points...)
    numm++;
    scale = 1.0/(double)numm;
    for (i = 0; i < blsiz2; i++) {
        fprintf(stdout,"%.*lf,%.*lf,%.6lf\n", prec, out[i][0], prec, out[i][1], vspec[i]*scale);
    }
}

void usage( const char *argv0 )
{
    printf("Usage:"
    "\n%s [N]"
    "\nWhere:"
    "\n\tN:\tFFT length in points (default 1024)"
    "\nDemonstrates the input loading and output interpretation of complex"
    "\nand real value streams as inputs to a complex FFT.  The program dumps"
    "\nthe array indices used depending on complex or real-valued input and"
    "\nthe formula needed to compute the power spectral density estimate."
    "\n\n", argv0);
    exit(1);
}

void main(int argc, char *argv[], char *envp[])
{
    fftw_complex *in, *out;
    fftw_plan p;
    size_t N = 1024;
    char tag[80];
    char *endp = NULL;

    if(argc>1 && ( *argv[1]=='?' || *argv[1]=='-' ) ) usage(argv[0]);

    errno = 0;
    if(argc>1) N = strtoul(argv[1],&endp,0); if(errno) perror("N"); errno = 0;
         if(endp && tolower(*endp)=='k') N *= 1024;
    else if(endp && tolower(*endp)=='m') N *= 1024*1024;

    // DEMONSTRATE taking the N-pt FFT of complex value input
    // and extracting an N-pt power spectrum estimate from the result

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // vspectra_fftw.c (taking input from COMPLEX samples from RTL-SDR dongle)
    snprintf(tag,sizeof(tag),"FFT INPUT (%d I+Q values, %d pt FFT)", N, N);
    input_complex(in, N, tag);
    //dump(in,N,"INPUT-COMPLEX");
    fftw_execute(p);
    //dump(out,N,"OUTPUT-COMPLEX");
    snprintf(tag,sizeof(tag),"FFT OUTPUT (%d spectral bins, %d pt FFT)", N, N);
    output_complex(out, N, tag, -1);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    // DEMONSTRATE taking the 2*N-pt FFT of real value input
    // and extracting an N-pt power spectrum estimate from the result

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * 2);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * 2);
    p = fftw_plan_dft_1d(N * 2, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // vspectra_pci_fftw.c (taking input from REAL samples from a single ADC)
    snprintf(tag,sizeof(tag),"FFT INPUT (%d real values, %d pt FFT)", 4*N, 2*N);
    input_real(in, N * 2, tag );
    //dump(in,N,"INPUT-REAL");
    fftw_execute(p);
    //dump(out,N,"OUTPUT-REAL");
    snprintf(tag,sizeof(tag),"FFT OUTPUT (%d spectral bins, %d pt FFT)", N, 2*N);
    output_real(out, N * 2, tag, -1);

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
}
