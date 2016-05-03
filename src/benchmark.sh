#!/bin/sh
FS="1e6 2e6 10e6 28e6"
N="1024 2048 4096 8192 16384 1M"
OVL="0.0 0.1 0.5"
OUTF=${1:-results.txt}
SEP="-------------------------------------------------------------------------"

cat uname -a > ${OUTF}
cat /proc/cpuinfo >> ${OUTF}
cat /proc/meminfo >> ${OUTF}
echo ${SEP} >> ${OUTF}
for fs in ${FS} ; do
	for n in ${N} ; do
		for ovl in ${OVL} ; do
			echo $fs, $n, $ovl
			echo $fs samples/sec, ${n}-point FFT, ${ovl} overlap >> ${OUTF}
			./fftw-calc $n $fs $ovl 0 >> ${OUTF}
			echo "" >> ${OUTF}
		done
	done
done
echo wrote ${OUTF}

