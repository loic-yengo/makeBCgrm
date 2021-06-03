makeBCgrm: 	makeBCgrm.cpp
	g++ -I"/afm01/UQ/Q3007/loic_on_imbshare/ReadPLINKBin/eigen.3.3.4/" -fopenmp -Wall -g -O3 makeBCgrm.cpp -o makeBCgrm -lm
clean:
	rm makeBCgrm

