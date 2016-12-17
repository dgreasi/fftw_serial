#include "libraries.h"

typedef struct comp_comp {
	float Re;
	float Im;
} complex;

//Appropriate exponential power of omega calculation: w^ki = e^(2*PI*k*i/n)
complex omega(int n, int i, int k) {
	complex omeg;

	omeg.Re = cos(k*i*2*PI/n);
	omeg.Im = sin(k*i*2*PI/n);

	return omeg;
}

// addition of 2 complex numbers c1,c2
complex csum(complex c1, complex c2){
	complex sum;

	sum.Re = c1.Re + c2.Re;
	sum.Im = c1.Im + c2.Im;
	return sum;
}

// multiplication of 2 complex numbers c1,c2
complex cmul(complex c1, complex c2){
	complex mul;

	mul.Re = c1.Re * c2.Re - c1.Im * c2.Im;
	mul.Im = c1.Re * c2.Im + c1.Im * c2.Re;
	return mul;
}

// subtraction of 2 complex numbers c1,c2
complex csub(complex c1, complex c2){
	complex sub;

	sub.Re = c1.Re - c2.Re;
	sub.Im = c1.Im - c2.Im;
	return sub;
}

void serial_FFT(complex *X, complex *Y, long n) {
	long start, end, ii, shift,prev_shift, l;
	long r, temp, i, j, k, m, w, y, t, i2;
	complex temp2, temp3, omeg;
	complex t1;
		
	/* Calculate r=logn with n=2^r */	
	r=0;
	temp=n;
	while ( (n /= 2 ) != 0 ){
		r++;}
	n=temp;
	
	for (i=0; i<n; i++){
		Y[i] = X[i];	}	
	
	/* Compute the FFT */
	l=1;
	shift = n;
	prev_shift = shift;
	for (m=r-1; m>=0; m--) {
		shift /= 2;
		start = 0;
		end = shift;
		
		for (ii=0;ii<l;ii++) {
	/* Compute the appropriate power of omega based on bitwise arithmetic of indices */

			w=0;
			for (y=0; y<r; y++) {
				w= w << 1;
				t = (start >> y) & 1;
				w = (w+t);				
			}
			w = (w << m) & (n-1);
			omeg = omega(n,-1,w);
		
			for (i=start; i<end; i++) {			
				t1 = Y[i];
				temp2 = cmul(Y[i+shift], omeg);
				Y[i + shift] = csub(t1, temp2);	
				Y[i] = csum(t1, temp2);	
			}
			start += prev_shift;
			end += prev_shift;
		}
		prev_shift=shift;
		l *= 2;
	}
	
	/* Scale and reverse indices to obtain the correct transformed values */
	for (i=0; i<n; i++) {
		Y[i].Re /= n;
		Y[i].Im /= n;
	}
	
	i2 = n >> 1;
	j = 0;
	for (i=0;i<n-1;i++) {
		if (i < j) {	
			temp3 = Y[i];
			Y[i] = Y[j];
			Y[j] = temp3;
		}
		k = i2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	
}

int main (int argc, char *argv[]) {
	int size, i;
	complex *X, *Y;
	
	//Provide power of size at input
    	size = pow(2,atoi(argv[1]));
	
	Y = (complex *) malloc(size*sizeof(complex));
	X = (complex *) malloc(size*sizeof(complex));
	for(i=0;i<size;i++){
	
		X[i].Re=(float)rand();
		X[i].Im=(float)rand();
	}
	
	clock_t begin=clock();
	serial_FFT(X,Y,size);
	clock_t end=clock();
	
	double time_spent=(double)(end-begin) / CLOCKS_PER_SEC;
	printf("Time spent in serial_FFT for %d elements: %lf\n",size,time_spent);
	
	//Printing the values of Y for debugging
	//for(i=0;i<size;i++) {	
	//	printf("%fl %fl",Y[i].Re,Y[i].Im); }		

	free(X);
	free(Y);
	return 0;

}
