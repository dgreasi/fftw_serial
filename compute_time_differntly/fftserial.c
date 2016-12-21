#include "libraries.h"

typedef struct comp_comp {
	float Re;
	float Im;
} complex;

int bitrev(int inp, int numbits)
{
  int i, rev=0;
  //printf("Original int %d\n",inp);
  for (i=0; i < numbits; i++)
  {
    rev = (rev << 1) | (inp & 1);
    inp >>= 1;
  }
 // printf("Reversed int %d\n",rev);
  return rev;
}

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

double serial_FFT(complex *X, complex *Y, long n) {
	long ii;
	long r, temp, i, j, k, m, w;
	complex omeg;
	complex  *R,*S;
	double ext;
	
	R = (complex *) malloc(n*sizeof(complex));
	S = (complex *) malloc(n*sizeof(complex));
	
	clock_t begin=clock();	
	/* Calculate r=logn with n=2^r */	
	r=0;
	temp=n;
	while ( (n /= 2 ) != 0 ){
		r++;}
	n=temp;
	
	for (i=0; i<n; i++){
		R[i].Re = X[i].Re;	
		R[i].Im = X[i].Im;
	}
	
	/* Compute the FFT */
	
	for (m=0; m<r; m++) {
		
		for (i=0; i<n; i++){
			S[i].Re = R[i].Re;	
			S[i].Im = R[i].Im;
			//printf("step %ld : S has %fl \n",m,S[i].Re);
		}
		
		for (ii=0;ii<n;ii++) {
	/* Compute the appropriate power of omega based on bitwise arithmetic of indices */
		j=(ii & (~(1 << (r-m-1)))) | (0 << (r-m-1));
		k=(ii & (~(1 << (r-m-1)))) | (1 << (r-m-1));
		w=bitrev(ii,r);
		w =w << (r-1-m);
		omeg=omega(n,-1,w);
		
		R[ii]=csum(S[j],cmul(omeg,S[k]));
	}
}
	
	/* Scale and reverse indices to obtain the correct transformed values */
	
	
	for (i=0;i<n-1;i++) {
		Y[i]=R[bitrev(i,r)];
	}
	clock_t end=clock();
	ext=(double)(end-begin)/CLOCKS_PER_SEC;
	free(R);
	free(S);
	return ext;
}

int main(int argc, char** argv) {	
	int size, i;
	complex *X, *Y;
	double sum,mean,ext;
	//Provide power of size at input
    	size = pow(2,atoi(argv[1]));
	
	Y = (complex *) malloc(size*sizeof(complex));
	X = (complex *) malloc(size*sizeof(complex));
	for(i=0;i<size;i++){
	
		X[i].Re=(float)rand();
		X[i].Im=(float)rand();
	}
	sum=0;
	ext=0;
	for(i=0; i<1000; i++){
	//clock_t begin=clock();
	ext=serial_FFT(X,Y,size);
	//clock_t end=clock();
	sum=sum+ext;
	//double time_spent=(double)(end-begin) / CLOCKS_PER_SEC;
	//sum=sum+time_spent;
	}
	mean=(double) sum/ (double) 1000;
	//double time_spent=(double)(end-begin) / CLOCKS_PER_SEC;
	printf("Time spent in serial_FFT for %d elements: %lf\n",size,mean);
	
	//Printing the values of Y for debugging
	//for(i=0;i<size;i++) {	
	//	printf("%fl %fl",Y[i].Re,Y[i].Im); }		

	free(X);
	free(Y);
	return 0;
}