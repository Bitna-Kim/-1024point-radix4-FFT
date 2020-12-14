
int16_t* FFT4(int16_t *x, int16_t *y) 
{ 
	 int32_t* x32, *pcumulre, *pcumulim; 
	 __m128i x128, mmcumulre1, mmcumulim1, mmcumulre2, mmcumulim2, mmtmpre, mmtmpim, *result_FFT4;  
	 __m128i *Wptr; 


	 x32 = (int32_t*)x;
	 result_FFT4 = (__m128i*)y;
	 
	 mmcumulre1= _mm_setzero_si128();
	 mmcumulim1= _mm_setzero_si128();
	 
	 for(int i=0;i<4;i++)
	 { 
		  x128 = _mm_set_epi32(x32[i], x32[i], x32[i], x32[i]);
		  Wptr = (__m128i*)minusW4[i];
		  mmtmpre=_mm_madd_epi16(x128,Wptr[0]); 


		  Wptr = (__m128i*)shuffleW4[i];
		  mmtmpim=_mm_madd_epi16(x128,Wptr[0]); 

		  mmtmpre=_mm_srai_epi32(mmtmpre, 15);
		  mmtmpim=_mm_srai_epi32(mmtmpim, 15);

		  /*
		  printf("madd result  real (%i) : ", i);
		  dumpSimdData(&mmtmpre, 1);
		  printf("madd result  imag (%i) : ", i);
		  dumpSimdData(&mmtmpim, 1);
		  printf("\n");
		  */

		   
		  mmcumulre1 = _mm_add_epi32(mmcumulre1, mmtmpre);
		  mmcumulim1 = _mm_add_epi32(mmcumulim1, mmtmpim);

		  /*
		  printf("add result  real (%i) : ", i);
		  dumpSimdData(&mmcumulre1, 1);
		  printf("add result  imag (%i) : ", i);
		  dumpSimdData(&mmcumulim1, 1);
		  printf("\n");
		  */
  
 	 }

	 mmcumulre2=_mm_packs_epi32(mmcumulre1, mmcumulre1); 
	 mmcumulim2=_mm_packs_epi32(mmcumulim1, mmcumulim1);


	 /*
	 printf("pack result  real : ");
	 dumpSimdData(&mmcumulre2, 0);
	 printf("pack result  imag : ");
	 dumpSimdData(&mmcumulim2, 0);
	 printf("\n");
	 */

	 result_FFT4[0] =_mm_unpacklo_epi16(mmcumulre2, mmcumulim2);

	 /*
	 printf("pack result  final : ");
	 dumpSimdData(&result_FFT4[0], 0);
	 */
	 
	 
	 _mm_empty();
	 _m_empty();
	 
	 return((int16_t*) result_FFT4);


}


int16_t* bflyRadix4_16(int16_t *x, int16_t *y, int16_t *minus_W0, int16_t *minus_W1, int16_t *minus_W2, int16_t *minus_W3, 
						int16_t *shuffle_W0, int16_t *shuffle_W1, int16_t *shuffle_W2, int16_t *shuffle_W3)
{
	__m128i *pW16, *x128, *y128;
	__m128i mmtmpre0, mmtmpre1, mmtmpre2, mmtmpre3, mmtmpim0, mmtmpim1, mmtmpim2, mmtmpim3, mmre, mmim;

	
	y128 = (__m128i*)y;


	for(int i=0;i<4;i++)
	{
		
		x128 = (__m128i*)x;
		
		//1st term's (i)th output 4개
		pW16 = (__m128i*)&minus_W0[i*8]; 
		mmtmpre0 = _mm_madd_epi16(*x128, pW16[0]);	
				
		pW16= (__m128i*)&shuffle_W0[i*8];
		mmtmpim0 = _mm_madd_epi16(*x128, pW16[0]);
			
		mmtmpre0 = _mm_srai_epi32(mmtmpre0, 15);
		mmtmpim0 = _mm_srai_epi32(mmtmpim0, 15);

		x128++;
		
		
		//2nd term's (i)th output 4개
		pW16 = (__m128i*)&minus_W1[i*8];

		
		mmtmpre1= _mm_madd_epi16(*x128, pW16[0]);
		
		pW16 = (__m128i*)&shuffle_W1[i*8];

		
		mmtmpim1 = _mm_madd_epi16(*x128, pW16[0]);


		mmtmpre1 = _mm_srai_epi32(mmtmpre1, 15);
		mmtmpim1 = _mm_srai_epi32(mmtmpim1, 15);


		x128++;


		//3rd term's (i)th output 4개
		pW16 = (__m128i*)&minus_W2[i*8];
		mmtmpre2= _mm_madd_epi16(*x128, pW16[0]);

		
		pW16 = (__m128i*)&shuffle_W2[i*8];
		mmtmpim2 = _mm_madd_epi16(*x128, pW16[0]);


		mmtmpre2 = _mm_srai_epi32(mmtmpre2, 15);
		mmtmpim2 = _mm_srai_epi32(mmtmpim2, 15);

			
		x128++;

		
		//4th term's (i)th output 4개
		pW16 = (__m128i*)&minus_W3[i*8];
		mmtmpre3= _mm_madd_epi16(*x128, pW16[0]);
		
		pW16 = (__m128i*)&shuffle_W3[i*8];
		mmtmpim3 = _mm_madd_epi16(*x128, pW16[0]);
		

		mmtmpre3 = _mm_srai_epi32(mmtmpre3, 15);
		mmtmpim3 = _mm_srai_epi32(mmtmpim3, 15);

	
		//final result_ cumulation
		mmre = _mm_add_epi32(mmtmpre0, mmtmpre1);
		mmre = _mm_add_epi32(mmre, mmtmpre2);
		mmre = _mm_add_epi32(mmre, mmtmpre3); //real 4개가 나와야. 루프 돌 때마다 그 다음 real이 나와야. 

		mmim = _mm_add_epi32(mmtmpim0, mmtmpim1);
		mmim = _mm_add_epi32(mmim
, mmtmpim2);
		mmim = _mm_add_epi32(mmim, mmtmpim3); //im 4개가 나와야. i가 증가할 때마다 그다음 4개 imaginary수가 나와야.
		
		
		//bit control
		mmre = _mm_packs_epi32(mmre, mmre);
		mmim = _mm_packs_epi32(mmim, mmim);	


		y128[0] = _mm_unpacklo_epi16(mmre, mmim);

		y128++;

	}

				 
	return(y);

}


	
int16_t* bflyRadix4_Nbit(int16_t *x, int16_t *y, int16_t *minus_W0, int16_t *minus_W1, int16_t *minus_W2, int16_t *minus_W3, 
						int16_t *shuffle_W0, int16_t *shuffle_W1, int16_t *shuffle_W2, int16_t *shuffle_W3, int N)
{
	__m128i *pW_N, *x128, *y128;
	__m128i __attribute__((aligned(16))) mmtmpre0[10000], __attribute__((aligned(16))) mmtmpre1[10000], __attribute__((aligned(16))) mmtmpre2[10000], __attribute__((aligned(16))) mmtmpre3[10000], 
		__attribute__((aligned(16))) mmtmpim0[10000], __attribute__((aligned(16))) mmtmpim1[10000], __attribute__((aligned(16))) mmtmpim2[10000], __attribute__((aligned(16))) mmtmpim3[10000], mmre, mmim;
	
	
	y128 = (__m128i*)y;


	for(int i=0; i<4; i++)
	{
		
		x128 = (__m128i*)x;

		//1st term's (i)th output 4개
		for(int j=0; j<(N/16); j++) //-----------------------add
		{
			pW_N = (__m128i*)&minus_W0[(N/2)*i+8*j]; 				//if i=0, minus_W0[0], minus_W0[8], minus_W0[16], minus_W0[24] 	--re	
			
			mmtmpre0[j] = _mm_madd_epi16(*x128, pW_N[0]);

			pW_N= (__m128i*)&shuffle_W0[(N/2)*i+8*j];             //if i=0, shuffle_W0[0], shuffle_W0[8], shuffle_W0[16], shuffle_W0[24] 	--im
			mmtmpim0[j] = _mm_madd_epi16(*x128, pW_N[0]);

			mmtmpre0[j] = _mm_srai_epi32(mmtmpre0[j], 15);
			mmtmpim0[j] = _mm_srai_epi32(mmtmpim0[j], 15);

			pW_N++;
			x128++;
		}

		
		
		//2nd term's (i)th output 4개 
		for(int j=0; j<(N/16); j++) //-----------------------add
		{
			pW_N = (__m128i*)&minus_W1[(N/2)*i+8*j]; 				//if i=0, minus_W0[0], minus_W0[8], minus_W0[16], minus_W0[24] 	--re	
			
			mmtmpre1[j] = _mm_madd_epi16(*x128, pW_N[0]);

			pW_N= (__m128i*)&shuffle_W1[(N/2)*i+8*j];             //if i=0, shuffle_W0[0], shuffle_W0[8], shuffle_W0[16], shuffle_W0[24] 	--im
			mmtmpim1[j] = _mm_madd_epi16(*x128, pW_N[0]);

			mmtmpre1[j] = _mm_srai_epi32(mmtmpre1[j], 15);
			mmtmpim1[j] = _mm_srai_epi32(mmtmpim1[j], 15);
 
			pW_N++;
			x128++;
		}



		//3rd term's (i)th output 4개
		for(int j=0; j<(N/16); j++) //-----------------------add
		{
			pW_N = (__m128i*)&minus_W2[(N/2)*i+8*j]; 				//if i=0, minus_W0[0], minus_W0[8], minus_W0[16], minus_W0[24] 	--re	
			
			mmtmpre2[j] = _mm_madd_epi16(*x128, pW_N[0]);

			pW_N= (__m128i*)&shuffle_W2[(N/2)*i+8*j];             //if i=0, shuffle_W0[0], shuffle_W0[8], shuffle_W0[16], shuffle_W0[24] 	--im
			mmtmpim2[j] = _mm_madd_epi16(*x128, pW_N[0]);

			mmtmpre2[j] = _mm_srai_epi32(mmtmpre2[j], 15);
			mmtmpim2[j] = _mm_srai_epi32(mmtmpim2[j], 15);

			pW_N++;
			x128++;
		}


		
		//4th term's (i)th output 4개
		for(int j=0; j<(N/16); j++) //-----------------------add
		{
			pW_N = (__m128i*)&minus_W3[(N/2)*i+8*j]; 				//if i=0, minus_W0[0], minus_W0[8], minus_W0[16], minus_W0[24] 	--re	
			
			mmtmpre3[j] = _mm_madd_epi16(*x128, pW_N[0]);

			pW_N= (__m128i*)&shuffle_W3[(N/2)*i+8*j];             //if i=0, shuffle_W0[0], shuffle_W0[8], shuffle_W0[16], shuffle_W0[24] 	--im
			mmtmpim3[j] = _mm_madd_epi16(*x128, pW_N[0]);

			mmtmpre3[j] = _mm_srai_epi32(mmtmpre3[j], 15);
			mmtmpim3[j] = _mm_srai_epi32(mmtmpim3[j], 15);

			pW_N++;
			x128++;
		}




		//final result_ cumulation
		for(int j=0; j<(N/16); j++)
		{

			mmre = _mm_add_epi32(mmtmpre0[j], mmtmpre1[j]);
			mmre = _mm_add_epi32(mmre, mmtmpre2[j]);
			mmre = _mm_add_epi32(mmre, mmtmpre3[j]); 

			mmim = _mm_add_epi32(mmtmpim0[j], mmtmpim1[j]);
			mmim = _mm_add_epi32(mmim
, mmtmpim2[j]);
			mmim = _mm_add_epi32(mmim, mmtmpim3[j]); 
		
			
			//bit control
			mmre = _mm_packs_epi32(mmre, mmre);
			mmim = _mm_packs_epi32(mmim, mmim);	

			y128[0] = _mm_unpacklo_epi16(mmre, mmim);


			y128++;
		}
		
	}

				 
	return(y);
	
}




int16_t* FFT16(int16_t *x, int16_t *y) //ssss
{
	int16_t __attribute__((aligned(16))) y_4[10000]; //*outputy;
	int32_t *x32;
	__m128i x128, *y128_4;

	x32 = (int32_t*)x;
	y128_4 = (__m128i*)y_4;

	for(int i=0;i<4;i++)
	{
		
		x128 =_mm_set_epi32(x32[i+12], x32[i+8], x32[i+4], x32[i]); 
		FFT4((int16_t*)&x128, (int16_t*)y128_4);
		y128_4++;
	}
	

	bflyRadix4_16(y_4, y, minus_W16_0, minus_W16_1, minus_W16_2, minus_W16_3, shuffle_W16_0, shuffle_W16_1, shuffle_W16_2, shuffle_W16_3);


	_mm_empty();
	_m_empty();

	return (y);
}



int16_t* FFT64(int16_t *x, int16_t *y) 
	{
		int16_t __attribute__((aligned(16))) y_16[10000]; 
		int32_t *x32;
		int16_t __attribute__((aligned(16)))x128[10000];
		__m128i *px128, *y128_16;
	
		x32 = (int32_t*)x;
		y128_16 = (__m128i*)y_16;
		px128 = (__m128i*)x128;
		
		for(int i=0;i<4;i++)
			{
				for(int j=0; j<4; j++)
					{		
						//printf("FFT64 : (%i, %i)\n",i,j);
						//fflush(stdout);
						px128[0] =_mm_set_epi32(x32[i+16*j+12], x32[i+16*j+8], x32[i+16*j+4], x32[i+16*j]); //if i=0, x(0), x(4), x(8), x(12), ... , x(60) , N=64's 'input
						px128++;	
					}
		
				FFT16((int16_t*)(px128-4), (int16_t*)y128_16);
				y128_16 = y128_16+4;
		
			}

		bflyRadix4_Nbit(y_16, y,  minus_W64_0, minus_W64_1, minus_W64_2, minus_W64_3, shuffle_W64_0, shuffle_W64_1, shuffle_W64_2, shuffle_W64_3, 64);

		
		_mm_empty();
		_m_empty();

		return (y);
	
	}



int16_t* FFT256(int16_t *x, int16_t *y) //ssss
	{
		int16_t __attribute__((aligned(16))) y_64[10000]; //*outputy;
		int32_t *x32;
		int16_t __attribute__((aligned(16))) x128[10000];
		__m128i *px128, *y128_64;

		x32 = (int32_t*)x;
		px128 = (__m128i*)x128;
		y128_64 = (__m128i*)y_64;

		for(int i=0;i<4;i++)
			{
				for(int j=0; j<16; j++)
					{		
						px128[0] =_mm_set_epi32(x32[i+16*j+12], x32[i+16*j+8], x32[i+16*j+4], x32[i+16*j]); //if i=0, x(0), x(4), x(8), x(12), ... , x(60) , N=64's 'input
						px128++;
						
					}
				
					
				FFT64((int16_t*)(px128-16), (int16_t*)y128_64);
				
				y128_64 = y128_64+16;
				

			}

		bflyRadix4_Nbit(y_64, y,  minus_W256_0, minus_W256_1, minus_W256_2, minus_W256_3, shuffle_W256_0, shuffle_W256_1, shuffle_W256_2, shuffle_W256_3, 256);
		
		_mm_empty();
		_m_empty();

		return (y);
		


	
	}


int16_t* FFT1024(int16_t *x, int16_t *y) //ssss
	{
		int16_t __attribute__((aligned(16))) y_256[10000]; //*outputy;
		int32_t *x32;
		int16_t __attribute__((aligned(16))) x128[10000];
		__m128i *px128, *y128_256;

		x32 = (int32_t*)x;
		px128 = (__m128i*)x128;
		y128_256 = (__m128i*)y_256;

		for(int i=0;i<4;i++)
			{
				for(int j=0; j<64; j++)
					{		
						px128[0] =_mm_set_epi32(x32[i+16*j+12], x32[i+16*j+8], x32[i+16*j+4], x32[i+16*j]); //if i=0, x(0), x(4), x(8), x(12), ... , x(60) , N=64's 'input
						px128++;
						
					}
				
					
				FFT256((int16_t*)(px128-64), (int16_t*)y128_256);
				
				y128_256 = y128_256+64;

			}

		bflyRadix4_Nbit(y_256, y,  minus_W1024_0, minus_W1024_1, minus_W1024_2, minus_W1024_3, shuffle_W1024_0, shuffle_W1024_1, shuffle_W1024_2, shuffle_W1024_3, 1024);

		_mm_empty();
		_m_empty();


		return (y);


	
	}

	