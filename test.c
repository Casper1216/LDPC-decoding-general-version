#include "LDPC_decoder.h"
#ifndef _STDIO_H_
#define _STDIO_H_
#endif
#ifndef _STDLIB_H_
#define _STDLIB_H_
#endif

int main(){
	srand(time(0));
	//LDPC
	//BPSK 1 or -1
	
	FILE *fp1 = fopen("H_1944_972.txt", "r");
	if (fp1 == NULL) {
        fprintf(stderr, "fopen() failed.\n");
        exit(EXIT_FAILURE);
    }
    
    int n,m,dv,dc;
    
    
    fscanf(fp1,"%d ",&n);
    fscanf(fp1,"%d ",&m);
    fscanf(fp1,"%d ",&dv);
    fscanf(fp1,"%d ",&dc);
    
    printf("%d %d %d %d\n",n,m,dv,dc);
    
    const double R=(double)(n-m)/n; //coderate
    
    //建立 number of col row nonzero location
	 
    int* col = (int*)malloc(sizeof(int)*n);		//number of col nonzero 
    
    for(int i=0;i<n;i++){
    	fscanf(fp1,"%d ",&col[i]);
    	
	}
	
	int* row = (int*)malloc(sizeof(int)*m);		//number of row nonzero 
    
    for(int i=0;i<m;i++){
    	fscanf(fp1,"%d ",&row[i]);	
	}	
	
	
	//Construct adjacency list
	//VN 
	int** VN_set = (int **)malloc(sizeof(int*) * n);
	for(int i=0;i<n;i++){
		VN_set[i] = (int *)malloc(sizeof(int) * col[i]);
	}
	//CN
	int** CN_set = (int **)malloc(sizeof(int*) * m);
	for(int i=0;i<m;i++){
		CN_set[i] = (int *)malloc(sizeof(int) * row[i]);
	}	

	
	for(int i=0;i<n;i++){
		for(int j=0;j<col[i];j++){
			fscanf(fp1,"%d ",&VN_set[i][j]);
			VN_set[i][j]--;
		}
		for(int k=col[i];k<dv;k++){
			int buffer;
			fscanf(fp1,"%d ",&buffer);	//把多餘的 col idex 用 buffer 消耗掉 
		}
	}
	
	
	for(int i=0;i<m;i++){
		for(int j=0;j<row[i];j++){
			fscanf(fp1,"%d ",&CN_set[i][j]);
			CN_set[i][j]--;
		}
		for(int k=row[i];k<dc;k++){
			int buffer;
			fscanf(fp1,"%d ",&buffer);	//把多餘的 col idex 用 buffer 消耗掉 
		}
	}
	
	fclose(fp1);
	//***************************************************************************
	//SNR
	const int SNR_L = 1;
	double *SNR_dB = (double *)malloc(sizeof(double) * SNR_L);
	SNR_dB[0] = 1.2;

	//**********************************************
	//平均 iteration 
	double *avgIter = (double *)malloc(sizeof(double) * SNR_L);
	
    int numtime = 10000;
	int frameerror 	= 100;
	int iteration 	= 50;	
	
	double *BER = (double*)malloc(sizeof(double) * SNR_L);
	double *FER = (double*)malloc(sizeof(double) * SNR_L);
	// Start Record the time
    time_t  start = clock();
	

    //***************************************************************************
	//LDPC_SPA_until_frame(BER,FER, n, m, dv, dc, R,CN_set, VN_set, row, col, avgIter, SNR_dB, SNR_L, iteration, frameerror);
    //LDPC_SPA(BER,FER, n, m, dv, dc, R,CN_set, VN_set, row, col, avgIter, SNR_dB, SNR_L, iteration, numtime);
    //LDPC_MSA(BER,FER, n, m, dv, dc, R,CN_set, VN_set, row, col, avgIter, SNR_dB, SNR_L, iteration, numtime);
    //LDPC_LBP(BER,FER, n, m, dv, dc, R,CN_set, VN_set, row, col, avgIter, SNR_dB, SNR_L, iteration, numtime);
    double normalize_factor = 0.75;
    LDPC_Layered_NMSA(BER,FER, n, m, dv, dc, R,CN_set, VN_set, row, col, avgIter, SNR_dB, SNR_L, iteration, numtime,normalize_factor);


    //***************************************************************************


	// Record the end time
    time_t end = clock();

    double diff = end - start; // ms
    printf(" %f  sec\n", diff / CLOCKS_PER_SEC );
    
	
	
	//寫入檔案 CSV
	FILE *fp = fopen("LDPC_SPA.csv", "w");
    
    //避免開啟失敗 
    if (fp == NULL) {
        fprintf(stderr, "fopen() failed.\n");
        exit(EXIT_FAILURE);
    }
    
    
    for(int i=0;i<4;i++){
    	//printf("%d %d\n",q++,H[i][0]);
		for(int j=0;j<SNR_L;j++){
			if(i==0)
				fprintf(fp, "%f,",SNR_dB[j]);
			else if(i==1)
				fprintf(fp, "%E,",BER[j]);
			else if(i==2)
				fprintf(fp, "%E,",FER[j]);
				
		}	
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	
	return 0;
}