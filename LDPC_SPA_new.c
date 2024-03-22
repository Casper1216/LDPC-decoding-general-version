#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<stdbool.h>
#define pi acos(-1)

//fixed execution time: numtime
//Reduce VN update complexity

void LDPC_SPA(double* BER,double* FER,int n,int m,int dv,int dc,double R,int** CN_set,int** VN_set,int* row,int* col,double* avgIter,double* SNR_dB,int SNR_L
,int iteration, int numtime);

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
	SNR_dB[0] = 0.4;
	//SNR_dB[1] = 3;
	//SNR_dB[2] = 2;
	// SNR_dB[3] = 2.5;
	// SNR_dB[4] = 3;
	
	//***************************************************************************
	//平均 iteration 
	double *avgIter = (double *)malloc(sizeof(double) * SNR_L);
	
	int numtime =1;
	int iteration = 50;	
	
	double *BER = (double*)malloc(sizeof(double) * SNR_L);
	double *FER = (double*)malloc(sizeof(double) * SNR_L);

	// Start Record the time
    time_t  start = clock();
	
	LDPC_SPA(BER,FER, n, m, dv, dc, R,CN_set, VN_set, row, col, avgIter, SNR_dB, SNR_L, iteration, numtime);

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
			//else
				//fprintf(fp, "%E,",avgIter[j]/numtime);
				
		}	
		
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	
	return 0;
}

void LDPC_SPA(double* BER,double* FER,int n,int m,int dv,int dc,double R,int** CN_set,int** VN_set,int* row,int* col,double* avgIter,double* SNR_dB,int SNR_L
,int iteration,int numtime){
	//channel information
	double* Fn = (double *)malloc(sizeof(double) * n);
	//tau
	double** tau= (double **)malloc(sizeof(double*) * m);
	for(int i=0;i<m;i++){
		tau[i] = (double *)malloc(sizeof(double) * dc);
	}
	
	//CN update //CN[j][i]  CN j to CN i
	double** CN = (double **)malloc(sizeof(double*) * m);
	for(int i=0;i<m;i++){
		CN[i] = (double *)malloc(sizeof(double) * dc);
	}
	
	//VN update //VN[i][j]  VN i to CN j
	double** VN = (double **)malloc(sizeof(double*) * n);
	for(int i=0;i<n;i++){
		VN[i] = (double *)malloc(sizeof(double) * dv);
	}
	
	//total LLR VN[i]
	double* VN_total = (double *)malloc(sizeof(double) * n);
	
	//decoded sequence
	int* u_hat = (int *)malloc(sizeof(int) * n);
	
	double *y = (double*)malloc(sizeof(double) * n);
	//noise
	double* noise = (double *)malloc(sizeof(double) * n);
	
	//判斷symdrome
	int* s = (int *)malloc(sizeof(int) * m);
	
	//sigal;
	//假設全為 0
	int *u = (int*)malloc(sizeof(int) * n);	//binary sequence
	int *x = (int*)malloc(sizeof(int) * n);
	
	for(int i=0;i<n;i++){
		u[i] = 0;
		x[i] = 1;		//Eb=1	Eav = R*1
	} 
	
	for(int q=0;q<SNR_L;q++){
		
		avgIter[q] = 0;
		long long error=0;
		long long frameerror=0;
		for(long long num=0;num<numtime;num++){
			
			
			double U ,V ;
			double sigma = sqrt((1/(2*R))*pow(10,-(SNR_dB[q]/10)));
			
			
			for(int i=0;i<n;i++){
				U = (double) rand() / (double)(RAND_MAX);	//uniform (0~1) 0<U<1
				while(U==0||U==1)
					U = (double) rand() / (double)(RAND_MAX);	
				
				
				V = (double) rand() / (double)(RAND_MAX);
				while(V==0||V==1)
					V = (double) rand() / (double)(RAND_MAX);	
				
				noise[i] = sqrt(-2*log(U))*cos(2*pi*V);
				noise[i] = sigma* noise[i] ;
				
				
				y[i] = x[i] + noise[i];
				
			} 
			
			//initialization
			for(int i=0;i<n;i++){
				Fn[i] = 2*y[i]/(pow(sigma,2));	
				VN_total[i] =  Fn[i]; 	
				for(int j=0;j<col[i];j++){	
					VN[i][j] = 2*y[i]/(pow(sigma,2));	
				}	
				
			}
			
			//iterative decoding
			int iter=0;
			for(iter=0;iter<iteration;iter++){
				

				//CN update	
				for(int j=0;j<m;j++){	
					for(int i=0;i<row[j];i++){
						tau[j][i]=1;	
					}
				}
			
				
				for(int j=0;j<m;j++){						//go through all CN	
					for(int i=0;i<row[j];i++){				//相連之 VN 

						//printf("CN:%d to VN:%d\n",j,CN_set[j][i]);
						for(int np=0 ; np<row[j] ; np++){			//n'
							if(CN_set[j][np]>=0&&CN_set[j][i]!=CN_set[j][np]){	//n != n'
								
								//printf("VN:%d to CN:%d\n",CN_set[j][np],j); //N(m) set
								for(int f=0;f<col[CN_set[j][np]];f++){
									//找到VN 中之index 
									if(VN_set[CN_set[j][np]][f]==j){
										
										tau[j][i] *=tanh(VN[CN_set[j][np]][f]/2);	
										
									}
								}
							}		
						}
						
						//計算完tau 
						//CN[j][i] = log((1.0+tau[j][i])/(1.0-tau[j][i]));
						if(tau[j][i]==1)
							CN[j][i] = DBL_MAX;
						else if(tau[j][i]==-1)
							CN[j][i] = -DBL_MAX;
						else
							CN[j][i] = 2*atanh(tau[j][i]);
						
					}
				}
				

				//total LLR
				//decode
				for(int i=0;i<n;i++){ 				//go through all VN	
					VN_total[i] =  Fn[i]; 			
					for(int j=0;j<col[i];j++){				
						for(int m=0;m<row[VN_set[i][j]];m++){
							if(CN_set[VN_set[i][j]][m]==i){		//找到與VN相連之 CN  idex 
								//printf("VN %d CN %d to %d\n",i,VN_set[i][j],CN_set[VN_set[i][j]][m]);
								VN_total[i] += CN[VN_set[i][j]][m];
								
							}
						}
	
					}
					//printf("VN_toatl %f\n",VN_total[i]);
					
					if(VN_total[i]>=0)
						u_hat[i] = 0;
					else
						u_hat[i] = 1;
				}
				//VN update
				for(int i=0;i<n;i++){ 				//go through all VN	
					for(int j=0;j<col[i];j++){			//相連之 CN 
						
						//printf("VN:%d to CN:%d\n",i,VN_set[i][j]);

						for(int p=0;p<row[VN_set[i][j]];p++){
							//找到 CN 中之index 
							if(CN_set[VN_set[i][j]][p]==i){
			
								VN[i][j] = VN_total[i] - CN[VN_set[i][j]][p];		
							}
						}
					}
				}
				
				bool iszerovector = true;
				
				for(int j=0;j<m;j++){
					s[j] =0;
					for(int i=0;i<row[j];i++){
						if(CN_set[j][i]>=0&&u_hat[CN_set[j][i]]==1)
							s[j] += 1;
					}
					
					s[j] %=2;
					
					if(s[j]!=0){
						iszerovector = false;
						break;
					}	
				}
				
				if(iszerovector){
					break;
				}
				
			}
			avgIter[q] +=iter;
			
			for(int i=0;i<n;i++){	
				if(u[i]!=u_hat[i])
					error++;
			}
			
			for(int i=0;i<n;i++){
				if(u[i]!=u_hat[i]){
					frameerror++;
					break;
				}
					
			}
			
			printf("error: %lld, num: %lld, BER: %E, FER: %E Average iteration: %f\r",error,num,((double)error)/((double)(n*num)),((double)frameerror)/(num),(double)avgIter[q]/num);			
		}
		
		BER[q] = ((double)error)/(double)n/(double)numtime;
		FER[q] = ((double)frameerror)/((double)numtime);
		printf("BER: %E, FER: %E Average iteration: %f\n",BER[q],FER[q],avgIter[q]/numtime);
	
	} 
}






