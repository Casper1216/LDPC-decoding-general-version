#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<stdbool.h>
#include<iostream>
#include<vector>
#include<unordered_map>
#define pi acos(-1)

using namespace std;

//fixed execution time: numtime
//add
//vector<unordered_map<int,int>> VN_2_CN(n);	//VN i mapping: key: connect CN j ,value: index in CN_set 
//vector<unordered_map<int,int>> CN_2_VN(m);	//CN j mapping: key: connect VN i ,value: index in VN_set

void LDPC_SPA(double* BER,double* FER,int n,int m,int dv,int dc,double R,int** CN_set,int** VN_set,int* row,int* col,vector<unordered_map<int,int>>& VN_2_CN,vector<unordered_map<int,int>>& CN_2_VN,
double* avgIter,double* SNR_dB,int SNR_L,int iteration,int numtime);


int main(){
	srand(time(0));
	//LDPC
	//BPSK 1 or -1
	
	FILE *fp1 = fopen("H_36736_3584.txt", "r");
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
    
    for(int i=0;i<n;i++)
    	fscanf(fp1,"%d ",&col[i]);	
	
	int* row = (int*)malloc(sizeof(int)*m);		//number of row nonzero 
    
    for(int i=0;i<m;i++)
    	fscanf(fp1,"%d ",&row[i]);	
	
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
	vector<unordered_map<int,int>> VN_2_CN(n);	//VN i mapping: key: connect CN j ,value: index in CN_set 
	vector<unordered_map<int,int>> CN_2_VN(m);	//CN j mapping: key: connect VN i ,value: index in VN_set 
	for(int i=0;i<n;i++){
		for(int j=0;j<col[i];j++){
			fscanf(fp1,"%d ",&VN_set[i][j]);
			VN_set[i][j]--;
			VN_2_CN[i][VN_set[i][j]] = j;
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
			CN_2_VN[i][CN_set[i][j]] = j;
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

	
	//***************************************************************************
	//平均 iteration 
	double *avgIter = (double *)malloc(sizeof(double) * SNR_L);
	
	int numtime =1;
	int iteration = 50;	
	
	double *BER = (double*)malloc(sizeof(double) * SNR_L);
	double *FER = (double*)malloc(sizeof(double) * SNR_L);

	// Start Record the time
    time_t  start = clock();

	LDPC_SPA(BER,FER, n, m, dv, dc, R,CN_set, VN_set, row, col, VN_2_CN, CN_2_VN, avgIter, SNR_dB, SNR_L, iteration, numtime);

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

void LDPC_SPA(double* BER,double* FER,int n,int m,int dv,int dc,double R,int** CN_set,int** VN_set,int* row,int* col,vector<unordered_map<int,int>>& VN_2_CN,vector<unordered_map<int,int>>& CN_2_VN,
double* avgIter,double* SNR_dB,int SNR_L,int iteration,int numtime){
	//channel information
	double* Fn = (double *)malloc(sizeof(double) * n);
	//tau
	double** tau= (double **)malloc(sizeof(double*) * m);
	for(int i=0;i<m;i++){
		tau[i] = (double *)malloc(sizeof(double) * dc);
	}
	
	//CN update //CN[j][i]  CN j to VN i
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
		int count = 1;
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
					for(int i=0;i<row[j];i++){				//CN_j to VN_{CN_set[j][i]}
	
						//printf("CN:%d to VN:%d\n",j,CN_set[j][i]);
						for(int np=0 ; np<row[j] ; np++){			//n'
							if(CN_set[j][i]!=CN_set[j][np]){		//n != n'
								tau[j][i] *=tanh(VN[CN_set[j][np]][VN_2_CN[CN_set[j][np]][j]]/2);	
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

				//VN update
	
				for(int i=0;i<n;i++){ 				//go through all VN	
					for(int j=0;j<col[i];j++){		//VN_i to CN_{VN_set[i][j]} 
						
						VN[i][j] = Fn[i];	
						//printf("VN:%d to CN:%d\n",i,VN_set[i][j]);
						for(int mp=0 ; mp<col[i] ; mp++){			//m'
							if(VN_set[i][mp]>=0&&VN_set[i][j]!=VN_set[i][mp]){	//m != m'
								
								VN[i][j] += CN[VN_set[i][mp]][CN_2_VN[VN_set[i][mp]][i]];	
								
										
							}		
						}	
			
					}
					
				}

				
				
				
				//total LLR
				//decode
					
				for(int i=0;i<n;i++){ 				//go through all VN	
					VN_total[i] =  Fn[i]; 			
					for(int j=0;j<col[i];j++){		

						VN_total[i] += CN[VN_set[i][j]][CN_2_VN[VN_set[i][j]][i]];
	
					}
					//printf("VN_toatl %f\n",VN_total[i]);
					
					if(VN_total[i]>=0)
						u_hat[i] = 0;
					else
						u_hat[i] = 1;
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
			//printf("%d\n",iter);
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
			
			if(count==10000){
				printf("error: %lld, num: %lld, BER: %E, FER: %E Average iteration: %f\n",error,num,((double)error)/((double)(n*num)),((double)frameerror)/(num),(double)avgIter[q]/num);
				count=0;
			}
			count++; 
			
		}
		
		BER[q] = ((double)error)/(double)n/(double)numtime;
		FER[q] = ((double)frameerror)/((double)numtime);
		printf("BER: %E, FER: %E Average iteration: %f\n",BER[q],FER[q],avgIter[q]/numtime);
	
	} 
}






