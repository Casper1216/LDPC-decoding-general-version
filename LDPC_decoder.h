#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<float.h>
#include<limits.h>
#include<stdbool.h>
#define pi acos(-1)

void Add_AWGN(int n,double* y,int* x,double sigma){
    double U ,V ;
    for(int i=0;i<n;i++){
        U = (double) rand() / (double)(RAND_MAX);	//uniform (0~1) 0<U<1
        while(U==0||U==1)
            U = (double) rand() / (double)(RAND_MAX);	
        V = (double) rand() / (double)(RAND_MAX);
        while(V==0||V==1)
            V = (double) rand() / (double)(RAND_MAX);	
        //AWGN(0,1) =  sqrt(-2*log(U))*cos(2*pi*V)
        y[i] = (double)x[i] + sigma* sqrt(-2*log(U))*cos(2*pi*V);
    } 
}

//Run until the number of execution time == numtime
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
    //receive signal
	double *y = (double*)malloc(sizeof(double) * n);
	//symdrome
	int* s = (int *)malloc(sizeof(int) * m);
	//transmitted sigal;
	//All 0 codeword
	int *u = (int*)malloc(sizeof(int) * n);	//binary sequence
	int *x = (int*)malloc(sizeof(int) * n);
	
	for(int i=0;i<n;i++){
		u[i] = 0;
		x[i] = 1;		//Eb=1	Eavg = R*1
	} 
	
	for(int q=0;q<SNR_L;q++){
		avgIter[q] = 0;
		long long error=0;
		int frameerror=0;

		for(int num=0;num<numtime;num++){
			double sigma = sqrt((1/(2*R))*pow(10,-(SNR_dB[q]/10)));
            Add_AWGN(n,y,x,sigma);
			
			//initialization
			for(int i=0;i<n;i++){
				Fn[i] = 2*y[i]/(pow(sigma,2));	
				VN_total[i] =  Fn[i]; 	
				for(int j=0;j<col[i];j++)
					VN[i][j] = 2*y[i]/(pow(sigma,2));	
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
				for(int j=0;j<m;j++){						//go through all CN	j
					for(int i=0;i<row[j];i++){				//VN_{CN_set[j][i]} connected by CN j
						for(int np=0 ; np<row[j] ; np++){			//n'
							if(CN_set[j][np]>=0&&CN_set[j][i]!=CN_set[j][np]){	//n != n'
								
								for(int f=0;f<col[CN_set[j][np]];f++){
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
							if(CN_set[VN_set[i][j]][m]==i){		
								//找到與VN相連之 CN  idex 
								VN_total[i] += CN[VN_set[i][j]][m];
								
							}
						}
					}					
					if(VN_total[i]>=0)
						u_hat[i] = 0;
					else
						u_hat[i] = 1;
				}
				//VN update
				for(int i=0;i<n;i++){ 					//go through all VN	i
					for(int j=0;j<col[i];j++){			//CN_{VN_set[j][i]} connected by VN i
						for(int p=0;p<row[VN_set[i][j]];p++){
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
			printf("error: %lld, num: %d, BER: %E, FER: %E Average iteration: %f\r",error,num,((double)error)/((double)(n*num)),((double)frameerror)/(num),(double)avgIter[q]/num);
		}
		
		BER[q] = ((double)error)/(double)n/(double)numtime;
		FER[q] = ((double)frameerror)/((double)numtime);
		printf("BER: %E, FER: %E Average iteration: %f\n",BER[q],FER[q],avgIter[q]/numtime);
	} 
}
// void LDPC_SPA(double* BER,double* FER,int n,int m,int dv,int dc,double R,int** CN_set,int** VN_set,int* row,int* col,double* avgIter,double* SNR_dB,int SNR_L
// ,int iteration,int numtime){
// 	//channel information
// 	double* Fn = (double *)malloc(sizeof(double) * n);
// 	//tau
// 	double** tau= (double **)malloc(sizeof(double*) * m);
// 	for(int i=0;i<m;i++){
// 		tau[i] = (double *)malloc(sizeof(double) * dc);
// 	}
// 	//CN update //CN[j][i]  CN j to CN i
// 	double** CN = (double **)malloc(sizeof(double*) * m);
// 	for(int i=0;i<m;i++){
// 		CN[i] = (double *)malloc(sizeof(double) * dc);
// 	}
// 	//VN update //VN[i][j]  VN i to CN j
// 	double** VN = (double **)malloc(sizeof(double*) * n);
// 	for(int i=0;i<n;i++){
// 		VN[i] = (double *)malloc(sizeof(double) * dv);
// 	}
// 	//total LLR VN[i]
// 	double* VN_total = (double *)malloc(sizeof(double) * n);
// 	//decoded sequence
// 	int* u_hat = (int *)malloc(sizeof(int) * n);
//     //receive signal
// 	double *y = (double*)malloc(sizeof(double) * n);
// 	//symdrome
// 	int* s = (int *)malloc(sizeof(int) * m);
// 	//transmitted sigal;
// 	//All 0 codeword
// 	int *u = (int*)malloc(sizeof(int) * n);	//binary sequence
// 	int *x = (int*)malloc(sizeof(int) * n);
	
// 	for(int i=0;i<n;i++){
// 		u[i] = 0;
// 		x[i] = 1;		//Eb=1	Eavg = R*1
// 	} 
	
// 	for(int q=0;q<SNR_L;q++){
// 		avgIter[q] = 0;
// 		long long error=0;
// 		int frameerror=0;

// 		for(int num=0;num<numtime;num++){
// 			double sigma = sqrt((1/(2*R))*pow(10,-(SNR_dB[q]/10)));
//             Add_AWGN(n,y,x,sigma);
			
// 			//initialization
// 			for(int i=0;i<n;i++){
// 				Fn[i] = 2*y[i]/(pow(sigma,2));	
// 				for(int j=0;j<col[i];j++){	
					
// 					VN[i][j] = 2*y[i]/(pow(sigma,2));	
// 				}	
// 			}
// 			//iterative decoding
// 			int iter=0;
// 			for(iter=0;iter<iteration;iter++){
				
// 				//CN update	
// 				for(int j=0;j<m;j++){	
// 					for(int i=0;i<row[j];i++){
// 						tau[j][i]=1;	
// 					}
// 				}
				
// 				for(int j=0;j<m;j++){						//go through all CN	j
// 					for(int i=0;i<row[j];i++){				//VN i connected by CN j

// 							for(int np=0 ; np<row[j] ; np++){			//n'
// 								if(CN_set[j][np]>=0&&CN_set[j][i]!=CN_set[j][np]){	//n != n'
									
// 									for(int f=0;f<col[CN_set[j][np]];f++){
// 										//找到VN 中之index 
// 										if(VN_set[CN_set[j][np]][f]==j){
											
// 											tau[j][i] *=tanh(VN[CN_set[j][np]][f]/2);	
// 										}
// 									}
// 								}		
// 							}
// 						//tau is computed.
// 						//CN[j][i] = log((1.0+tau[j][i])/(1.0-tau[j][i]));
// 						if(tau[j][i]==1)
// 							CN[j][i] = DBL_MAX;
// 						else if(tau[j][i]==-1)
// 							CN[j][i] = -DBL_MAX;
// 						else
// 							CN[j][i] = 2*atanh(tau[j][i]);
// 					}
// 				}
				
// 				//VN update
// 				for(int i=0;i<n;i++){ 				    //go through all VN	i
// 					for(int j=0;j<col[i];j++){			//CN j connected by VN i
						
//                         VN[i][j] = Fn[i];	
//                         for(int mp=0 ; mp<col[i] ; mp++){			//m'
//                             if(VN_set[i][mp]>=0&&VN_set[i][j]!=VN_set[i][mp]){	//m != m'
                                
//                                 for(int p=0;p<row[VN_set[i][mp]];p++){
//                                     //找到 CN 中之index 
//                                     if(CN_set[VN_set[i][mp]][p]==i){
                                    
//                                         VN[i][j] += CN[VN_set[i][mp]][p];		
//                                     }
//                                 }
//                             }		
//                         }					
// 					}
// 				}
				
// 				//total LLR
// 				//decode
// 				for(int i=0;i<n;i++){ 				//go through all VN	
// 					VN_total[i] =  Fn[i]; 			
// 					for(int j=0;j<col[i];j++){				
// 						for(int m=0;m<row[VN_set[i][j]];m++){
// 							if(CN_set[VN_set[i][j]][m]==i){		//找到與VN相連之 CN  idex 
// 								VN_total[i] += CN[VN_set[i][j]][m];
								
// 							}
// 						}
// 					}
					
// 					if(VN_total[i]>=0)
// 						u_hat[i] = 0;
// 					else
// 						u_hat[i] = 1;
// 				}

// 				bool iszerovector = true;
// 				for(int j=0;j<m;j++){
// 					s[j] =0;
// 					for(int i=0;i<row[j];i++){
// 						if(CN_set[j][i]>=0&&u_hat[CN_set[j][i]]==1)
// 							s[j] += 1;
// 					}
// 					s[j] %=2;
// 					if(s[j]!=0){
// 						iszerovector = false;
// 						break;
// 					}	
// 				}
// 				if(iszerovector){
// 					break;
// 				}
// 			}
// 			avgIter[q] +=iter;
			
// 			for(int i=0;i<n;i++){	
// 				if(u[i]!=u_hat[i])
// 					error++;
// 			}
// 			for(int i=0;i<n;i++){
// 				if(u[i]!=u_hat[i]){
// 					frameerror++;
// 					break;
// 				}
// 			}
// 			printf("error: %lld, num: %d, BER: %E, FER: %E Average iteration: %f\r",error,num,((double)error)/((double)(n*num)),((double)frameerror)/(num),(double)avgIter[q]/num);
// 		}
// 		BER[q] = ((double)error)/(double)n/(double)numtime;
// 		FER[q] = ((double)frameerror)/((double)numtime);
// 		printf("BER: %E, FER: %E Average iteration: %f\n",BER[q],FER[q],avgIter[q]/numtime);
// 	} 
// }

//Run until the number of frame errors==frameerror
void LDPC_SPA_until_frame(double* BER,double* FER,int n,int m,int dv,int dc,double R,int** CN_set,int** VN_set,int* row,int* col,double* avgIter,double* SNR_dB,int SNR_L
,int iteration,int frameerror){
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
    //receive signal
	double *y = (double*)malloc(sizeof(double) * n);
	//symdrome
	int* s = (int *)malloc(sizeof(int) * m);
	//transmitted sigal;
	//All 0 codeword
	int *u = (int*)malloc(sizeof(int) * n);	//binary sequence
	int *x = (int*)malloc(sizeof(int) * n);
	for(int i=0;i<n;i++){
		u[i] = 0;
		x[i] = 1;		//Eb=1	Eavg = R*1
	} 
	for(int q=0;q<SNR_L;q++){
		avgIter[q] = 0;
		long long error=0;
		int temp_frameerror=0;
		long long numtime = 0;

		while(temp_frameerror<frameerror){
			
			double sigma = sqrt((1/(2*R))*pow(10,-(SNR_dB[q]/10)));
            Add_AWGN(n,y,x,sigma);
			
			//initialization
			for(int i=0;i<n;i++){
				Fn[i] = 2*y[i]/(pow(sigma,2));	
				VN_total[i] =  Fn[i]; 	
				for(int j=0;j<col[i];j++)
					VN[i][j] = 2*y[i]/(pow(sigma,2));	
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
				for(int j=0;j<m;j++){						//go through all CN	j
					for(int i=0;i<row[j];i++){				//VN_{CN_set[j][i]} connected by CN j
						for(int np=0 ; np<row[j] ; np++){			//n'
							if(CN_set[j][np]>=0&&CN_set[j][i]!=CN_set[j][np]){	//n != n'
								
								for(int f=0;f<col[CN_set[j][np]];f++){
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
							if(CN_set[VN_set[i][j]][m]==i){		
								//找到與VN相連之 CN  idex 
								VN_total[i] += CN[VN_set[i][j]][m];
							}
						}
					}					
					if(VN_total[i]>=0)
						u_hat[i] = 0;
					else
						u_hat[i] = 1;
				}

				//VN update
				for(int i=0;i<n;i++){ 					//go through all VN	i
					for(int j=0;j<col[i];j++){			//CN_{VN_set[j][i]} connected by VN i
						for(int p=0;p<row[VN_set[i][j]];p++){
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
					temp_frameerror++;
					break;
				}
			}
			printf("num: %lld, error: %lld, framerror: %d, BER: %E, FER: %E, Average iteration: %f\r",numtime,error,temp_frameerror,((double)error/(double)n/(double)numtime),(double)temp_frameerror/(double)numtime,avgIter[q]/(double)numtime);
				
			numtime++;
		}
		BER[q] = ((double)error)/(double)n/(double)numtime;
		FER[q] = ((double)frameerror)/((double)numtime);
		printf("BER: %E, FER: %E Average iteration: %f\n",BER[q],FER[q],avgIter[q]/numtime);
	} 
}
// void LDPC_SPA_until_frame(double* BER,double* FER,int n,int m,int dv,int dc,double R,int** CN_set,int** VN_set,int* row,int* col,double* avgIter,double* SNR_dB,int SNR_L
// ,int iteration,int frameerror){
//     //channel information
// 	double* Fn = (double *)malloc(sizeof(double) * n);
// 	//tau
// 	double** tau= (double **)malloc(sizeof(double*) * m);
// 	for(int i=0;i<m;i++){
// 		tau[i] = (double *)malloc(sizeof(double) * dc);
// 	}
// 	//CN update //CN[j][i]  CN j to CN i
// 	double** CN = (double **)malloc(sizeof(double*) * m);
// 	for(int i=0;i<m;i++){
// 		CN[i] = (double *)malloc(sizeof(double) * dc);
// 	}
// 	//VN update //VN[i][j]  VN i to CN j
// 	double** VN = (double **)malloc(sizeof(double*) * n);
// 	for(int i=0;i<n;i++){
// 		VN[i] = (double *)malloc(sizeof(double) * dv);
// 	}
// 	//total LLR VN[i]
// 	double* VN_total = (double *)malloc(sizeof(double) * n);
// 	//decoded sequence
// 	int* u_hat = (int *)malloc(sizeof(int) * n);
//     //receive signal
// 	double *y = (double*)malloc(sizeof(double) * n);
// 	//symdrome
// 	int* s = (int *)malloc(sizeof(int) * m);
// 	//transmitted sigal;
// 	//All 0 codeword
// 	int *u = (int*)malloc(sizeof(int) * n);	//binary sequence
// 	int *x = (int*)malloc(sizeof(int) * n);
// 	for(int i=0;i<n;i++){
// 		u[i] = 0;
// 		x[i] = 1;		//Eb=1	Eavg = R*1
// 	} 
// 	for(int q=0;q<SNR_L;q++){
// 		avgIter[q] = 0;
// 		int error=0;
// 		int temp_frameerror=0;
// 		long long numtime = 0;

// 		while(temp_frameerror<frameerror){
			
// 			double sigma = sqrt((1/(2*R))*pow(10,-(SNR_dB[q]/10)));
//             Add_AWGN(n,y,x,sigma);
			
// 			//initialization
// 			for(int i=0;i<n;i++){
// 				Fn[i] = 2*y[i]/(pow(sigma,2));	
// 				for(int j=0;j<col[i];j++){	
					
// 					VN[i][j] = 2*y[i]/(pow(sigma,2));	
// 				}	
// 			}
			
// 			//iterative decoding
// 			int iter=0;
// 			for(iter=0;iter<iteration;iter++){
// 				//CN update	
// 				for(int j=0;j<m;j++){	
// 					for(int i=0;i<row[j];i++){
// 						tau[j][i]=1;	
// 					}
// 				}

// 				for(int j=0;j<m;j++){						//go through all CN	j
// 					for(int i=0;i<row[j];i++){				//VN i connected by CN j

// 							for(int np=0 ; np<row[j] ; np++){			//n'
// 								if(CN_set[j][np]>=0&&CN_set[j][i]!=CN_set[j][np]){	//n != n'
									
// 									//N(m) set
// 									for(int f=0;f<col[CN_set[j][np]];f++){
// 										//找到VN 中之index 
// 										if(VN_set[CN_set[j][np]][f]==j){
											
// 											tau[j][i] *=tanh(VN[CN_set[j][np]][f]/2);	
// 										}
// 									}
// 								}		
// 							}
// 						//tau is computed.
// 						//CN[j][i] = log((1.0+tau[j][i])/(1.0-tau[j][i]));
// 						if(tau[j][i]==1)
// 							CN[j][i] = DBL_MAX;
// 						else if(tau[j][i]==-1)
// 							CN[j][i] = -DBL_MAX;
// 						else
// 							CN[j][i] = 2*atanh(tau[j][i]);
// 					}
// 				}
// 				//VN update
// 				for(int i=0;i<n;i++){ 				    //go through all VN	i
// 					for(int j=0;j<col[i];j++){			//CN j connected by VN i
						
//                         VN[i][j] = Fn[i];	
//                         for(int mp=0 ; mp<col[i] ; mp++){			//m'
//                             if(VN_set[i][mp]>=0&&VN_set[i][j]!=VN_set[i][mp]){	//m != m'
                                
//                                 for(int p=0;p<row[VN_set[i][mp]];p++){
//                                     //找到 CN 中之index 
//                                     if(CN_set[VN_set[i][mp]][p]==i){
                                    
//                                         VN[i][j] += CN[VN_set[i][mp]][p];		
//                                     }
                                    
//                                 }
                                        
//                             }		
//                         }					
// 					}
// 				}
				
// 				//total LLR
// 				//decode
// 				for(int i=0;i<n;i++){ 				//go through all VN	
// 					VN_total[i] =  Fn[i]; 			
// 					for(int j=0;j<col[i];j++){				
// 						for(int m=0;m<row[VN_set[i][j]];m++){
// 							if(CN_set[VN_set[i][j]][m]==i){		//找到與VN相連之 CN  idex 
// 								VN_total[i] += CN[VN_set[i][j]][m];
								
// 							}
// 						}
	
// 					}					
// 					if(VN_total[i]>=0)
// 						u_hat[i] = 0;
// 					else
// 						u_hat[i] = 1;
// 				}
				
// 				bool iszerovector = true;
// 				for(int j=0;j<m;j++){
// 					s[j] =0;
// 					for(int i=0;i<row[j];i++){
// 						if(CN_set[j][i]>=0&&u_hat[CN_set[j][i]]==1)
// 							s[j] += 1;
// 					}
// 					s[j] %=2;
// 					if(s[j]!=0){
// 						iszerovector = false;
// 						break;
// 					}	
// 				}
				
// 				if(iszerovector){
// 					break;
// 				}
// 			}
// 			avgIter[q] +=iter;

// 			for(int i=0;i<n;i++){	
// 				if(u[i]!=u_hat[i])
// 					error++;
// 			}
// 			for(int i=0;i<n;i++){
// 				if(u[i]!=u_hat[i]){
// 					temp_frameerror++;
// 					break;
// 				}
// 			}
// 			printf("num: %lld, error: %d, framerror: %d, BER: %E, FER: %E, Average iteration: %f\r",numtime,error,temp_frameerror,((double)error/(double)n/(double)numtime),(double)temp_frameerror/(double)numtime,avgIter[q]/(double)numtime);
				
// 			numtime++;
// 		}
// 		BER[q] = ((double)error)/(double)n/(double)numtime;
// 		FER[q] = ((double)frameerror)/((double)numtime);
// 		printf("BER: %E, FER: %E Average iteration: %f\n",BER[q],FER[q],avgIter[q]/numtime);
// 	} 
// }


void LDPC_MSA(double* BER,double* FER,int n,int m,int dv,int dc,double R,int** CN_set,int** VN_set,int* row,int* col,double* avgIter,double* SNR_dB,int SNR_L
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
    //receive signal
	double *y = (double*)malloc(sizeof(double) * n);
	//symdrome
	int* s = (int *)malloc(sizeof(int) * m);
	//transmitted sigal;
	//All 0 codeword
	int *u = (int*)malloc(sizeof(int) * n);	//binary sequence
	int *x = (int*)malloc(sizeof(int) * n);
	for(int i=0;i<n;i++){
		u[i] = 0;
		x[i] = 1;		//Eb=1	Eavg = R*1
	} 
	
	for(int q=0;q<SNR_L;q++){
		avgIter[q] = 0;
		long long error=0;
		int frameerror=0;

		for(int num=0;num<numtime;num++){
			double sigma = sqrt((1/(2*R))*pow(10,-(SNR_dB[q]/10)));
            Add_AWGN(n,y,x,sigma);
			
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
				for(int j=0;j<(m);j++){	
					for(int i=0;i<row[j];i++){
						CN[j][i]=1;	
					}
				}
				for(int j=0;j<m;j++){				//go through all CN	j
					for(int i=0;i<row[j];i++){		//VN i connected by CN j
						
                        double min_beta = DBL_MAX;
                        for(int np=0 ; np<row[j] ; np++){			//n'
                            if(CN_set[j][np]>=0&&CN_set[j][i]!=CN_set[j][np]){	//n != n'
                                
                                for(int f=0;f<col[CN_set[j][np]];f++){
                                    //找到VN 中之index 
                                    if(VN_set[CN_set[j][np]][f]==j){
                                        double temp_min = fabs(VN[CN_set[j][np]][f]);
                                        if(temp_min<min_beta)
                                            min_beta = temp_min;    //find min L VNi to CNj
                                        if(VN[CN_set[j][np]][f]<0)
                                            CN[j][i]*=-1;
                                    }
                                }
                            }		
                        }
						//計算完tau 
						CN[j][i] *= min_beta;
					}
				}

				//VN update
				for(int i=0;i<n;i++){ 				//go through all VN	i
					for(int j=0;j<col[i];j++){		//CN j connected by VN i
						
                        VN[i][j] = Fn[i];	
                        for(int mp=0 ; mp<col[i] ; mp++){			//m'
                            if(VN_set[i][mp]>=0&&VN_set[i][j]!=VN_set[i][mp]){	//m != m'
                                
                                for(int p=0;p<row[VN_set[i][mp]];p++){
                                    //找到 CN 中之index 
                                    if(CN_set[VN_set[i][mp]][p]==i){
                                    
                                        VN[i][j] += CN[VN_set[i][mp]][p];		
                                    }
                                }
                            }		
                        }					
					}
				}
				//total LLR
				//decode
				for(int i=0;i<n;i++){ 				//go through all VN	
					VN_total[i] =  Fn[i]; 			
					for(int j=0;j<col[i];j++){				
						for(int m=0;m<row[VN_set[i][j]];m++){
							if(CN_set[VN_set[i][j]][m]==i){		//找到與VN相連之 CN  idex 
								VN_total[i] += CN[VN_set[i][j]][m];
								
							}
						}
					}					
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
			printf("error: %lld, num: %d, BER: %E, FER: %E Average iteration: %f\r",error,num,((double)error)/((double)(n*num)),((double)frameerror)/(num),(double)avgIter[q]/num);
				
		}
		
		BER[q] = ((double)error)/(double)n/(double)numtime;
		FER[q] = ((double)frameerror)/((double)numtime);
		printf("BER: %E, FER: %E Average iteration: %f\n",BER[q],FER[q],avgIter[q]/numtime);
	} 
}

void LDPC_LBP(double* BER,double* FER,int n,int m,int dv,int dc,double R,int** CN_set,int** VN_set,int* row,int* col,double* avgIter,double* SNR_dB,int SNR_L
,int iteration, int numtime){
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
	//receive signal
	double *y = (double*)malloc(sizeof(double) * n);
	//transmitted sigal;
	//All 0 codeword
	int *u = (int*)malloc(sizeof(int) * n);	//binary sequence
	int *x = (int*)malloc(sizeof(int) * n);
	for(int i=0;i<n;i++){
		u[i] = 0;
		x[i] = 1;		//Eb=1	Eavg = R*1
	} 

	for(int q=0;q<SNR_L;q++){
		avgIter[q] = 0;
		long long error=0;
		int frameerror=0;

		for(int num=0;num<numtime;num++){
			
			double sigma = sqrt((1/(2*R))*pow(10,-(SNR_dB[q]/10)));
			Add_AWGN(n,y,x,sigma);

			//initialization
			for(int i=0;i<n;i++){
				VN_total[i] = 2*y[i]/(pow(sigma,2));
			}
			
			for(int j=0;j<m;j++){	
				for(int i=0;i<row[j];i++){	
					CN[j][i] = 0;	
				}
			}
			//iterative decoding
			int iter=0;
			for(iter=0;iter<iteration;iter++){
				
				//CN update	
				for(int j=0;j<m;j++){	
					for(int i=0;i<row[j];i++){
						tau[j][i] = 1;	
					}
				}
				for(int j=0;j<m;j++){					//go through all CN	j
				
					//對此CN 相連之VN 先update 用上一個 layer 所得之 VN_total 
					for(int i=0;i<row[j];i++){				//VN i connected by CN j
					
						for(int f=0;f<col[CN_set[j][i]];f++){
				
							if(VN_set[CN_set[j][i]][f]==j){
								//找到VN 中之index 
								VN[CN_set[j][i]][f] = VN_total[CN_set[j][i]] - CN[j][i];			
							}
						}
					}
					
					for(int i=0;i<row[j];i++){				//VN i connected by CN j

							int VN_node,VN_idex;
							for(int np=0 ; np<row[j] ; np++){			//n'
								if(CN_set[j][i]==CN_set[j][np]){		//n == n'
								
									VN_node = CN_set[j][np];
									
									for(int f=0;f<col[CN_set[j][np]];f++){
										//找到VN 中之index 										
										if(VN_set[VN_node][f]==j){
											//找到VN 中之index 
											VN_idex = f;
											break;
										}
									}	
								}
								//else if(CN_set[j][i]!=CN_set[j][np])
								else{	//n != n'
									
									for(int f=0;f<col[CN_set[j][np]];f++){
										//找到VN 中之index 
										if(VN_set[CN_set[j][np]][f]==j){
											
											tau[j][i] *=tanh(VN[CN_set[j][np]][f]/2);	
										}
									}	
								}	
							}
							//tau is computed.							
							if(tau[j][i]==1)
								CN[j][i] = DBL_MAX;
							else if(tau[j][i]==-1)
								CN[j][i] = -DBL_MAX;
							else
								CN[j][i] = 2*atanh(tau[j][i]);
							
							VN_total[VN_node] = VN[VN_node][VN_idex] + CN[j][i];
					}
				}
				//total LLR
				//decode
				for(int i=0;i<n;i++){ 				//go through all VN	
					if(VN_total[i]>=0)
						u_hat[i] = 0;
					else
						u_hat[i] = 1;
				}
				//判斷symdrome
				int* s = (int *)malloc(sizeof(int) * m);
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
			printf("error: %lld, num: %d, BER: %E, FER: %E Average iteration: %f\r",error,num,((double)error)/((double)(n*num)),((double)frameerror)/(num),(double)avgIter[q]/num);				
		}
		BER[q] = ((double)error)/(double)n/(double)numtime;
		FER[q] = ((double)frameerror)/((double)numtime);
		printf("BER: %E, FER: %E Average iteration: %f\n",BER[q],FER[q],avgIter[q]/numtime);
	} 
}

void LDPC_Layered_NMSA(double* BER,double* FER,int n,int m,int dv,int dc,double R,int** CN_set,int** VN_set,int* row,int* col,double* avgIter,double* SNR_dB,int SNR_L
,int iteration, int numtime,double normalfactor){
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
	//receive signal
	double *y = (double*)malloc(sizeof(double) * n);
	//noise
	double* noise = (double *)malloc(sizeof(double) * n);
	//transmitted sigal;
	//All 0 codeword
	int *u = (int*)malloc(sizeof(int) * n);	//binary sequence
	int *x = (int*)malloc(sizeof(int) * n);
	for(int i=0;i<n;i++){
		u[i] = 0;
		x[i] = 1;		//Eb=1	Eavg = R*1
	} 

	for(int q=0;q<SNR_L;q++){
		avgIter[q] = 0;
		long long error=0;
		int frameerror=0;
		
		for(int num=0;num<numtime;num++){

			double sigma = sqrt((1/(2*R))*pow(10,-(SNR_dB[q]/10)));
			Add_AWGN(n,y,x,sigma);
			
			//initialization
			for(int i=0;i<n;i++){
				VN_total[i] = 2*y[i]/(pow(sigma,2));	
				for(int j=0;j<col[i];j++){
						
					VN[i][j] =  2*y[i]/(pow(sigma,2));
				}
			}
			for(int j=0;j<m;j++){	
				for(int i=0;i<row[j];i++){
						
					CN[j][i] = 0;	
				}
			}
			//iterative decoding
			int iter=0;
			for(iter=0;iter<iteration;iter++){
				
				//CN update	
				for(int j=0;j<m;j++){	
					for(int i=0;i<row[j];i++){
						
						tau[j][i] = 1;	
					}
				}
				for(int j=0;j<m;j++){					//go through all CN	j

					//對此CN 相連之VN 先update 用上一個 layer 所得之 VN_total 
					for(int i=0;i<row[j];i++){				//VN i connected by CN j
					
						for(int f=0;f<col[CN_set[j][i]];f++){
							
							if(VN_set[CN_set[j][i]][f]==j){
								//找到VN 中之index 
								VN[CN_set[j][i]][f] = VN_total[CN_set[j][i]] - CN[j][i];
							}
						}
					}
				
					for(int i=0;i<row[j];i++){		//相連之 VN 
						
						double min_beta = DBL_MAX;
						int sign = 1;				//代表VN_update 正負號之相乘 
						
						int VN_node,VN_idex;
						for(int np=0 ; np<row[j] ; np++){			//n'
							if(CN_set[j][i]==CN_set[j][np]){		//n == n'
							
								VN_node = CN_set[j][np];
								
								for(int f=0;f<col[CN_set[j][np]];f++){
									//找到VN 中之index 									
									if(VN_set[VN_node][f]==j){
										//找到VN 中之index 
										VN_idex = f;
										break;
									}
								}	
							}
							else{									//n != n'
								
								for(int f=0;f<col[CN_set[j][np]];f++){
									//找到VN 中之index 
									if(VN_set[CN_set[j][np]][f]==j){
										
										//printf("VN idex: %d\n",f);
										if(fabs(VN[CN_set[j][np]][f])<min_beta)
											min_beta = fabs(VN[CN_set[j][np]][f] );
										if(VN[CN_set[j][np]][f]<0)
											sign*=-1;
										
									}
								}	
							}
						}
						//計算完tau 
						CN[j][i] =  sign*min_beta*normalfactor;
						
						VN_total[VN_node] = VN[VN_node][VN_idex] + CN[j][i];
					}
				}
		
				//total LLR
				//decode
				for(int i=0;i<n;i++){ 				//go through all VN	
					
					if(VN_total[i]>=0)
						u_hat[i] = 0;
					else
						u_hat[i] = 1;
				}

				//判斷symdrome
				int* s = (int *)malloc(sizeof(int) * m);
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
			printf("error: %lld, num: %d, BER: %E, FER: %E Average iteration: %f\r",error,num,((double)error)/((double)(n*num)),((double)frameerror)/(num),(double)avgIter[q]/num);
		}
		
		BER[q] = ((double)error)/(double)n/(double)numtime;
		FER[q] = ((double)frameerror)/((double)numtime);
		printf("BER: %E, FER: %E Average iteration: %f\n",BER[q],FER[q],avgIter[q]/numtime);
	} 
}
