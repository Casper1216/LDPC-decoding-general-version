#include <iostream>
#include <iomanip>
using namespace std;
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "random_num_gen_unit.h"
#include "channel_randSeq.h"
unsigned int randSeedNoise[16] = {
            0x484d8043,
            0x884e9f85,
            0xc8e3843b,
            0x31838d2a,
            0x9abe6a26,
            0x66c65516,
            0x8932cfc8,
            0xff8c8081,
            0x67c66973,
            0x51ff4aec,
            0x29cdbaab,
            0xf2fbe346,
            0x7cc254f8,
            0x1be8e78d,
            0x765a2e63,
            0x339fc99a};

channel_randSeq::channel_randSeq(unsigned int max_cwdLgIn)
{
    int i;
    randNoise.resize(16);
    transferLg=max_cwdLgIn;
    for(i=0;i<16;i++) randNoise[i].setSeed(randSeedNoise[i]);
    recvHD.resize(transferLg);
    recvLv.resize(transferLg);
    recvSym.resize(transferLg);
    receivedHD = &(recvHD[0]);
    receivedLevel = &(recvLv[0]);
    receivedEncoded  = &(recvSym[0]);
    LLRTableNum = 0;
    RBER = 0;
    //randSeqChanMode = 0;
    //noiseLevel =0;

    double defalut_N2RBER[15]={0.0100,0.0090,0.0080,0.0070,0.0060,0.0055,0.0050,0.0045,0.0040,0.0035,0.0030,0.0025,0.0020,0.0015,0.0010};
    double defalut_FN4RBER[15]={0.0150,0.0140,0.0130,0.0120,0.0110,0.0100,0.0090,0.0080,0.0070,0.0060,0.0050,0.0040,0.0030,0.0020,0.0010};
    double defalut_N4RBER[19]={0.0200,0.0180,0.0170,0.0160,0.0150,0.0140,0.0135,0.0130,0.0125,0.0120,0.0115,0.0110,0.0100,0.0090,0.0080,0.0070,0.0060,0.0050,0.0040};
    double defalut_N6RBER[19]={0.0220,0.0200,0.0180,0.0170,0.0160,0.0155,0.0150,0.0145,0.0140,0.0135,0.0130,0.0125,0.0120,0.0110,0.0100,0.0090,0.0080,0.0070,0.0060};
    double defalut_N8RBER[19]={0.0220,0.0200,0.0180,0.0170,0.0160,0.0155,0.0150,0.0145,0.0140,0.0135,0.0130,0.0125,0.0120,0.0110,0.0100,0.0090,0.0080,0.0070,0.0060};
    vtModeN2RBER=vector<double>(defalut_N2RBER,defalut_N2RBER+15);
    vtModeFN4RBER=vector<double>(defalut_FN4RBER,defalut_FN4RBER+15);
    vtModeN4RBER=vector<double>(defalut_N4RBER,defalut_N4RBER+19);
    vtModeN6RBER=vector<double>(defalut_N6RBER,defalut_N6RBER+19);
    vtModeN8RBER=vector<double>(defalut_N8RBER,defalut_N8RBER+19);
    vtModeC3.resize(9);
    vtModeC3LongInt.resize(9);
}
channel_randSeq::~channel_randSeq() {}

void channel_randSeq::set_randSeed(unsigned int seed)
{
    /* random random seed : 32-bit integer */
    int i;
    srand(seed);
    for(i=0;i<16;i++){
        randSeedNoise[i]=((rand()&0xff)<<24) ^ ((rand()&0xff)<<16) ^ ((rand()&0xff)<<8) ^ (rand()&0xff) ;
    }
    for(i=0;i<16;i++) randNoise[i].setSeed(randSeedNoise[i]);
}

void channel_randSeq::set_mode(unsigned int randSeqChanModeIn,unsigned int noiseLevelIn)
{
    //noiseLevel = noiseLevelIn;
    switch(randSeqChanModeIn) {
        case 0:
        case 1:
            RBER = vtModeN2RBER[noiseLevelIn];
            break;
        case 2:
            RBER = vtModeFN4RBER[noiseLevelIn];
            break;
        case 3:
            RBER = vtModeN4RBER[noiseLevelIn];
            break;
        case 5:
            RBER = vtModeN6RBER[noiseLevelIn];
            break;
        case 7:
            RBER = vtModeN8RBER[noiseLevelIn];
            break;
        default:
            cout<<"randSeqChanMode = "<<randSeqChanModeIn<<" is invalid"<<endl; getchar(); exit(-1);
    }
    set_mode_RBER(randSeqChanModeIn,RBER);
}

void channel_randSeq::set_mode_RBER(unsigned int randSeqChanModeIn,double RBERIn) {
    unsigned int i;
    double LLRTableFlt[8];
    //randSeqChanMode = randSeqChanModeIn;
    LLRTableNum = (randSeqChanModeIn%2==0)? randSeqChanModeIn+2:randSeqChanModeIn+1;
    RBER = RBERIn;
    CDFgen(vtModeC3,RBER);
    for(i=0;i<LLRTableNum+1;i++) {
        cout<<vtModeC3[i]<<"\t";
        vtModeC3LongInt[i] = (unsigned long long) (vtModeC3[i]*pow(2,32));
        cout<<vtModeC3LongInt[i]<<endl;
    }
    cout<<"RBER = "<<RBER<<endl;
    for(i=0;i<LLRTableNum;i++){
        LLRTableFlt[i] = log( (vtModeC3[i+1]-vtModeC3[i])/(vtModeC3[LLRTableNum-i]-vtModeC3[LLRTableNum-1-i]));
        cout<<"LLRTableNum = "<<LLRTableNum<<", LLRTableFlt["<<i<<"] = "<<LLRTableFlt[i]<<endl;
    }
    //memset(receivedLevel,63,sizeof(unsigned char)*max_cwdLg);
}

void channel_randSeq::set_mode_biasedRBER(unsigned int randSeqChanModeIn,double RBERin,double biasedRBER)
{
    unsigned int i;
    vector<double> cdf_orig;
    //randSeqChanMode = randSeqChanModeIn;
    LLRTableNum = (randSeqChanModeIn==0)? 2:randSeqChanModeIn+1;
    RBER = biasedRBER;
    cdf_orig.resize(9);
    ISweight.resize(8);
    CDFgen(cdf_orig,RBERin);
    biasedCDFgen(vtModeC3,RBERin,biasedRBER);
    for(i=0;i<LLRTableNum+1;i++) {
        vtModeC3LongInt[i] = (unsigned long long) (vtModeC3[i]*pow(2,32));
    }
    for(i=0;i<LLRTableNum;i++){
        ISweight[i] = (cdf_orig[i+1]-cdf_orig[i])/(vtModeC3[i+1]-vtModeC3[i]);
        cout<<"LLRTableNum = "<<LLRTableNum<<", ISweight["<<i<<"] = "<<ISweight[i]<<endl;
    }
    //memset(receivedLevel,63,sizeof(unsigned char)*max_cwdLg);
}

void channel_randSeq::add_AWGN(unsigned char *cwd,unsigned int transfer_length) {
    unsigned int j=0,k=0;
    unsigned int addNoiseCnt = 0;
    if(transfer_length>transferLg){
        transferLg=transfer_length;
        recvHD.resize(transferLg);
        recvLv.resize(transferLg);
        recvSym.resize(transferLg);
        receivedHD = &(recvHD[0]);
        receivedLevel = &(recvLv[0]);
        receivedEncoded  = &(recvSym[0]);
    }
    while (addNoiseCnt<transfer_length) {
        j=15-(addNoiseCnt%16);
        randNumOut[j] = randNoise[j].genRand();
        for(k=0;k<LLRTableNum;k++){
            if(randNumOut[j]<vtModeC3LongInt[k+1]){
                recvLv.at(addNoiseCnt) = (cwd[addNoiseCnt]==0)? k:LLRTableNum-1-k;
                break;
            }
        }
        if(k==LLRTableNum){
            cout<<"rand"<<hex<<randNumOut[j]<<dec<<"add_AWGN k=="<<LLRTableNum<<endl;
            recvLv.at(addNoiseCnt) = (cwd[addNoiseCnt]==0)? LLRTableNum-1:0;
        }
        recvSym[addNoiseCnt] = LLREncode(recvLv[addNoiseCnt]);
        recvHD[addNoiseCnt] = (recvSym[addNoiseCnt])>>2;
        addNoiseCnt++;
    }
    for(k=0;k<j;k++)   randNumOut[j] = randNoise[j].genRand();
}

void channel_randSeq::add_fixerr(unsigned char *cwd,unsigned int err_num,unsigned int transfer_length) {
    unsigned int i=15,addNoiseCnt = 0;
    static unsigned char Lv2EncMap2[2]={3,7};
    if(transfer_length>transferLg){
        transferLg=transfer_length;
        recvHD.resize(transferLg);
        recvLv.resize(transferLg);
        recvSym.resize(transferLg);
        receivedHD = &(recvHD[0]);
        receivedLevel = &(recvLv[0]);
        receivedEncoded  = &(recvSym[0]);
    }
    recvHD.assign(transfer_length,0);
    while (addNoiseCnt<err_num){
        randNumOut[i] = randNoise[i].genRand() % transfer_length;
        if(recvHD.at(randNumOut[i])==0){
            recvHD.at(randNumOut[i])=1;
            addNoiseCnt++;
        }
        i=15-((16-i)%16);
    }
    for(i=0;i<transfer_length;i++){
        recvLv[i] = cwd[i] ^ recvHD[i];
        recvSym[i]= Lv2EncMap2[ recvLv[i] ];
        recvHD[i] = recvLv[i];
    }
}

void channel_randSeq::add_N4fixerr(unsigned char *cwd,unsigned int err_num,double SCR,double SER,unsigned int transfer_length)
{
    unsigned int i=15,addNoiseCnt = 0;
    unsigned int SEbit,WEbit,WCbit;//,SCbit;
    static unsigned char Lv2EncMap4[4]={3,0,4,7};
    if(transfer_length>transferLg){
        transferLg=transfer_length;
        recvHD.resize(transferLg);
        recvLv.resize(transferLg);
        recvSym.resize(transferLg);
        receivedHD = &(recvHD[0]);
        receivedLevel = &(recvLv[0]);
        receivedEncoded  = &(recvSym[0]);
    }
    WCbit=(unsigned int)floor((1.0-SCR)*(transfer_length-err_num));
    //SCbit=transfer_length-err_num-WCbit;
    //SCbit=(unsigned int)(SCR*(transfer_length-err_num));
    //WCbit=transfer_length-err_num-SCbit;
    SEbit=(unsigned int)(SER*err_num);
    WEbit=err_num-SEbit;
    recvHD.assign(transfer_length,0);   // default strong correct
    while (addNoiseCnt<SEbit) {         // strong error location
        randNumOut[0] = randNoise[0].genRand() % transfer_length;
        if(recvHD[randNumOut[0]]==0){
            recvHD[randNumOut[0]]=3;
            addNoiseCnt++;
        }
    }
    addNoiseCnt=0;
    while (addNoiseCnt<WEbit) {         // weak error location
        randNumOut[0] = randNoise[0].genRand() % transfer_length;
        if(recvHD[randNumOut[0]]==0){
            recvHD[randNumOut[0]]=2;
            addNoiseCnt++;
        }
    }
    addNoiseCnt=0;
    while (addNoiseCnt<WCbit) {         // weak correct location
        randNumOut[0] = randNoise[0].genRand() % transfer_length;
        if(recvHD[randNumOut[0]]==0){
            recvHD[randNumOut[0]]=1;
            addNoiseCnt++;
        }
    }

    for(i=0;i<transfer_length;i++){
        recvLv[i] = (cwd[i]==0)? recvHD[i]:3-recvHD[i];
        recvSym[i]= Lv2EncMap4[ recvLv[i] ];
        recvHD[i] = (recvSym[i])>>2;
    }
}


void channel_randSeq::add_HRE(unsigned char *cwd,unsigned int err_num,unsigned int transfer_length)
{
    unsigned int j,k,addNoiseCnt = 0;
    unsigned char val;
    if(transfer_length>transferLg){
        transferLg=transfer_length;
        recvHD.resize(transferLg);
        recvLv.resize(transferLg);
        recvSym.resize(transferLg);
        receivedHD = &(recvHD[0]);
        receivedLevel = &(recvLv[0]);
        receivedEncoded  = &(recvSym[0]);
    }
    recvHD.assign(transferLg,0);
    while (addNoiseCnt<err_num) {
        randNumOut[0] = randNoise[0].genRand() % transfer_length;
        if(recvHD[randNumOut[0]]==0){
            recvHD[randNumOut[0]]=1;
            addNoiseCnt++;
        }
    }
    addNoiseCnt = 0;
    while (addNoiseCnt<transfer_length) {
        j=15-(addNoiseCnt%16);
        randNumOut[j] = randNoise[j].genRand();
        val = cwd[addNoiseCnt] ^ recvHD[addNoiseCnt];
        for(k=0;k<LLRTableNum;k++){
            if(randNumOut[j]<vtModeC3LongInt[k+1]){
                recvLv[addNoiseCnt] = (val==0)? k:LLRTableNum-1-k;
                break;
            }
        }
        recvSym[addNoiseCnt] = LLREncode(recvLv[addNoiseCnt]);
        recvHD[addNoiseCnt] = (recvSym[addNoiseCnt])>>2;
        addNoiseCnt++;
    }
}

void channel_randSeq::setFullSeedNoise(int *seedReg){
    int i,j;

    for(i = 0; i < 16; i ++) {
        for(j = 0; j < 43; j++)
            randNoise[i].LFSRReg[j] = seedReg[1317-i*43+j];
        for(j = 0; j < 37; j++)
            randNoise[i].CASRReg[j] = seedReg[592-i*37+j];
    }
}

void channel_randSeq::getPRNG(string& seed)
{
    unsigned int i;
    string randUnitSeed;
    seed.clear();
    for(i=0;i<16;i++){
        randNoise[i].getSeed(randUnitSeed);
        seed=seed+randUnitSeed;
    }
}
void channel_randSeq::setPRNG(string& seed)
{
    unsigned int i;
    string randUnitSeed;
    for(i=0;i<16;i++){
        randUnitSeed=seed.substr(i*20,20);
        randNoise[i].setSeed(randUnitSeed);
    }
}

void channel_randSeq::CDFgen(vector<double> &cdf,double RBER)
{
    double sigma,var_high,var_low,var_test;
    double RBER_c=1;

    var_high=4/pow(-log2(1-pow(0.00000001,1/1.1064))/0.3073,1/0.8935);
    var_low=4/pow(-log2(1-pow(0.99999999,1/1.1064))/0.3073,1/0.8935);
    if(RBER==0 && LLRTableNum==2){
        cdf[0]=0;
        cdf[1]=1;
        cdf[2]=1;
        return;
    }

    while(fabs(RBER_c-RBER)>1e-15){
        var_test=(var_high+var_low)/2;
        RBER_c=erfc(1/sqrt(2*var_test))/2;

        if(RBER_c>RBER){
            var_high=var_test;
        }else{
            var_low=var_test;
        }
    }

    var_test=var_low;
    sigma=sqrt(var_test);
    if(LLRTableNum==2){
        cdf[0]=0;
        cdf[1]=1-erfc(1/sqrt(2*var_test))/2;
        cdf[2]=1;
    }else if(LLRTableNum==4){
        cdf[0]=0;
        cdf[1]=1-erfc((1-1*sigma/2)/sqrt(2*var_test))/2;
        cdf[2]=1-erfc(1/sqrt(2*var_test))/2;
        cdf[3]=1-erfc((1+1*sigma/2)/sqrt(2*var_test))/2;
        cdf[4]=1;
    }else if(LLRTableNum==6){
        cdf[0]=0;
        cdf[1]=1-erfc((1-2*sigma/2)/sqrt(2*var_test))/2;
        cdf[2]=1-erfc((1-1*sigma/2)/sqrt(2*var_test))/2;
        cdf[3]=1-erfc(1/sqrt(2*var_test))/2;
        cdf[4]=1-erfc((1+1*sigma/2)/sqrt(2*var_test))/2;
        cdf[5]=1-erfc((1+2*sigma/2)/sqrt(2*var_test))/2;
        cdf[6]=1;
    }else if(LLRTableNum==8){
        cdf[0]=0;
        cdf[1]=1-erfc((1-3*sigma/2)/sqrt(2*var_test))/2;
        cdf[2]=1-erfc((1-2*sigma/2)/sqrt(2*var_test))/2;
        cdf[3]=1-erfc((1-1*sigma/2)/sqrt(2*var_test))/2;
        cdf[4]=1-erfc(1/sqrt(2*var_test))/2;
        cdf[5]=1-erfc((1+1*sigma/2)/sqrt(2*var_test))/2;
        cdf[6]=1-erfc((1+2*sigma/2)/sqrt(2*var_test))/2;
        cdf[7]=1-erfc((1+3*sigma/2)/sqrt(2*var_test))/2;
        cdf[8]=1;
    }
    return;
}

void channel_randSeq::biasedCDFgen(vector<double> &cdf,double RBER,double biasedRBER)
{
    double sigma,var_high,var_low,var_test;
    double mean_high,mean_low,mean;
    double RBER_c=1;

    var_high=4/pow(-log2(1-pow(0.00000001,1/1.1064))/0.3073,1/0.8935);
    var_low=4/pow(-log2(1-pow(0.99999999,1/1.1064))/0.3073,1/0.8935);
    while(fabs(RBER_c-RBER)>1e-15){
        var_test=(var_high+var_low)/2;
        RBER_c=erfc(1/sqrt(2*var_test))/2;
        if(RBER_c>RBER){
            var_high=var_test;
        }else{
            var_low=var_test;
        }
    }
    var_test=var_low;
    sigma=sqrt(var_test);

    if(biasedRBER==0.5) mean=0;
    else{
        RBER_c=1;
        mean_high=1.60;
        mean_low=-1.60;
        while(fabs(RBER_c-biasedRBER)>1e-10){
            mean=(mean_high+mean_low)/2;
            RBER_c=erfc(-mean/sqrt(2*var_test))/2;
            if(RBER_c>biasedRBER){
                mean_high=mean;
            }else{
                mean_low=mean;
            }
        }
        mean=mean_low;
    }
    //cout<<"mean = "<<mean<<endl;

    if(LLRTableNum==2){
        cdf[0]=0;
        cdf[1]=1-erfc(-mean/sqrt(2*var_test))/2;
        cdf[2]=1;
    }
    if(LLRTableNum==4){
        cdf[0]=0;
        cdf[1]=1-erfc((-mean-1*sigma/2)/sqrt(2*var_test))/2;
        cdf[2]=1-erfc(-mean/sqrt(2*var_test))/2;
        cdf[3]=1-erfc((-mean+1*sigma/2)/sqrt(2*var_test))/2;
        cdf[4]=1;
    }
    if(LLRTableNum==6){
        cdf[0]=0;
        cdf[1]=1-erfc((-mean-2*sigma/2)/sqrt(2*var_test))/2;
        cdf[2]=1-erfc((-mean-1*sigma/2)/sqrt(2*var_test))/2;
        cdf[3]=1-erfc(-mean/sqrt(2*var_test))/2;
        cdf[4]=1-erfc((-mean+1*sigma/2)/sqrt(2*var_test))/2;
        cdf[5]=1-erfc((-mean+2*sigma/2)/sqrt(2*var_test))/2;
        cdf[6]=1;
    }
    if(LLRTableNum==8){
        cdf[0]=0;
        cdf[1]=1-erfc((-mean-3*sigma/2)/sqrt(2*var_test))/2;
        cdf[2]=1-erfc((-mean-2*sigma/2)/sqrt(2*var_test))/2;
        cdf[3]=1-erfc((-mean-1*sigma/2)/sqrt(2*var_test))/2;
        cdf[4]=1-erfc(-mean/sqrt(2*var_test))/2;
        cdf[5]=1-erfc((-mean+1*sigma/2)/sqrt(2*var_test))/2;
        cdf[6]=1-erfc((-mean+2*sigma/2)/sqrt(2*var_test))/2;
        cdf[7]=1-erfc((-mean+3*sigma/2)/sqrt(2*var_test))/2;
        cdf[8]=1;
    }
}
