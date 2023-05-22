#ifndef channel_randSeq_h
#define channel_randSeq_h
#ifdef  channel_randSeq_h

#include "random_num_gen_unit.h"
#include <vector>

using namespace std;

class channel_randSeq {
    public:
        channel_randSeq(unsigned int max_cwdLgIn);
        ~channel_randSeq();
        vector<random_num_gen_unit> randNoise;
        unsigned char      *receivedHD;
        unsigned char      *receivedLevel;
        unsigned char      *receivedEncoded;
        vector<unsigned char>   recvHD,recvLv,recvSym;
        unsigned int        LLRTableNum;
        double              RBER;
        void set_randSeed(unsigned int seed);
        void set_mode(unsigned int randSeqChanModeIn,unsigned int noiseModeIn);
        void set_mode_RBER(unsigned int randSeqChanModeIn,double RBERIn);
        void set_mode_biasedRBER(unsigned int randSeqChanModeIn,double RBER,double biasedRBER);
        void add_AWGN(unsigned char *cwd,unsigned int transfer_length);
        void add_fixerr(unsigned char *cwd,unsigned int err_num,unsigned int transfer_length);
        void add_N4fixerr(unsigned char *cwd,unsigned int err_num,double SCR,double SER,unsigned int transfer_length);
        void add_HRE(unsigned char *cwd,unsigned int err_num,unsigned int transfer_length);
        void setFullSeedNoise(int *seedReg);
        void getPRNG(string& seed);
        void setPRNG(string& seed);
        vector<double>          vtModeC3;
        vector<double>          ISweight;
    protected:
        unsigned int            transferLg;
        vector<double>          vtModeN2RBER,vtModeN4RBER,vtModeN6RBER,vtModeN8RBER;
        vector<double>          vtModeFN4RBER;
        //unsigned int            randSeqChanMode;
        //unsigned int            noiseLevel;
        unsigned long long      randNumOut[16];
        unsigned char           LLREncode(unsigned char LevelIn);
        vector<unsigned long long>  vtModeC3LongInt;
        void                    CDFgen(vector<double> &cdf,double RBER);
        void                    biasedCDFgen(vector<double> &cdf,double RBER,double biasedRBER);
};


inline unsigned char channel_randSeq::LLREncode(unsigned char LevelIn)
{
    unsigned char Lv2EncMap2[2]={3,7};
    unsigned char Lv2EncMap4[4]={3,0,4,7};
    unsigned char Lv2EncMap6[6]={3,1,0,4,5,7};
    unsigned char Lv2EncMap8[8]={3,2,1,0,4,5,6,7};
    if(LevelIn>=LLRTableNum)    cout<<"(channel_randSeq::LLREncode) Error: LevelIn out of range."<<endl;
    if(LLRTableNum==2) {
        return Lv2EncMap2[LevelIn];
    } else if(LLRTableNum==4) {
        return Lv2EncMap4[LevelIn];
    } else if(LLRTableNum==6) {
        return Lv2EncMap6[LevelIn];
    } else if(LLRTableNum==8) {
        return Lv2EncMap8[LevelIn];
    } else {
        cout<<"error: LLRTableNum = "<<LLRTableNum<<" is invalid"<<endl;
        getchar();
    }
    return 0;
}

#endif
#endif
