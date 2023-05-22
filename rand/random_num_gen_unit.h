#ifndef random_num_gen_unit_h
#define random_num_gen_unit_h

#include <bitset>
#include <cstdlib>
#include <cassert>
#include <string>

using std::bitset;
using std::string;

static uint32_t prngGenRand() {
	FILE *f = fopen("/dev/urandom", "r");
	if(f == NULL) {
		perror("There is no urandom in the OS");
		return 0;
	}
	uint32_t ret[2] = {0};
	if(fgets((char*)ret, sizeof(uint32_t) + 1, f) == NULL) {
		perror("Fail to read random number from urandom.");
		fclose(f);
		return 0;
	}
	fclose(f);
	return *ret;
}
static char *gen_uuid(char *buf) {
	FILE *f = fopen("/proc/sys/kernel/random/uuid", "r");
	if(f == NULL) {
		perror("There is no random/uuid in the OS");
		buf[0] = '\0';
		return NULL;
	}
	if(fgets(buf, 37, f) == NULL) {
		perror("Fail to read uuid from random/uuid.");
		fclose(f);
		buf[0] = '\0';
		return NULL;
	}
	fclose(f);
	return buf;
}

class random_num_gen_unit {
    public:
        random_num_gen_unit(){
            LFSRReg.reset();
            CASRReg.reset();
            CASRVarReg.reset();
        };
        ~random_num_gen_unit(){};
        void setSeed(unsigned int seedIn){
            LFSRReg = bitset<43>(seedIn);
            CASRReg = bitset<37>(seedIn);
        };
        void setFullSeed(bitset<43> & LFSR, bitset<37> & CASR){
            LFSRReg=LFSR;
            CASRReg=CASR;
        };
        void getSeed(string &seed) const;
        void setSeed(string &seed);
        unsigned long long genRand();
        unsigned long long operator()() { return genRand(); }
        bitset<43> LFSRReg;
        bitset<37> CASRReg,CASRVarReg;
    protected:
};

inline void random_num_gen_unit::getSeed(string &seed) const
{
    unsigned int i,num;
    bitset<80> buff;
    seed.clear();
    for(i=0;i<43;i++)buff[i]=LFSRReg[i];
    for(i=0;i<37;i++)buff[43+i]=CASRReg[i];
    for(i=0;i<80;i+=4){
        num=buff[i]*8+buff[i+1]*4+buff[i+2]*2+buff[i+3];
        if(num>9)   seed.push_back(num-10+'a');
        else        seed.push_back(num+'0');
    }
}
inline void random_num_gen_unit::setSeed(string &seed)
{
    assert(seed.length() >= 20);
    unsigned int i,num;
    bitset<80> buff;
    for(i=0;i<20;i++){
        num=(seed[i]>='a')?seed[i]-'a'+10:seed[i]-'0';
        buff[4*i]=(num>>3) & 1;
        buff[4*i+1]=(num>>2) & 1;
        buff[4*i+2]=(num>>1) & 1;
        buff[4*i+3]=num & 1;
    }
    for(i=0;i<43;i++)LFSRReg[i]=buff[i];
    for(i=0;i<37;i++)CASRReg[i]=buff[43+i];
}

inline unsigned long long random_num_gen_unit::genRand() {
    int i;
    bool LFSRBuf;
    bool CASRBuf;
    unsigned long long int ret = 0;
//    cout<<"LFSRReg init "<<LFSRReg<<endl;
//    cout<<"CASRReg init "<<CASRReg<<endl;
    for(i=31; i>=0; i--){
        ret = (ret<<1)^(LFSRReg[i]^CASRReg[i]);
    }
//    cout<<"ret = "<<ret.to_ulong()<<endl;
//    getchar();

    LFSRBuf = LFSRReg[42];
    LFSRReg=LFSRReg<<1;
    LFSRReg[0] = LFSRBuf;
    LFSRReg[1] = LFSRReg[1]^LFSRBuf;
    LFSRReg[20] = LFSRReg[20]^LFSRBuf;
    LFSRReg[41] = LFSRReg[41]^LFSRBuf;

    CASRVarReg=CASRReg;
    CASRBuf = CASRReg[36];
    CASRReg=CASRReg<<1;
    CASRReg[0] = CASRBuf;
    CASRReg[27] = CASRReg[27]^CASRVarReg[27];

    CASRBuf = CASRVarReg[0];
    CASRVarReg=CASRVarReg>>1;
    CASRVarReg[36] = CASRBuf;
    CASRReg ^= CASRVarReg;
    return ret;
}
#endif
