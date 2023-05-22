#ifndef random_num_gen_msg_h
#define random_num_gen_msg_h
#ifdef  random_num_gen_msg_h

#include <vector>
class random_num_gen_msg: public random_num_gen_unit {
    public:
        random_num_gen_msg();
        ~random_num_gen_msg();
        void genMsg(unsigned char* msg,unsigned int msgLg);
        void genZeroMsg(unsigned char* msg,unsigned int msgLg);
};



#endif
#endif
