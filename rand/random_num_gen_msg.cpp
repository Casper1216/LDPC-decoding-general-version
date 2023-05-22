#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "random_num_gen_unit.h"
#include "random_num_gen_msg.h"
random_num_gen_msg::random_num_gen_msg()
{
    setSeed(0x484d8043);    // default seed
    //nextIn.assign(20,0);
}
random_num_gen_msg::~random_num_gen_msg(){}

void random_num_gen_msg::genMsg(unsigned char* msg,unsigned int msgLg)
{
    unsigned long long ret;
    unsigned int genLg = 0,i;
    while (genLg<=msgLg-16) {
        ret = genRand();
        for(i=0; i<16; i++) {
            msg[genLg+i] = ( ret>>(i+4) ) & 1;
        }
        genLg+=16;
    }
    if(genLg<msgLg) {
        ret = genRand();
        for(i=0; i<msgLg-genLg; i++) {
            msg[genLg+i] = ( ret>>i ) & 1;
        }
        genLg = msgLg;
    }
}

void random_num_gen_msg::genZeroMsg(unsigned char* msg,unsigned int msgLg)
{
    memset(msg,0,sizeof(unsigned char)*msgLg);
}
