Date: 2023-05-17
# Rand API document
## class
### random_num_gen_unit
It is a pseudo-random number generator.
## member functions
### void random_num_gen_unit::setSeed(unsigned int seedIn)
Set seed for the generator.
### void random_num_gen_unit::setSeed(string &seed)
Set the internal state of the generator.
### void random_num_gen_unit::getSeed(string &seed) const
Get the internal state of the generator.
### unsigned long long random_num_gen_unit::operator()()
Generate a random number. It outputs a 32-bits unsigned integer.
## class
### channel_randSeq
It is a channel model. It can generate hard bit and soft bit data.
## member functions
### channel_randSeq::channel_randSeq(unsigned int max_cwdLgIn)
The parameter *max_cwdLgIn* is the length of the codeword.
### void channel_randSeq::set_mode_RBER(unsigned int randSeqChanModeIn, double RBERIn)
Use hard bit data if *randSeqChanModeIn*=1. Use soft bit data if *randSeqChanModeIn*=7. And set the RBER of the channel to the value of *RBERIn*.
### void channel_randSeq::add_AWGN(unsigned char *cwd, unsigned int transfer_length)
Add AWGN to the codeword *cwd* of length *transfer_length*. See the data member section for the response to the function.
## data member
### unsigned char *channel_randSeq::receivedEncoded
After calling `channel_randSeq::add_AWGN`, it points to the response to the `channel_randSeq::add_AWGN`. Its length equals to the *transfer_length* which was given by the user. Each of its elements is in the range [0, 7]. The first bit is the sign bit. The last two lsbs are represented as symbols(levels).
