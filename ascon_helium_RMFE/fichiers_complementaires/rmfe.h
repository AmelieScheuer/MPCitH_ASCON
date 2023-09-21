/******************************************************************************
* Copyright (c) 2023 Thales
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* *
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
* *
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
-
******************************************************************************/



#include <stdint.h>
#include "field2.h"


/*-------------------RMFE (4,9)-------------------*/

uint16_t phi_matrix_4_9[] = {
0x1e2, 0x1e1, 0xb3, 0xb1
};

uint8_t psi_matrix_4_9[] = {
0xf, 0x3, 0x3, 0x7, 0x4, 0x3, 0x7, 0x9, 0x6
};


field::GF2E phi_4_9(std::vector<uint8_t> input){
    uint16_t output=0;
    uint8_t tmp;
    for(size_t i=0;i<input.size();i++){
        tmp=input[i];
        for(size_t j=0;j<8;j++){
            if((tmp&1) != 0){
                output^=phi_matrix_4_9[3-j-(i<<3)];
            }
            tmp>>=1;
        }
    }
    return field::GF2E(output);
}

std::vector<uint8_t> psi_4_9(field::GF2E input){
    std::vector<uint8_t> vec=input.to_bytes();
    uint8_t tmp=0;
    uint8_t output=0;
    for(size_t k=0;k<vec.size();k++){
        tmp=vec[k];
        for(size_t i=0;i<8;i++){
            if((tmp&1) != 0){
                output^=psi_matrix_4_9[(k<<3)+i];
            }
            tmp>>=1;
        }
    }
    std::vector<uint8_t> res(1);
    res[0]=output;
    return res;
}



/*-------------------RMFE (8,21)-------------------*/

uint32_t phi_matrix_8_21[] = {
0x135398, 0x1f1732, 0x1d9d10, 0x15e5dc,
0x14d9f, 0x160e05, 0x1e76cb, 0xd0936
};

uint8_t psi_matrix_8_21[] = {
0xff, 0x36, 0x36, 0x39, 0x33, 0x3d, 0x38,
0x4d, 0x94, 0xfa, 0x54, 0xfa, 0xc8, 0xe1,
0x4a, 0xf7, 0xc0, 0x9c, 0x72, 0xd3, 0x77
};


field::GF2E phi_8_21(std::vector<uint8_t> input){
    uint32_t output=0;
    uint8_t tmp;
    for(size_t i=0;i<input.size();i++){
        tmp=input[i];
        for(size_t j=0;j<8;j++){
            if((tmp&1) != 0){
                output^=phi_matrix_8_21[7-j-(i<<3)];
            }
            tmp>>=1;
        }
    }
    return field::GF2E(output);
}

std::vector<uint8_t> psi_8_21(field::GF2E input){
    std::vector<uint8_t> vec=input.to_bytes();
    uint8_t tmp=0;
    uint8_t output=0;
    for(size_t k=0;k<vec.size();k++){
        tmp=vec[k];
        for(size_t i=0;i<8;i++){
            if((tmp&1) != 0){
                output^=psi_matrix_8_21[(k<<3)+i];
            }
            tmp>>=1;
        }
    }
    std::vector<uint8_t> res(1);
    res[0]=output;
    return res;
}

/*-------------------RMFE (8,24)-------------------*/

uint32_t phi_matrix_8_24[] = {
0xe1fd33, 0x37c792, 0x41f445, 0xda277b,
0xc709fa, 0xe17ff5, 0x7aacc9, 0x113358
};

uint8_t psi_matrix_8_24[] = {
0xff, 0x36, 0x36, 0x39, 0x33, 0x3d, 0x38, 0x3f,
0x61, 0x45, 0xb0, 0xd0, 0xa8, 0x78, 0x67, 0xc6,
0xe5, 0xf6, 0x53, 0x11, 0x22, 0x2d, 0xf5, 0xa
};


field::GF2E phi_8_24(std::vector<uint8_t> input){
    uint32_t output=0;
    uint8_t tmp;
    for(size_t i=0;i<input.size();i++){
        tmp=input[i];
        for(size_t j=0;j<8;j++){
            if((tmp&1) != 0){
                output^=phi_matrix_8_24[7-j-(i<<3)];
            }
            tmp>>=1;
        }
    }
    return field::GF2E(output);
}

std::vector<uint8_t> psi_8_24(field::GF2E input){
    std::vector<uint8_t> vec=input.to_bytes();
    uint8_t tmp=0;
    uint8_t output=0;
    for(size_t k=0;k<vec.size();k++){
        tmp=vec[k];
        for(size_t i=0;i<8;i++){
            if((tmp&1) != 0){
                output^=psi_matrix_8_24[(k<<3)+i];
            }
            tmp>>=1;
        }
    }
    std::vector<uint8_t> res(1);
    res[0]=output;
    return res;
}

/*-------------------RMFE (8,32)-------------------*/

uint32_t phi_matrix_8_32[] = {
0xcce7b98e, 0x768a302, 0x77455083, 0xfa4f438b,
0x68e4c61e, 0x3db2e96b, 0xb0b8fa61, 0xa36bdc91
};

uint8_t psi_matrix_8_32[] = {
0xff, 0x36, 0x36, 0x36, 0x39, 0x30, 0x30, 0x3f,
0x97, 0x12, 0x1c, 0x12, 0x19, 0x1d, 0x40, 0x52,
0xbf, 0x62, 0x7a, 0x36, 0x68, 0x2d, 0x2c, 0x29,
0x2f, 0x23, 0x87, 0x0e, 0xaa, 0x73, 0x34, 0x68
};


field::GF2E phi_8_32(std::vector<uint8_t> input){
    uint32_t output=0;
    uint8_t tmp;
    for(size_t i=0;i<input.size();i++){
        tmp=input[i];
        for(size_t j=0;j<8;j++){
            if((tmp&1) != 0){
                output^=phi_matrix_8_32[7-j-(i<<3)];
            }
            tmp>>=1;
        }
    }
    return field::GF2E(output);
}

std::vector<uint8_t> psi_8_32(field::GF2E input){
    std::vector<uint8_t> vec=input.to_bytes();
    uint8_t tmp=0;
    uint8_t output=0;
    for(size_t k=0;k<vec.size();k++){
        tmp=vec[k];
        for(size_t i=0;i<8;i++){
            if((tmp&1) != 0){
                output^=psi_matrix_8_32[(k<<3)+i];
            }
            tmp>>=1;
        }
    }
    std::vector<uint8_t> res(1);
    res[0]=output;
    return res;
}

/*-------------------RMFE (16,45)-------------------*/

uint64_t phi_matrix_16_45[] = {
0x1f3b22539ead, 0x1f3b22539e2c, 0x15d23c62eb36, 0x15d23c62ebc8,
0x1f8c9cb1cd35, 0xaf7172f3ee7, 0x1e7b7502d9fc, 0x1146cf83b583,
0x1d910a6477c7, 0x1f3a768c9ef, 0x13ce3a882bb2, 0x1b682a1af2a6,
0x1bb7b11ab81b, 0x12ae9788f577, 0x8a4a59dcf88, 0xf3f0f8e7a0b
};

uint16_t psi_matrix_16_45[] = {
0xffff, 0x3663, 0x3636, 0x398d, 0x3366, 0x3d98, 0x38d9, 0x3fff, 0x3663,
0x3636, 0x398d, 0x3366, 0x3d98, 0x38d9, 0x3fff, 0x5a86, 0x11cd, 0xa87c,
0x8d4f, 0xd3fc, 0x71a5, 0x8865, 0xadd0, 0x3e00, 0x4c30, 0xb81, 0x7456,
0xe342, 0x4ec5, 0x698f, 0x3633, 0x5a46, 0x726c, 0x935b, 0x717a, 0x53ec,
0x885, 0xb4d4, 0x41cc, 0xc16a, 0xb702, 0x2070, 0x5911, 0x6607, 0x88c4
};


field::GF2E phi_16_45(std::vector<uint8_t> input){
    uint64_t output=0;
    uint8_t tmp;
    for(size_t i=0;i<input.size();i++){
        tmp=input[i];
        for(size_t j=0;j<8;j++){
            if((tmp&1) != 0){
                output^=phi_matrix_16_45[15-j-(i<<3)];
            }
            tmp>>=1;
        }
    }
    return field::GF2E(output);
}

std::vector<uint8_t> psi_16_45(field::GF2E input){
    std::vector<uint8_t> vec=input.to_bytes();
    uint8_t tmp=0;
    uint16_t output=0;
    for(size_t k=0;k<vec.size();k++){
        tmp=vec[k];
        for(size_t i=0;i<8;i++){
            if((tmp&1) != 0){
                output^=psi_matrix_16_45[(k<<3)+i];
            }
            tmp>>=1;
        }
    }
    std::vector<uint8_t> res(2);
    res[0]=output & 0xff;
    res[1]=output>>8;
    return res;
}

/*-------------------RMFE (16,48)-------------------*/

uint64_t phi_matrix_16_48[] = {
0xdef1809e3aa4, 0xdef1809e3a25, 0x94a100ebd338, 0x94a100ebd3c6,
0x19b63c96eeef, 0xeedbd69abf33, 0x5adc0cf562c, 0x77b20a8fdbb,
0xf28233bf665d, 0x27b0980087cf, 0x421cf370e839, 0xaaaaae3b8bd1,
0xcfde2047e2c8, 0xed8161f45210, 0xcd7e940f94b5, 0x271e29235c22
};

uint16_t psi_matrix_16_48[] = {
0xffff, 0x3663, 0x3636, 0x398d, 0x3366, 0x3d98, 0x38d9, 0x3fff,
0x3663, 0x3636, 0x398d, 0x3366, 0x3d98, 0x38d9, 0x3fff, 0x3663,
0x4afc, 0x0ee5, 0x057b, 0x070b, 0x5ed1, 0x83b7, 0x982b, 0x85fe,
0x2cd5, 0x006b, 0xcd85, 0x2fda, 0xe466, 0x954c, 0x67ca, 0xb089,
0xa50a, 0x967d, 0x3132, 0xf21a, 0xd733, 0x3e92, 0xfbcf, 0x8c93,
0x50ac, 0xd26b, 0x358f, 0x10d3, 0x38c7, 0x67b5, 0xc365, 0xe411
};

field::GF2E phi_16_48(std::vector<uint8_t> input){
    uint64_t output=0;
    uint8_t tmp;
    for(size_t i=0;i<input.size();i++){
        tmp=input[i];
        for(size_t j=0;j<8;j++){
            if((tmp&1) != 0){
                output^=phi_matrix_16_48[15-j-(i<<3)];
            }
            tmp>>=1;
        }
    }
    return field::GF2E(output);
}

std::vector<uint8_t> psi_16_48(field::GF2E input){
    std::vector<uint8_t> vec=input.to_bytes();
    uint8_t tmp=0;
    uint16_t output=0;
    for(size_t k=0;k<vec.size();k++){
        tmp=vec[k];
        for(size_t i=0;i<8;i++){
            if((tmp&1) != 0){
                output^=psi_matrix_16_48[(k<<3)+i];
            }
            tmp>>=1;
        }
    }
    std::vector<uint8_t> res(2);
    res[0]=output & 0xff;
    res[1]=output>>8;
    return res;
}

/*-------------------RMFE (16,64)-------------------*/

uint64_t phi_matrix_16_64[] = {
0x69f2e570f66763a3, 0x69f2e570f6676328, 0x4ea3465f5bba42c2, 0x4ea3465f5bba4230,
0x7119bc6201946faf, 0xa1eed7bdc3b45e55, 0x8382243130ebc8fe, 0xc2433629a89e2c0a,
0xea35c434c2790daf, 0x8c4e820a852943bc, 0x59a11cd439c3baa8, 0xaeec232d24c621a7,
0xa680da1fa828224b, 0x9de4b0100f45604d, 0xe0cc9989fd20dc00, 0xdba8f3865a4d9e6c
};

uint16_t psi_matrix_16_64[] = {
0xffff, 0x3663, 0x369c, 0x360f, 0x3993, 0x30fc, 0x309f, 0x3f93, 0x396c, 0x306f, 0x3f03, 0x36fc, 0x390f, 0x3f63, 0x39fc, 0x3fff,
0x4c47, 0x5437, 0x24f2, 0x7476, 0xdb6f, 0x8b89, 0xeaf, 0x59b2, 0x8a92, 0xeebf, 0xdc86, 0x7491, 0x88ef, 0xca47, 0xde65, 0x7a2c,
0xdd95, 0x8259, 0x2815, 0x592e, 0x1f70, 0xb323, 0xf74b, 0x1c3c, 0x786a, 0xa5ae, 0xcbeb, 0xf75e, 0xcf2f, 0x9300, 0xa214, 0xe4dd,
0x5129, 0xb6e2, 0xd973, 0x6a57, 0x3420, 0x42b2, 0x9bfb, 0x41f4, 0x5f16, 0xecf6, 0xed35, 0xedd1, 0x477b, 0xcc26, 0x7ddd, 0x9d8d
};


field::GF2E phi_16_64(std::vector<uint8_t> input){
    uint64_t output=0;
    uint8_t tmp;
    for(size_t i=0;i<input.size();i++){
        tmp=input[i];
        for(size_t j=0;j<8;j++){
            if((tmp&1) != 0){
                output^=phi_matrix_16_64[15-j-(i<<3)];
            }
            tmp>>=1;
        }
    }
    return field::GF2E(output);
}

std::vector<uint8_t> psi_16_64(field::GF2E input){
    std::vector<uint8_t> vec=input.to_bytes();
    uint8_t tmp=0;
    uint16_t output=0;
    for(size_t k=0;k<vec.size();k++){
        tmp=vec[k];
        for(size_t i=0;i<8;i++){
            if((tmp&1) != 0){
                output^=psi_matrix_16_64[(k<<3)+i];
            }
            tmp>>=1;
        }
    }
    std::vector<uint8_t> res(2);
    res[0]=output & 0xff;
    res[1]=output>>8;
    return res;
}







class RMFE {
    public :
        static size_t vecsize;
        static size_t fieldsize;
        static std::function<field::GF2E(std::vector<uint8_t>)> phi;
        static std::function<std::vector<uint8_t>(field::GF2E)> psi;
        static void init_RMFE(const banquet_instance_t &instance);
};

std::function<field::GF2E(std::vector<uint8_t>)> RMFE::phi = nullptr;
std::function<std::vector<uint8_t>(field::GF2E)> RMFE::psi = nullptr;
size_t RMFE::fieldsize = 0;
size_t RMFE::vecsize=0;

void RMFE::init_RMFE(const banquet_instance_t &instance){
    fieldsize=instance.RMFE_fieldsize;
    vecsize=instance.RMFE_vecsize;
    switch (vecsize){
        case 4: {
            if(fieldsize==9){
                phi=phi_4_9;
                psi=psi_4_9;
            }
            else{
                throw std::runtime_error("RMFE is not implemented for this specific field size and vector size.");
            }
        } break;
        case 8: {
            switch (fieldsize){
                case 21: {
                    phi=phi_8_21;
                    psi=psi_8_21;
                } break;
                case 24: {
                    phi=phi_8_24;
                    psi=psi_8_24;
                } break;
                case 32: {
                    phi=phi_8_32;
                    psi=psi_8_32;
                } break;
                default:
                    throw std::runtime_error("RMFE is not implemented for this specific field size and vector size.");
            }
        } break;
        case 16: {
            switch (fieldsize){
                case 45: {
                    phi=phi_16_45;
                    psi=psi_16_45;
                } break;
                case 48: {
                    phi=phi_16_48;
                    psi=psi_16_48;
                } break;
                case 64: {
                    phi=phi_16_64;
                    psi=psi_16_64;
                } break;
                default:
                    throw std::runtime_error("RMFE is not implemented for this specific field size and vector size.");
            }
        } break;
        default:
            throw std::runtime_error("RMFE is not implemented for this specific field size and vector size.");
    } 
}
