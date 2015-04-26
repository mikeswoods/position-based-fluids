/*******************************************************************************
 * PrefixSum.h
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#ifndef PREFIX_SUM_H
#define PREFIX_SUM_H

#include <iostream>
#include <cmath>
#include <vector>
#include "MSAOpenCL.h"

/******************************************************************************/

class PrefixSum
{
    private:
        // OpenCL manager
        msa::OpenCL& openCL;

        std::vector<msa::OpenCLBuffer*> ScanPartialSums;
        unsigned int ElementsAllocated;
        unsigned int LevelsAllocated;
    
        bool IsPowerOfTwo(int n)
        {
            return ((n&(n-1))==0) ;
        }
        
        int floorPow2(int n)
        {
            int exp;
            frexp((float)n, &exp);
            return 1 << (exp - 1);
        }
    
        void loadKernels();
    
    protected:
        int GROUP_SIZE;
    
        int CreatePartialSumBuffers(unsigned int count);

        void ReleasePartialSums();

        int PreScan(size_t *global
                   ,size_t *local
                   ,size_t shared
                   ,msa::OpenCLBuffer& output_data
                   ,msa::OpenCLBuffer& input_data
                   ,unsigned int n
                   ,int group_index
                   ,int base_index);

        int PreScanStoreSum(size_t *global
                           ,size_t *local
                           ,size_t shared
                           ,msa::OpenCLBuffer& output_data
                           ,msa::OpenCLBuffer& input_data
                           ,msa::OpenCLBuffer& partial_sums
                           ,unsigned int n
                           ,int group_index
                           ,int base_index);

        int PreScanStoreSumNonPowerOfTwo(size_t *global
                                        ,size_t *local
                                        ,size_t shared
                                        ,msa::OpenCLBuffer& output_data
                                        ,msa::OpenCLBuffer& input_data
                                        ,msa::OpenCLBuffer& partial_sums
                                        ,unsigned int n
                                        ,int group_index
                                        ,int base_index);
    
        int PreScanNonPowerOfTwo(size_t *global
                                ,size_t *local
                                ,size_t shared
                                ,msa::OpenCLBuffer& output_data
                                ,msa::OpenCLBuffer& input_data
                                ,unsigned int n
                                ,int group_index
                                ,int base_index);
    
        int UniformAdd(size_t *global
                      ,size_t *local
                      ,msa::OpenCLBuffer& output_data
                      ,msa::OpenCLBuffer& partial_sums
                      ,unsigned int n
                      ,unsigned int group_offset
                      ,unsigned int base_index);
        
        int PreScanBufferRecursive(msa::OpenCLBuffer& output_data
                                  ,msa::OpenCLBuffer& input_data
                                  ,int max_group_size
                                  ,int max_work_item_count
                                  ,int element_count
                                  ,int level);
        void PreScanBuffer(msa::OpenCLBuffer& output_data
                          ,msa::OpenCLBuffer& input_data
                          ,unsigned int max_group_size
                          ,unsigned int max_work_item_count
                          ,unsigned int element_count);
    
        public:
            PrefixSum(msa::OpenCL& openCL, int GROUP_SIZE = 256);

            void scan(msa::OpenCLBuffer& output_data
                     ,msa::OpenCLBuffer& input_data
                     ,unsigned int element_count);
};

/******************************************************************************/

#endif
