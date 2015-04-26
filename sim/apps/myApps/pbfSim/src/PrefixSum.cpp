/*******************************************************************************
 * PrefixSum.cpp
 * - Based off of Apple Inc.'s parallel prefix scan OpenCL implementation
 *
 * CIS563: Physically Based Animation final project
 * Class wrapper created by Michael Woods & Michael O'Meara
 * Original code courtesy of Apple Inc. 2008. All Rights Reserved.
 
 ******************************************************************************/

#include "PrefixSum.h"

/******************************************************************************/
//
// Abstract:   This example shows how to perform an efficient parallel prefix sum (aka Scan)
//             using OpenCL.  Scan is a common data parallel primitive which can be used for
//             variety of different operations -- this example uses local memory for storing
//             partial sums and avoids memory bank conflicts on architectures which serialize
//             memory operations that are serviced on the same memory bank by offsetting the
//             loads and stores based on the size of the local group and the number of
//             memory banks (see appropriate macro definition).  As a result, this example
//             requires that the local group size > 1.
//
// Version:    <1.0>
//
// Disclaimer: IMPORTANT:  This Apple software is supplied to you by Apple Inc. ("Apple")
//             in consideration of your agreement to the following terms, and your use,
//             installation, modification or redistribution of this Apple software
//             constitutes acceptance of these terms.  If you do not agree with these
//             terms, please do not use, install, modify or redistribute this Apple
//             software.
//
//             In consideration of your agreement to abide by the following terms, and
//             subject to these terms, Apple grants you a personal, non - exclusive
//             license, under Apple's copyrights in this original Apple software ( the
//             "Apple Software" ), to use, reproduce, modify and redistribute the Apple
//             Software, with or without modifications, in source and / or binary forms;
//             provided that if you redistribute the Apple Software in its entirety and
//             without modifications, you must retain this notice and the following text
//             and disclaimers in all such redistributions of the Apple Software. Neither
//             the name, trademarks, service marks or logos of Apple Inc. may be used to
//             endorse or promote products derived from the Apple Software without specific
//             prior written permission from Apple.  Except as expressly stated in this
//             notice, no other rights or licenses, express or implied, are granted by
//             Apple herein, including but not limited to any patent rights that may be
//             infringed by your derivative works or by other works in which the Apple
//             Software may be incorporated.
//
//             The Apple Software is provided by Apple on an "AS IS" basis.  APPLE MAKES NO
//             WARRANTIES, EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION THE IMPLIED
//             WARRANTIES OF NON - INFRINGEMENT, MERCHANTABILITY AND FITNESS FOR A
//             PARTICULAR PURPOSE, REGARDING THE APPLE SOFTWARE OR ITS USE AND OPERATION
//             ALONE OR IN COMBINATION WITH YOUR PRODUCTS.
//
//             IN NO EVENT SHALL APPLE BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL OR
//             CONSEQUENTIAL DAMAGES ( INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//             SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//             INTERRUPTION ) ARISING IN ANY WAY OUT OF THE USE, REPRODUCTION, MODIFICATION
//             AND / OR DISTRIBUTION OF THE APPLE SOFTWARE, HOWEVER CAUSED AND WHETHER
//             UNDER THEORY OF CONTRACT, TORT ( INCLUDING NEGLIGENCE ), STRICT LIABILITY OR
//             OTHERWISE, EVEN IF APPLE HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Copyright ( C ) 2008 Apple Inc. All Rights Reserved.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#define NUM_BANKS (16)
#define MAX_ERROR (1e-7)

/******************************************************************************/

using namespace std;

/******************************************************************************/

PrefixSum::PrefixSum(msa::OpenCL& _openCL
                    ,int _GROUP_SIZE) :
    openCL(_openCL),
    GROUP_SIZE(_GROUP_SIZE),
    ScanPartialSums(NULL),
    ElementsAllocated(0),
    LevelsAllocated(0)
{
    this->loadKernels();
}

/**
 * Loads and initializes the kernels used in the prefix sum scan
 */
void PrefixSum::loadKernels()
{
    this->openCL.loadProgramFromFile("kernels/Scan.cl");
    
    this->openCL.loadKernel("PreScanKernel");
    this->openCL.loadKernel("PreScanStoreSumKernel");
    this->openCL.loadKernel("PreScanStoreSumNonPowerOfTwoKernel");
    this->openCL.loadKernel("PreScanNonPowerOfTwoKernel");
    this->openCL.loadKernel("UniformAddKernel");
}

/**
 * CreatePartialSumBuffers
 */
int PrefixSum::CreatePartialSumBuffers(unsigned int count)
{
    ElementsAllocated          = count;
    unsigned int group_size    = GROUP_SIZE;
    unsigned int element_count = count;
    
    int level = 0;
    do {
        unsigned int group_count = (int)fmax(1, (int)ceil((float)element_count / (2.0f * group_size)));

        if (group_count > 1) {
            level++;
        }

        element_count = group_count;
        
    } while (element_count > 1);
    
    for (int i = 0; i < level; i++) {
        ScanPartialSums.push_back(new msa::OpenCLBuffer());
    }

    LevelsAllocated = level;
    
    element_count = count;
    level = 0;
    
    do {
        unsigned int group_count = (int)fmax(1, (int)ceil((float)element_count / (2.0f * group_size)));
       
        if (group_count > 1) {
            ScanPartialSums[level++]->initBuffer(group_count * sizeof(float));
        }
        
        element_count = group_count;
        
    } while (element_count > 1);
    
    return CL_SUCCESS;
}

/**
 * ReleasePartialSums
 */
void PrefixSum::ReleasePartialSums()
{
    unsigned int i;

    ElementsAllocated = 0;
    LevelsAllocated = 0;
}

/**
 * PrefixSum
 */
int PrefixSum::PreScan(size_t *global
                      ,size_t *local
                      ,size_t shared
                      ,msa::OpenCLBuffer& output_data
                      ,msa::OpenCLBuffer& input_data
                      ,unsigned int n
                      ,int group_index
                      ,int base_index)
{
    auto kernel = this->openCL.kernel("PreScanKernel");
    kernel->setArg(0, output_data);
    kernel->setArg(1, input_data);
    kernel->setArg(2, (void*)NULL, shared);
    kernel->setArg(3, group_index);
    kernel->setArg(4, base_index);
    kernel->setArg(5, static_cast<int>(n));

    kernel->run1D(*global, *local);

    return CL_SUCCESS;
}

/**
 * PrefixSum
 */
int PrefixSum::PreScanStoreSum(size_t *global
                              ,size_t *local
                              ,size_t shared
                              ,msa::OpenCLBuffer& output_data
                              ,msa::OpenCLBuffer& input_data
                              ,msa::OpenCLBuffer& partial_sums
                              ,unsigned int n
                              ,int group_index
                              ,int base_index)
{
    auto kernel = this->openCL.kernel("PreScanStoreSumKernel");
    kernel->setArg(0, output_data);
    kernel->setArg(1, input_data);
    kernel->setArg(2, partial_sums);
    kernel->setArg(3, (void*)NULL, shared);
    kernel->setArg(4, group_index);
    kernel->setArg(5, base_index);
    kernel->setArg(6, static_cast<int>(n));
    
    kernel->run1D(*global, *local);

    return CL_SUCCESS;
}

/**
 * PreScanStoreSumNonPowerOfTwo
 */
int PrefixSum::PreScanStoreSumNonPowerOfTwo(size_t *global
                                           ,size_t *local
                                           ,size_t shared
                                           ,msa::OpenCLBuffer& output_data
                                           ,msa::OpenCLBuffer& input_data
                                           ,msa::OpenCLBuffer& partial_sums
                                           ,unsigned int n
                                           ,int group_index
                                           ,int base_index)
{
    auto kernel = this->openCL.kernel("PreScanStoreSumNonPowerOfTwoKernel");
    kernel->setArg(0, output_data);
    kernel->setArg(1, input_data);
    kernel->setArg(2, partial_sums);
    kernel->setArg(3, (void*)NULL, shared);
    kernel->setArg(4, group_index);
    kernel->setArg(5, base_index);
    kernel->setArg(6, static_cast<int>(n));
    
    kernel->run1D(*global, *local);

    return CL_SUCCESS;
}

/**
 * PreScanNonPowerOfTwo
 */
int PrefixSum::PreScanNonPowerOfTwo(size_t *global
                                   ,size_t *local
                                   ,size_t shared
                                   ,msa::OpenCLBuffer& output_data
                                   ,msa::OpenCLBuffer& input_data
                                   ,unsigned int n
                                   ,int group_index
                                   ,int base_index)
{
    auto kernel = this->openCL.kernel("PreScanNonPowerOfTwoKernel");
    
    kernel->setArg(0, output_data);
    kernel->setArg(1, input_data);
    kernel->setArg(2, (void*)NULL, shared);
    kernel->setArg(3, group_index);
    kernel->setArg(4, base_index);
    kernel->setArg(5, static_cast<int>(n));
    
    kernel->run1D(*global, *local);

    return CL_SUCCESS;
}

/**
 * UniformAdd
 */
int PrefixSum::UniformAdd(size_t *global
                         ,size_t *local
                         ,msa::OpenCLBuffer& output_data
                         ,msa::OpenCLBuffer& partial_sums
                         ,unsigned int n
                         ,unsigned int group_offset
                         ,unsigned int base_index)
{
    auto kernel = this->openCL.kernel("UniformAddKernel");
    
    kernel->setArg(0, output_data);
    kernel->setArg(1, partial_sums);
    kernel->setArg(2, (void*)NULL, sizeof(float));
    kernel->setArg(3, static_cast<int>(group_offset));
    kernel->setArg(4, static_cast<int>(base_index));
    kernel->setArg(5, static_cast<int>(n));
    
    kernel->run1D(*global, *local);

    return CL_SUCCESS;
}

/**
 * PreScanBufferRecursive
 */
int PrefixSum::PreScanBufferRecursive(msa::OpenCLBuffer& output_data
                                     ,msa::OpenCLBuffer& input_data
                                     ,int max_group_size
                                     ,int max_work_item_count
                                     ,int element_count
                                     ,int level)
{
    unsigned int group_size = max_group_size;
    unsigned int group_count = (int)fmax(1.0f, (int)ceil((float)element_count / (2.0f * group_size)));
    unsigned int work_item_count = 0;
    
    if (group_count > 1) {
        work_item_count = group_size;
    } else if (IsPowerOfTwo(element_count)) {
        work_item_count = element_count / 2;
    } else {
        work_item_count = floorPow2(element_count);
    }

    work_item_count = (work_item_count > max_work_item_count) ? max_work_item_count : work_item_count;
    
    unsigned int element_count_per_group = work_item_count * 2;
    unsigned int last_group_element_count = element_count - (group_count-1) * element_count_per_group;
    unsigned int remaining_work_item_count = (int)fmax(1.0f, last_group_element_count / 2);
    remaining_work_item_count = (remaining_work_item_count > max_work_item_count) ? max_work_item_count : remaining_work_item_count;
    unsigned int remainder = 0;
    size_t last_shared = 0;
    
    
    if (last_group_element_count != element_count_per_group) {
        remainder = 1;
        
        if (!IsPowerOfTwo(last_group_element_count)) {
            remaining_work_item_count = floorPow2(last_group_element_count);
        }

        remaining_work_item_count = (remaining_work_item_count > max_work_item_count) ? max_work_item_count : remaining_work_item_count;
        unsigned int padding = (2 * remaining_work_item_count) / NUM_BANKS;
        last_shared = sizeof(float) * (2 * remaining_work_item_count + padding);
    }
    
    remaining_work_item_count = (remaining_work_item_count > max_work_item_count) ? max_work_item_count : remaining_work_item_count;
    size_t global[] = { (int)fmax(1, group_count - remainder) * work_item_count, 1 };
    size_t local[]  = { work_item_count, 1 };
    
    unsigned int padding = element_count_per_group / NUM_BANKS;
    size_t shared = sizeof(float) * (element_count_per_group + padding);
    
    auto partial_sums = ScanPartialSums[level];
    int err = CL_SUCCESS;
    
    if (group_count > 1) {

        err = PreScanStoreSum(global, local, shared, output_data, input_data, *partial_sums, work_item_count * 2, 0, 0);

        if (err != CL_SUCCESS) {
            return err;
        }
        
        if (remainder) {

            size_t last_global[] = { 1 * remaining_work_item_count, 1 };
            size_t last_local[]  = { remaining_work_item_count, 1 };
            
            err = PreScanStoreSumNonPowerOfTwo(last_global
                                              ,last_local
                                              ,last_shared
                                              ,output_data
                                              ,input_data
                                              ,*partial_sums
                                              ,last_group_element_count
                                              ,group_count - 1
                                              ,element_count - last_group_element_count);
            
            if (err != CL_SUCCESS) {
                return err;
            }
        }
        
        err = PreScanBufferRecursive(*partial_sums, *partial_sums, max_group_size, max_work_item_count, group_count, level + 1);

        if (err != CL_SUCCESS) {
            return err;
        }
        
        err = UniformAdd(global, local, output_data, *partial_sums, element_count - last_group_element_count, 0, 0);

        if (err != CL_SUCCESS) {
            return err;
        }
        
        if (remainder) {

            size_t last_global[] = { 1 * remaining_work_item_count, 1 };
            size_t last_local[]  = { remaining_work_item_count, 1 };
            
            err = UniformAdd(last_global, last_local
                            ,output_data, *partial_sums
                            ,last_group_element_count
                            ,group_count - 1
                            ,element_count - last_group_element_count);
            
            if (err != CL_SUCCESS) {
                return err;
            }
        }

    } else if (IsPowerOfTwo(element_count)) {

        err = PreScan(global, local, shared, output_data, input_data, work_item_count * 2, 0, 0);

        if (err != CL_SUCCESS) {
            return err;
        }

    } else {

        err = PreScanNonPowerOfTwo(global, local, shared, output_data, input_data, element_count, 0, 0);

        if (err != CL_SUCCESS) {
            return err;
        }
    }
    
    return CL_SUCCESS;
}

/**
 * PreScanBuffer
 */
void PrefixSum::PreScanBuffer(msa::OpenCLBuffer& output_data
                             ,msa::OpenCLBuffer& input_data
                             ,unsigned int max_group_size
                             ,unsigned int max_work_item_count
                             ,unsigned int element_count)
{
    PreScanBufferRecursive(output_data
                          ,input_data
                          ,max_group_size
                          ,max_work_item_count
                          ,element_count
                          ,0);
}

void PrefixSum::scan(msa::OpenCLBuffer& output_data
                    ,msa::OpenCLBuffer& input_data
                    ,unsigned int element_count)
{
    CreatePartialSumBuffers(element_count);
    PreScanBufferRecursive(output_data
                          ,input_data
                          ,this->GROUP_SIZE
                          ,this->GROUP_SIZE
                          ,element_count
                          ,0);
}

