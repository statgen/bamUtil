/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <iostream>
#include <stdio.h>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <unistd.h>

#include "MemoryMap.h"

#if defined(WIN32)
SYSTEM_INFO MemoryMap::system_info;
static bool need_system_info = 1;
#else
#include <sys/mman.h>
#endif

#ifndef MAP_POPULATE
#define MAP_POPULATE 0x0000
#endif
#ifndef MAP_NONBLOCK
#define MAP_NONBLOCK 0x0000
#endif

MemoryMap::MemoryMap()
{
    constructor_clear();
#if defined(WIN32)
    if (need_system_info)
    {
        GetSystemInfo(&system_info);
        need_system_info = false;
    }
    page_size = system_info.dwAllocationGranularity;
#else
    page_size = sysconf(_SC_PAGE_SIZE);
#endif
};

MemoryMap::~MemoryMap()
{
    destructor_clear();
};

void MemoryMap::debug_print()
{
#if defined(WIN32)
    std::cout << "fd = " << file_handle << std::endl;
#else
    std::cout << "fd = " << fd << std::endl;
#endif
    std::cout << "data = 0x" << std::hex << data << std::endl;
    std::cout << "offset = 0x" << std::hex << offset << std::endl;
    std::cout << "mapped_length = 0x" << std::hex << mapped_length << std::endl;
    std::cout << "total_length = 0x" << std::hex << total_length << std::endl;
    std::cout << "page_size = 0x" << std::hex << page_size << std::endl;
};

void MemoryMap::constructor_clear()
{
#if defined(WIN32)
    file_handle = 0;
#else
    fd = -1;
#endif
    data = (void *) NULL;
    offset = 0;
    mapped_length = 0;
    total_length = 0;
    useMemoryMapFlag = true;
};

void MemoryMap::destructor_clear()
{
    if (data!=NULL)
    {
#if defined(WIN32)
        UnmapViewOfFile((LPVOID) data);
#else
        munmap(data, mapped_length);
#endif
    }
#if !defined(WIN32)
    if (fd!=-1)
    {
        ::close(fd);
    }
#endif
    constructor_clear();
};

#if defined(WIN32)
// XXX
// TODO implement win32 ::open and ::create
#else
bool MemoryMap::open(const char * file, int flags)
{

    struct stat buf;

    int mmap_prot_flag = PROT_READ;
    if (flags != O_RDONLY) mmap_prot_flag = PROT_WRITE;

    fd = ::open(file, flags);
    if (fd==-1)
    {
        fprintf(stderr, "MemoryMap::open: file %s not found\n", (const char *) file);
        constructor_clear();
        return true;
    }
    if (fstat(fd, &buf))
    {
        perror("MemoryMap::open");
        constructor_clear();
        return true;
    }
    mapped_length = buf.st_size;
    total_length = mapped_length;

    if (useMemoryMapFlag)
    {

        int additionalFlags = 0;

        // try this for amusement ... not portable:
//        additionalFlags |= MAP_HUGETLB;
// #define USE_LOCKED_MMAP
#if defined(USE_LOCKED_MMAP)
        // MAP_POPULATE only makes sense if we are reading
        // the data
        if (flags != O_RDONLY)
        {
            // furthermore, according to Linux mmap page, populate only
            // works if the map is private.
            additionalFlags |= MAP_POPULATE;
            additionalFlags |= MAP_PRIVATE;
        }
        else
        {
            additionalFlags |= MAP_SHARED;
        }
#else
additionalFlags |= MAP_SHARED;
#endif

        data = ::mmap(
                   NULL,           // start
                   mapped_length,
                   mmap_prot_flag, // protection flags
                   additionalFlags,
                   fd,
                   offset
               );
        if (data == MAP_FAILED)
        {
            ::close(fd);
            std::cerr << "Error: Attempting to map " << mapped_length << " bytes of file "
                      << file << ":" << std::endl;
            perror("MemoryMap::open");
            constructor_clear();
            return true;
        }

#if defined(USE_LOCKED_MMAP)
        //
        // non-POSIX, so non portable.
        // This call is limited by the RLIMIT_MEMLOCK resource.
        //
        // In bash, "ulimit -l" shows the limit.
        //
        if (mlock(data, mapped_length))
        {
            std::cerr << "Warning: Attempting to lock " << mapped_length << " bytes of file " << file << ":" << std::endl;
            perror("unable to lock memory");
            // not a fatal error, so continue
        }
#endif

        // these really don't appear to have any greatly useful effect on
        // Linux.  Last time I checked, the effect on Solaris and AIX was
        // exactly what was documented and was significant.
        //
        madvise(data, mapped_length, MADV_WILLNEED);   // warning, this could hose the system

    }
    else
    {
        data = (void *) malloc(mapped_length);
        if (data==NULL)
        {
            ::close(fd);
            perror("MemoryMap::open");
            constructor_clear();
            return true;
        }
        ssize_t resultSize = read(fd, data, mapped_length);
        if (resultSize!=(ssize_t) mapped_length)
        {
            ::close(fd);
            perror("MemoryMap::open");
            constructor_clear();
            return true;
        }
    }
    return false;
}

bool MemoryMap::create(const char *file, size_t size)
{

    if (file==NULL)
    {
        data = calloc(size, 1);
        if (data==NULL) return true;
    }
    else
    {
        int mmap_prot_flag = PROT_READ | PROT_WRITE;

        fd = ::open(file, O_RDWR|O_CREAT|O_TRUNC, 0666);
        if (fd==-1)
        {
            fprintf(stderr, "MemoryMap::open: can't create file '%s'\n",(const char *) file);
            constructor_clear();
            return true;
        }

        lseek(fd, (off_t) size - 1, SEEK_SET);
        char ch = 0;
        if(write(fd, &ch, 1)!=1) {
            perror("MemoryMap::create:");
            throw std::logic_error("unable to write at end of file");
        }

        data = ::mmap(
                   NULL,           // start
                   size,
                   mmap_prot_flag, // protection flags
                   MAP_SHARED,     // share/execute/etc flags
                   fd,
                   offset
               );
        if (data == MAP_FAILED)
        {
            ::close(fd);
            unlink(file);
            perror("MemoryMap::open");
            constructor_clear();
            return true;
        }
        mapped_length = size;
        total_length = size;
    }
    return false;
}
#endif

bool MemoryMap::create(size_t size)
{
    return create(NULL, size);
}

bool MemoryMap::close()
{
    destructor_clear();
    return false;
};

void MemoryMap::test()
{
    int result;

    result = this->open("test/test_memmap_data.txt");
    assert(result == 0);
    assert(data!=NULL);
    assert(mapped_length == 183);   // length of the above file
    close();

    // now try non memory mapped (direct slow file I/O)
    useMemoryMap(false);
    result = this->open("test/test_memmap_data.txt");
    assert(result == 0);
    assert(data!=NULL);
    assert(mapped_length == 183);   // length of the above file
    close();
}

int MemoryMap::prefetch()
{
    int sum = 0;
    size_t i;

    for (i=0; i<mapped_length; i += page_size) sum += *(i + (char *) data);

    return sum;
}

#if defined(TEST)
//
// compile test using:
// g++ -DTEST -o testMemoryMap MemoryMap.cpp Error.o -lz
//

int main(int argc, const char *argv)
{
    MemoryMap map;

    map.test();

//  map.debug_print();

    exit(0);
}
#endif

