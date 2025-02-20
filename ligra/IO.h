// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#pragma once
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <sys/mman.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "parallel.h"
#include "blockRadixSort.h"
#include "quickSort.h"
#include "utils.h"
#include "graph.h"
#include "utils.h"
 #include <sys/mman.h>
using namespace std;

typedef pair<uintE, uintE> intPair;
typedef pair<uintE, pair<uintE, intE>> intTriple;

template <class E>
struct pairFirstCmp
{
  bool operator()(pair<uintE, E> a, pair<uintE, E> b)
  {
    return a.first < b.first;
  }
};

template <class E>
struct getFirst
{
  uintE operator()(pair<uintE, E> a) { return a.first; }
};

template <class IntType>
struct pairBothCmp
{
  bool operator()(pair<uintE, IntType> a, pair<uintE, IntType> b)
  {
    if (a.first != b.first)
      return a.first < b.first;
    return a.second < b.second;
  }
};

// A structure that keeps a sequence of strings all allocated from
// the same block of memory
struct words
{
  long n;         // total number of characters
  char *Chars;    // array storing all strings
  long m;         // number of substrings
  char **Strings; // pointers to strings (all should be null terminated)
  words() {}
  words(char *C, long nn, char **S, long mm)
      : Chars(C), n(nn), Strings(S), m(mm) {}
  void del()
  {
    free(Chars);
    free(Strings);
  }
};

inline bool isSpace(char c)
{
  switch (c)
  {
  case '\r':
  case '\t':
  case '\n':
  case 0:
  case ' ':
    return true;
  default:
    return false;
  }
}

_seq<char> mmapStringFromFile(const char *filename)
{
  struct stat sb;
  int fd = open(filename, O_RDONLY);
  if (fd == -1)
  {
    perror("open");
    exit(-1);
  }
  if (fstat(fd, &sb) == -1)
  {
    perror("fstat");
    exit(-1);
  }
  if (!S_ISREG(sb.st_mode))
  {
    perror("not a file\n");
    exit(-1);
  }
  char *p = static_cast<char *>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  if (p == MAP_FAILED)
  {
    perror("mmap");
    exit(-1);
  }
  if (close(fd) == -1)
  {
    perror("close");
    exit(-1);
  }
  size_t n = sb.st_size;
  //  char *bytes = newA(char, n);
  //  parallel_for(size_t i=0; i<n; i++) {
  //    bytes[i] = p[i];
  //  }
  //  if (munmap(p, sb.st_size) == -1) {
  //    perror("munmap");
  //    exit(-1);
  //  }
  //  cout << "mmapped" << endl;
  //  free(bytes);
  //  exit(0);
  return _seq<char>(p, n);
}

_seq<char> readStringFromFile(char *fileName)
{
  ifstream file(fileName, ios::in | ios::binary | ios::ate);
  if (!file.is_open())
  {
    std::cout << "Unable to open file: " << fileName << std::endl;
    abort();
  }
  long end = file.tellg();
  file.seekg(0, ios::beg);
  long n = end - file.tellg();
  char *bytes = newA(char, n + 1);
  file.read(bytes, n);
  file.close();
  return _seq<char>(bytes, n);
}

// parallel code for converting a string to words
words stringToWords(char *Str, long n)
{
  //这里还要把空格摘出来- -!，空格会被置零。
  {
    parallel_for(long i = 0; i < n; i++) if (isSpace(Str[i])) Str[i] = 0;
  }

  // mark start of words

  bool *FL = newA(bool, n);
  FL[0] = Str[0];
  {
    //从第一个字符开始搜索，FL==true(前一个字符[i-1]为0，字符[i]不为0)
    parallel_for(long i = 1; i < n; i++) FL[i] = Str[i] && !Str[i - 1];
  }

  // offset for each start of word
  _seq<long> Off = sequence::packIndex<long>(FL, n);
  long m = Off.n;
  long *offsets = Off.A;

  // pointer to each start of word
  char **SA = newA(char *, m);
  {
    parallel_for(long j = 0; j < m; j++) SA[j] = Str + offsets[j];
  }

  free(offsets);
  free(FL);
  return words(Str, n, SA, m);
}

template <class vertex>
graph<vertex> readGraphFromFile(char *fname, bool isSymmetric, bool mmap)
{
  words W;
  if (mmap)
  {
    _seq<char> S = mmapStringFromFile(fname);
    char *bytes = newA(char, S.n);
    // Cannot mutate the graph unless we copy.
    parallel_for(size_t i = 0; i < S.n; i++)
    {
      bytes[i] = S.A[i];
    }
    if (munmap(S.A, S.n) == -1)
    {
      perror("munmap");
      exit(-1);
    }
    S.A = bytes;
    W = stringToWords(S.A, S.n);
  }
  else
  {

    _seq<char> S = readStringFromFile(fname);
    W = stringToWords(S.A, S.n);
  }
#ifndef WEIGHTED
  if (W.Strings[0] != (string) "AdjacencyGraph")
  {
#else
  if (W.Strings[0] != (string) "WeightedAdjacencyGraph")
  {
#endif
    cout << "Bad input file" << endl;
    abort();
  }

  long len = W.m - 1;
  long n = atol(W.Strings[1]);
  long m = atol(W.Strings[2]);
#ifndef WEIGHTED
  if (len != n + m + 2)
  {
#else
  if (len != n + 2 * m + 2)
  {
#endif
    cout << "Bad input file" << endl;
    abort();
  }

  uintT *offsets = newA(uintT, n);
#ifndef WEIGHTED
  uintE *edges = newA(uintE, m);
  memset(edges,0,sizeof(uintE)*m);
  unsigned long start_addr = (unsigned long)(void*)edges ;
  for(;start_addr<start_addr+sizeof(uintE)*m;start_addr+=4096){
  	pbbs::print_addr("edges",start_addr);
  }
  
#else
  intE *edges = newA(intE, 2 * m);
#endif

  {
    parallel_for(long i = 0; i < n; i++) offsets[i] = atol(W.Strings[i + 3]);
  }
  {
    parallel_for(long i = 0; i < m; i++)
    {
#ifndef WEIGHTED
      edges[i] = atol(W.Strings[i + n + 3]);
#else
      edges[2 * i] = atol(W.Strings[i + n + 3]);
      edges[2 * i + 1] = atol(W.Strings[i + n + m + 3]);
#endif
    }
  }
  // W.del(); // to deal with performance bug in malloc

  vertex *v = newA(vertex, n);

  {
    parallel_for(uintT i = 0; i < n; i++)
    {
      uintT o = offsets[i];
      uintT l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
      v[i].setOutDegree(l);
#ifndef WEIGHTED
      v[i].setOutNeighbors(edges + o);
#else
      v[i].setOutNeighbors(edges + 2 * o);
#endif
    }
  }

  if (!isSymmetric)
  {
    uintT *tOffsets = newA(uintT, n);
    {
      parallel_for(long i = 0; i < n; i++) tOffsets[i] = INT_T_MAX;
    }
#ifndef WEIGHTED
    intPair *temp = newA(intPair, m);
#else
    intTriple *temp = newA(intTriple, m);
#endif
    {
      parallel_for(long i = 0; i < n; i++)
      {
        uintT o = offsets[i];
        for (uintT j = 0; j < v[i].getOutDegree(); j++)
        {
#ifndef WEIGHTED
          temp[o + j] = make_pair(v[i].getOutNeighbor(j), i);
#else
          temp[o + j] = make_pair(v[i].getOutNeighbor(j), make_pair(i, v[i].getOutWeight(j)));
#endif
        }
      }
    }
    free(offsets);

#ifndef WEIGHTED
#ifndef LOWMEM
    intSort::iSort(temp, m, n + 1, getFirst<uintE>());
#else
    quickSort(temp, m, pairFirstCmp<uintE>());
#endif
#else
#ifndef LOWMEM
    intSort::iSort(temp, m, n + 1, getFirst<intPair>());
#else
    quickSort(temp, m, pairFirstCmp<intPair>());
#endif
#endif

    tOffsets[temp[0].first] = 0;
#ifndef WEIGHTED
    uintE *inEdges = newA(uintE, m);
    inEdges[0] = temp[0].second;
#else
    intE *inEdges = newA(intE, 2 * m);
    inEdges[0] = temp[0].second.first;
    inEdges[1] = temp[0].second.second;
#endif
    {
      parallel_for(long i = 1; i < m; i++)
      {
#ifndef WEIGHTED
        inEdges[i] = temp[i].second;
#else
        inEdges[2 * i] = temp[i].second.first;
        inEdges[2 * i + 1] = temp[i].second.second;
#endif
        if (temp[i].first != temp[i - 1].first)
        {
          tOffsets[temp[i].first] = i;
        }
      }
    }

    free(temp);

    // fill in offsets of degree 0 vertices by taking closest non-zero
    // offset to the right
    sequence::scanIBack(tOffsets, tOffsets, n, minF<uintT>(), (uintT)m);

    {
      parallel_for(long i = 0; i < n; i++)
      {
        uintT o = tOffsets[i];
        uintT l = ((i == n - 1) ? m : tOffsets[i + 1]) - tOffsets[i];
        v[i].setInDegree(l);
#ifndef WEIGHTED
        v[i].setInNeighbors(inEdges + o);
#else
        v[i].setInNeighbors(inEdges + 2 * o);
#endif
      }
    }

    free(tOffsets);
    Uncompressed_Mem<vertex> *mem = new Uncompressed_Mem<vertex>(v, n, m, edges, inEdges);
    return graph<vertex>(v, n, m, mem);
  }
  else
  {
    free(offsets);
    Uncompressed_Mem<vertex> *mem = new Uncompressed_Mem<vertex>(v, n, m, edges);
    return graph<vertex>(v, n, m, mem);
  }
}

//[BFS]读取binary文件
template <class vertex>
graph<vertex> readGraphFromBinary(char *iFile, bool isSymmetric)
{
	int ret = mlockall(MCL_CURRENT|MCL_FUTURE);
		if(ret!=0){
				printf("mlockall failed!\n");
		}
  //三个文件
  char *config = (char *)".config";
  char *adj = (char *)".adj";
  char *idx = (char *)".idx";
  char configFile[strlen(iFile) + strlen(config) + 1];
  char adjFile[strlen(iFile) + strlen(adj) + 1];
  char idxFile[strlen(iFile) + strlen(idx) + 1];
  *configFile = *adjFile = *idxFile = '\0';
  strcat(configFile, iFile);
  strcat(adjFile, iFile);
  strcat(idxFile, iFile);
  strcat(configFile, config);
  strcat(adjFile, adj);
  strcat(idxFile, idx);

  //.config文件只保存一个n
  ifstream in(configFile, ifstream::in);
  long n;
  in >> n;
  in.close();
  
  ifstream in2(adjFile, ifstream::in | ios::binary); // stored as uints
  in2.seekg(0, ios::end);
  long size = in2.tellg();
  in2.seekg(0);
#ifdef WEIGHTED
  long m = size / (2 * sizeof(uint));
#else
  // m记录了adj文件的uint的数量
  long m = size / sizeof(uint);
#endif
  //char *s = (char *)malloc(size);
  long page_to_allocate = size / 2097152 ;
  if(page_to_allocate==0)page_to_allocate=1;
  else{
  	page_to_allocate+=1;
  }
  long mmap_size = 2097152 * page_to_allocate;
  //printf("mmap_size is %ld",mmap_size);
  //printf("size is %ld",size);
  char *s = (char*)mmap(NULL, mmap_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB,-1, 0);
  if(s == MAP_FAILED){
	printf("mmap failed!\n")	;
	exit(-1);
  } 
  //printf("mmap success!\n");
  //char *s_end = s+size;
  memset(s,49,size);
  //printf("mmap memset\n");
  //char *p=s;
  unsigned long start_addr = (unsigned long)(void*)s ;
  unsigned long end_addr = (unsigned long)(void*)s+size;
  unsigned long p = start_addr;
  for(;p<end_addr;p+=2097152){
  	pbbs::print_addr("edges",(unsigned long)(void*)p);
  } 
  //printf("done\n");
  //pbbs::print_address("edges",start_addr,start_addr+size-4096);
  //printf("mmap address\n");
  //pbbs::print_address("edges range",start_addr,end_addr);
  in2.read(s, size);
  in2.close();
  //[内存]edges是邻接数组，例如: edges=[2,3,4,5,6], idx=[0,2,....]，则可以认为，vertex(0)的临点是2和3，vertex(1)的临点是4,5,6。
  uintE *edges = (uintE *)s;



  // idx文件保存的数据格式和n相同，并且数据类型是intT
  ifstream in3(idxFile, ifstream::in | ios::binary); // stored as longs
  in3.seekg(0, ios::end);
  size = in3.tellg();
  in3.seekg(0);
  if (n != size / sizeof(intT))
  {
    cout << "File size wrong\n";
    abort();
  }

  char *t = (char *)malloc(size);
  in3.read(t, size);
  in3.close();
  //这个offset保存了idx文件的所有信息，用一个连续数据保存，
  uintT *offsets = (uintT *)t;
  //[内存]这个点集需要监控的。
  //vertex *v = newA(vertex, n);
  int page_num = (n*sizeof(vertex))/2097152 + 1;
  long vertex_size = page_num * 2097152;
  vertex *v = (vertex*)mmap(NULL, vertex_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB,-1, 0);
  if(v == MAP_FAILED){
	printf("mmap failed!\n")	;
	exit(-1);
  } 
  memset(v,3,sizeof(vertex)*n);
  start_addr = (unsigned long)(void*)v ;
  end_addr = ((unsigned long)(void*)v) + sizeof(vertex)*n;
  p = start_addr;
  for(;p<end_addr;p+=2097152){
  	pbbs::print_addr("Vertices",(unsigned long)(void*)p);
  }
#ifdef WEIGHTED
  intE *edgesAndWeights = newA(intE, 2 * m);
  {
    parallel_for(long i = 0; i < m; i++)
    {
      edgesAndWeights[2 * i] = edges[i];
      edgesAndWeights[2 * i + 1] = edges[i + m];
    }
  }
  // free(edges);
#endif
  {
    //这个结构我们在graph500见过，idx表示点vertex_i，offset(也就是adj结构)的offset[i+1]-offset[i]里保存的都是这个vertex_i的临点。
    parallel_for(long i = 0; i < n; i++)
    {
      uintT o = offsets[i];
      uintT l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
      //通过edge数组的起始地址和临边数量可以直接在edges数组里找到某个点的所有临点。
      v[i].setOutDegree(l);
#ifndef WEIGHTED
      v[i].setOutNeighbors((uintE *)edges + o);
#else
      v[i].setOutNeighbors(edgesAndWeights + 2 * o);
#endif
    }
  }
  //[BFS]无向图都是Symmetric，所以我们的BFS不走这个逻辑
  if (!isSymmetric)
  {
    uintT *tOffsets = newA(uintT, n);
    {
      parallel_for(long i = 0; i < n; i++) tOffsets[i] = INT_T_MAX;
    }
#ifndef WEIGHTED
    intPair *temp = newA(intPair, m);
#else
    intTriple *temp = newA(intTriple, m);
#endif
    {
      parallel_for(intT i = 0; i < n; i++)
      {
        uintT o = offsets[i];
        for (uintT j = 0; j < v[i].getOutDegree(); j++)
        {
#ifndef WEIGHTED
          temp[o + j] = make_pair(v[i].getOutNeighbor(j), i);
#else
          temp[o + j] = make_pair(v[i].getOutNeighbor(j), make_pair(i, v[i].getOutWeight(j)));
#endif
        }
      }
    }
    free(offsets);
#ifndef WEIGHTED
#ifndef LOWMEM
    intSort::iSort(temp, m, n + 1, getFirst<uintE>());
#else
    quickSort(temp, m, pairFirstCmp<uintE>());
#endif
#else
#ifndef LOWMEM
    intSort::iSort(temp, m, n + 1, getFirst<intPair>());
#else
    quickSort(temp, m, pairFirstCmp<intPair>());
#endif
#endif
    tOffsets[temp[0].first] = 0;
#ifndef WEIGHTED
    uintE *inEdges = newA(uintE, m);
    inEdges[0] = temp[0].second;
#else
    intE *inEdges = newA(intE, 2 * m);
    inEdges[0] = temp[0].second.first;
    inEdges[1] = temp[0].second.second;
#endif
    {
      parallel_for(long i = 1; i < m; i++)
      {
#ifndef WEIGHTED
        inEdges[i] = temp[i].second;
#else
        inEdges[2 * i] = temp[i].second.first;
        inEdges[2 * i + 1] = temp[i].second.second;
#endif
        if (temp[i].first != temp[i - 1].first)
        {
          tOffsets[temp[i].first] = i;
        }
      }
    }
    free(temp);
    // fill in offsets of degree 0 vertices by taking closest non-zero
    // offset to the right
    sequence::scanIBack(tOffsets, tOffsets, n, minF<uintT>(), (uintT)m);
    {
      parallel_for(long i = 0; i < n; i++)
      {
        uintT o = tOffsets[i];
        uintT l = ((i == n - 1) ? m : tOffsets[i + 1]) - tOffsets[i];
        v[i].setInDegree(l);
#ifndef WEIGHTED
        v[i].setInNeighbors((uintE *)inEdges + o);
#else
        v[i].setInNeighbors((intE *)(inEdges + 2 * o));
#endif
      }
    }
    free(tOffsets);
#ifndef WEIGHTED
    Uncompressed_Mem<vertex> *mem = new Uncompressed_Mem<vertex>(v, n, m, edges, inEdges);
    return graph<vertex>(v, n, m, mem);
#else
    Uncompressed_Mem<vertex> *mem = new Uncompressed_Mem<vertex>(v, n, m, edgesAndWeights, inEdges);
    return graph<vertex>(v, n, m, mem);
#endif
  }
  free(offsets);
  printf("before return?\n");
#ifndef WEIGHTED
  //v是一个连续内存空间，记录点集；n是点的数量；m是adj数组的长度，edges就是adj数组。这两个结构都是以指针的形势传入，所以监控内存可以不用在里面监控
  Uncompressed_Mem<vertex> *mem = new Uncompressed_Mem<vertex>(v, n, m, edges);
  return graph<vertex>(v, n, m, mem);
#else
  Uncompressed_Mem<vertex> *mem = new Uncompressed_Mem<vertex>(v, n, m, edgesAndWeights);
  return graph<vertex>(v, n, m, mem);
#endif
}

//[BFS]从这里读取真个Graph的信息
template <class vertex>
graph<vertex> readGraph(char *iFile, bool compressed, bool symmetric, bool binary, bool mmap)
{
  //可以以binary或者file方式
  if (binary)
    return readGraphFromBinary<vertex>(iFile, symmetric);
  else
    return readGraphFromFile<vertex>(iFile, symmetric, mmap);
}

template <class vertex>
graph<vertex> readCompressedGraph(char *fname, bool isSymmetric, bool mmap)
{
  char *s;
  //这是使用mmap模式的时候，目前没用到
  if (mmap)
  {
    _seq<char> S = mmapStringFromFile(fname);
    // Cannot mutate graph unless we copy.
    char *bytes = newA(char, S.n);
    parallel_for(size_t i = 0; i < S.n; i++)
    {
      bytes[i] = S.A[i];
    }
    if (munmap(S.A, S.n) == -1)
    {
      perror("munmap");
      exit(-1);
    }
    s = bytes;
  }
  else
  {
    ifstream in(fname, ifstream::in | ios::binary);
    in.seekg(0, ios::end);
    long size = in.tellg();
    in.seekg(0);
    cout << "size = " << size << endl;
    s = (char *)malloc(size);
    in.read(s, size);
    in.close();
  }

  long *sizes = (long *)s;
  long n = sizes[0], m = sizes[1], totalSpace = sizes[2];

  cout << "n = " << n << " m = " << m << " totalSpace = " << totalSpace << endl;
  cout << "reading file..." << endl;

  uintT *offsets = (uintT *)(s + 3 * sizeof(long));
  long skip = 3 * sizeof(long) + (n + 1) * sizeof(intT);
  uintE *Degrees = (uintE *)(s + skip);
  skip += n * sizeof(intE);
  uchar *edges = (uchar *)(s + skip);

  uintT *inOffsets;
  uchar *inEdges;
  uintE *inDegrees;
  if (!isSymmetric)
  {
    skip += totalSpace;
    uchar *inData = (uchar *)(s + skip);
    sizes = (long *)inData;
    long inTotalSpace = sizes[0];
    cout << "inTotalSpace = " << inTotalSpace << endl;
    skip += sizeof(long);
    inOffsets = (uintT *)(s + skip);
    skip += (n + 1) * sizeof(uintT);
    inDegrees = (uintE *)(s + skip);
    skip += n * sizeof(uintE);
    inEdges = (uchar *)(s + skip);
  }
  else
  {
    inOffsets = offsets;
    inEdges = edges;
    inDegrees = Degrees;
  }

  vertex *V = newA(vertex, n);
  parallel_for(long i = 0; i < n; i++)
  {
    long o = offsets[i];
    uintT d = Degrees[i];
    V[i].setOutDegree(d);
    V[i].setOutNeighbors(edges + o);
  }

  if (sizeof(vertex) == sizeof(compressedAsymmetricVertex))
  {
    parallel_for(long i = 0; i < n; i++)
    {
      long o = inOffsets[i];
      uintT d = inDegrees[i];
      V[i].setInDegree(d);
      V[i].setInNeighbors(inEdges + o);
    }
  }

  cout << "creating graph..." << endl;
  Compressed_Mem<vertex> *mem = new Compressed_Mem<vertex>(V, s);

  graph<vertex> G(V, n, m, mem);
  return G;
}
