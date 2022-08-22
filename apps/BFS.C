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
#include "ligra.h"

struct BFS_F
{
  uintE *Parents; 
  //vertexSubset *Frontier  = NULL;
  BFS_F(uintE *_Parents) : Parents(_Parents) {}
  inline bool update(uintE s, uintE d)
  { // Update
    if (Parents[d] == UINT_E_MAX)
    {
      Parents[d] = s;
      return 1;
    }
    else
      return 0;
  }
  inline bool updateAtomic(uintE s, uintE d)
  { // atomic version of Update
    return (CAS(&Parents[d], UINT_E_MAX, s));
  }
  // cond function checks if vertex has been visited yet
  inline bool cond(uintE d) { return (Parents[d] == UINT_E_MAX); }
};

template <class vertex>
void Compute(graph<vertex> &GA, commandLine P)
{
  long start = P.getOptionLongValue("-r", 0);
  //应该是Graph的数量
  long n = GA.n;
  // GA.print_address();
  // creates Parents array, initialized to all -1, except for start
  //[内存]new A是对malloc的一个别名，这里创建了一个连续内存空间。
  // uintE* Parents = newA(uintE,n);
  static uintE *Parents = MAP_FAILED;
  static vertexSubset* Frontier = NULL;
  
  if (Parents == MAP_FAILED)
  {
    int page_num = (n * sizeof(uintE)) / 2097152 + 1;
    long parent_size = page_num * 2097152;
    uintE *Parents = (uintE *)mmap(NULL, parent_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
    if (Parents == MAP_FAILED)
    {
      printf("mmap failed!\n");
      exit(-1);
    }
    parallel_for(long i = 0; i < n; i++) Parents[i] = UINT_E_MAX;
    unsigned long start_addr = (unsigned long)(void *)Parents;
    unsigned long end_addr = ((unsigned long)(void *)Parents) + sizeof(uintE) * n;
    unsigned long p = start_addr;
    for (; p < end_addr; p += 2097152)
    {
      char c = *((char *)p);
      pbbs::print_addr("Parents", p);
    }
  }
  else
  {
    parallel_for(long i = 0; i < n; i++) Parents[i] = UINT_E_MAX;
  }

  Parents[start] = start;

  //以单个vertex初始化一个vertexSubset
  if (Frontier == NULL)
  {
    Frontier = new vertexSubset(n, start); // creates initial frontier
    // uintE* S = newA(uintE,n);
    int S_page_num = (n * sizeof(uintE)) / 2097152 + 1;
    long S_size = S_page_num * 2097152;
    uintE *S = (uintE *)mmap(NULL, S_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
    if (S == MAP_FAILED)
    {
      printf("mmap failed!\n");
      exit(-1);
    }
    memset(S, 0, sizeof(uintE) * n);
    unsigned long start_addr = (unsigned long)(void *)S;
    unsigned long end_addr = ((unsigned long)(void *)S) + sizeof(uintE) * n;
    unsigned long p = start_addr;
    for (; p < end_addr; p += 2097152)
    {
      char c = *((char *)p);
      pbbs::print_addr("Frontier vertices", p);
    }

    // bool* D = newA(bool,n);
    int D_page_num = (n * sizeof(bool)) / 2097152 + 1;
    long D_size = D_page_num * 2097152;
    bool *D = (bool *)mmap(NULL, D_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
    if (D == MAP_FAILED)
    {
      printf("mmap failed!\n");
      exit(-1);
    }
    memset(D, 0, sizeof(bool) * n);
    start_addr = (unsigned long)(void *)D;
    end_addr = ((unsigned long)(void *)D) + sizeof(bool) * n;
    p = start_addr;
    for (; p < end_addr; p += 2097152)
    {
      char c = *((char *)p);
      pbbs::print_addr("Frontier bool", p);
    }
    // pbbs::print_address("Fontier vertices", (unsigned long)(void*)S, (unsigned long)(void*)(S + n));
    // pbbs::print_address("Fontier bool array", (unsigned long)(void*)D, (unsigned long)(void*)(D + n));
    Frontier->setD(D, sizeof(bool) * n);
    Frontier->setS(S, sizeof(uintE) * n);
    Frontier->setShouldFree(false);
  }
  else
  {
    Frontier->refresh(n, start);
    uintE *S = Frontier->getS();
    bool* D = Frontier->getD();
    memset(S, 0, sizeof(uintE) * n);
    memset(D, 0, sizeof(bool) * n);
  }

  int times = 0;
  while (!Frontier.isEmpty())
  { // loop until frontier is empty
    vertexSubset output = bfsEdgeMap(GA, Frontier, BFS_F(Parents));
    // printf("return from bfsEdgeMap!\n");
    Frontier.del();
    // printf("Frontier.del()!\n");
    Frontier = output; // set new frontier
    // printf("=output\n");
    times++;
    printf("round :%d\n", times);
  }
  printf("pre end\n");
  Frontier.del();
  printf("post end\n");
  // free(Parents);
  // free(S);
  // free(D);
}
