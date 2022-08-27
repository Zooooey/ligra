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
#include <vector>

typedef double fType;

struct BC_F {
  fType* NumPaths;
  bool* Visited;
  BC_F(fType* _NumPaths, bool* _Visited) : 
    NumPaths(_NumPaths), Visited(_Visited) {}
  inline bool update(uintE s, uintE d){ //Update function for forward phase
    fType oldV = NumPaths[d];
    NumPaths[d] += NumPaths[s];
    //这里返回true的时候(处理前NumPaths[d]==0)，点d才会考虑进入下一次迭代的vertexSubset里。
    return oldV == 0.0;
  }
  inline bool updateAtomic (uintE s, uintE d) { //atomic Update, basically an add
    volatile fType oldV, newV; 
    do { 
      oldV = NumPaths[d]; newV = oldV + NumPaths[s];
    } while(!CAS(&NumPaths[d],oldV,newV));
    return oldV == 0.0;
  }
  //条件是Visited[D] == 0说明框架只处理Visited == 0的点，也就是没有被处理过的点。
  inline bool cond (uintE d) { return Visited[d] == 0; } //check if visited
};

struct BC_Back_F {
  fType* Dependencies;
  bool* Visited;
  BC_Back_F(fType* _Dependencies, bool* _Visited) : 
    Dependencies(_Dependencies), Visited(_Visited) {}
  inline bool update(uintE s, uintE d){ //Update function for backwards phase
    fType oldV = Dependencies[d];
    Dependencies[d] += Dependencies[s];
    return oldV == 0.0;
  }
  inline bool updateAtomic (uintE s, uintE d) { //atomic Update
    volatile fType oldV, newV;
    do {
      oldV = Dependencies[d];
      newV = oldV + Dependencies[s];
    } while(!CAS(&Dependencies[d],oldV,newV));
    return oldV == 0.0;
  }
  //这里还是一样，没有被访问过的点才会进入下一次迭代
  inline bool cond (uintE d) { return Visited[d] == 0; } //check if visited
};

//vertex map function to mark visited vertexSubset
struct BC_Vertex_F {
  bool* Visited;
  BC_Vertex_F(bool* _Visited) : Visited(_Visited) {}
  inline bool operator() (uintE i) {
    Visited[i] = 1;
    return 1;
  }
};

//vertex map function (used on backwards phase) to mark visited vertexSubset
//and add to Dependencies score
struct BC_Back_Vertex_F {
  bool* Visited;
  fType* Dependencies, *inverseNumPaths;
  BC_Back_Vertex_F(bool* _Visited, fType* _Dependencies, fType* _inverseNumPaths) : 
    Visited(_Visited), Dependencies(_Dependencies), inverseNumPaths(_inverseNumPaths) {}
  inline bool operator() (uintE i) {
    Visited[i] = 1;
    Dependencies[i] += inverseNumPaths[i];
    return 1; }};

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  long start = P.getOptionLongValue("-r",0);
  long n = GA.n;
	printf("???\n");
  //初始化一个连续空间NumPaths，并全部置0.0;
  static fType* NumPaths = NULL;
  if(NumPaths == NULL){
    NumPaths = pbbs::mmap_huge_page(sizeof(fType)*n, "NumPaths");
    //NumPaths = newA(fType,n);
    {parallel_for(long i=0;i<n;i++) NumPaths[i] = 0.0;}
    pbbs::print_hugepage_addr("NumPaths", (unsigned long)(void*)NumPaths, (unsigned long)(void*)(NumPaths + n));
  }else {
    {parallel_for(long i=0;i<n;i++) NumPaths[i] = 0.0;}
  }
  //起点置1.0
  NumPaths[start] = 1.0;

  //visited结构全部初始化0，起点置1
  static bool* Visited = NULL;
  if(Visited == NULL){
    Visited = pbbs::mmap_huge_page(sizeof(bool)*n,"Visited");
    //Visited = newA(bool,n);
    {parallel_for(long i=0;i<n;i++) Visited[i] = 0;}
    pbbs::print_hugepage_addr("Visited", (unsigned long)(void*)Visited, (unsigned long)(void*)(Visited + n));
  } else {
    {parallel_for(long i=0;i<n;i++) Visited[i] = 0;}
  }
  
  Visited[start] = 1;

  //构建一个vertexSubset，只有一个起点start。
  vertexSubset Frontier(n,start);
 
  //分层？
  vector<vertexSubset> Levels;
  Levels.push_back(Frontier);

  long round = 0;
  while(!Frontier.isEmpty()){ //first phase
    round++;
    //第一次返回的点集是起点的临边所指向的终点，这些终点的NumPaths被加上起点的Numapths值
    vertexSubset output = edgeMap(GA, Frontier, BC_F(NumPaths,Visited));
    //处理output里的点，第一次处理的话就是起点的所有后继节点，其Visited标记为1
    vertexMap(output, BC_Vertex_F(Visited)); //mark visited
    //推入Levels保存起当前output，
    Levels.push_back(output); //save frontier onto Levels
    //Frontier变为刚才的点，然后继续。
    Frontier = output;
    //到这里循环结束，如此往复之后，被遍历的点都是1(只有一种情况下其Numpaths值大于1，就是那些被同一个层级的点的出边指向的点。例如：vertexSubset={1,2},edges={1->3,2->3},这种情况下Numpaths[3]=2,因为Visited[3]=1发生在update(1,3)和update(2,3)后)
  }

  //又是一个连续结构，初始化置0
  static fType* Dependencies = NULL;
  if(Dependencies == NULL){
    Dependencies = pbbs::mmap_huge_page(sizeof(fType)*n,"Dependencies");
    //Dependencies = newA(fType,n);
    {parallel_for(long i=0;i<n;i++) Dependencies[i] = 0.0;}
    pbbs::print_hugepage_addr("Dependencies", (unsigned long)(void*)Dependencies, (unsigned long)(void*)(Dependencies + n));
  } else {
    {parallel_for(long i=0;i<n;i++) Dependencies[i] = 0.0;}
  }
  
  


  //刚才的NumPaths换了个指针变成inverse?
  //invert numpaths
  fType* inverseNumPaths = NumPaths;
  //值取倒数了。
  {parallel_for(long i=0;i<n;i++) inverseNumPaths[i] = 1/inverseNumPaths[i];}

  //????如果只迭代一次，round=1，只有Levels[0]有值吧？
  Levels[round].del();
  //reuse Visited
  //这个数据结构直接重用。全部置0。
  {parallel_for(long i=0;i<n;i++) Visited[i] = 0;}
  //vertexSubset等于最后一轮的点。
  Frontier = Levels[round-1];
  //使用vertex处理点集而非边集。
  //1. 这个Dependencies初始化全是0
  //2. 这个inverseNumPaths全是第一次正向迭代结果的倒数
  //3. Dependencies运行结束后得到inverseNumPaths里相应的值(虽然BC_Back_Vertex_f执行的是累加，但是Dependencies在这里的时候全是0)
  vertexMap(Frontier,BC_Back_Vertex_F(Visited,Dependencies,inverseNumPaths));

  //tranpose graph
  
  //出边和入边交换，倒转整个graph
  GA.transpose();
  for(long r=round-2;r>=0;r--) { //backwards phase
    //拿着刚才的倒数，倒着遍历一次累加回去。
    edgeMap(GA, Frontier, BC_Back_F(Dependencies,Visited), -1, no_output);
    //Levels被当成栈了，这里在出栈。
    Frontier.del();
    Frontier = Levels[r]; //gets frontier from Levels array
    //vertex map to mark visited and update Dependencies scores
    //标记Visited
    vertexMap(Frontier,BC_Back_Vertex_F(Visited,Dependencies,inverseNumPaths));
  }
  
  Frontier.del();

  //Update dependencies scores
  //特定公式，得到最终值
  parallel_for(long i=0;i<n;i++) {
    Dependencies[i]=(Dependencies[i]-inverseNumPaths[i])/inverseNumPaths[i];
  }
  //free(inverseNumPaths);
  //free(Visited);
  //free(Dependencies);
  munmap(inverseNumPaths,n*sizeof(fType));
  munmap(Visited,n*sizeof(bool));
  munmap(Dependencies,n*sizeof(fType));
}
