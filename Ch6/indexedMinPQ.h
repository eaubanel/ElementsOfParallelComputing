// Indexed min priority queue.
// Implementation can be modified from Sedgewick's 
// Algorithms in C, 3rd edition.
#ifndef INDEXEDMINPQ_H
#define INDEXEDMINPQ_H
// initialize heap of max size m, 
// which indexes into items array of lenth n
void minHeapInit(int *items, int n, int m);
// returns 1 if heap is empty, 0 otherwise
int minHeapEmpty();
// return 1 if heap is full, 0 otherwise
int minHeapFull();
// inserts element at index k into heap
void minHeapInsert(int k);
// extract index at top of heap
int heapExtractMin();
// change value of priority at index k then re-heapify
void minHeapChange(int k);
// returns 1 if element at index k is in heap, 0 otherwise
int isInHeap(int k);
// frees memory associated with heap
void minHeapDestroy();
#endif
