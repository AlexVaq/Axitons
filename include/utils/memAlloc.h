#ifndef	MemAllocGuard
#define	MemAllocGuard

#include "enum-vars.h"

void	trackFree    (void *);
void	trackAlloc   (void **, size_t);
void	deviceAlloc  (void **, size_t);
void	printMemStats();
#endif

