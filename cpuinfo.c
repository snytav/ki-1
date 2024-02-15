#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "lspci.h"

void writecpuinfo_(int *rank)
{
   FILE *cpuinfo = fopen("/proc/cpuinfo", "rb");
   FILE *f;
   char str[100];
   sprintf(str,"cpuinfo%05d.dat",*rank);
   f  = fopen(str,"wt");
   char *arg = 0;
   size_t size = 0;
   while(getdelim(&arg, &size, 0, cpuinfo) != -1)
   {
      fputs(arg,f);
   }
   char hostname[1024];
//   gethostname(hostname, 1023);
   fprintf(f,"hostname %s",hostname);

  
   free(arg);
   fclose(cpuinfo);
   fclose(f);

   lspci(*rank);
   meminfo(*rank);
}



/*int main(int argc, char **argv)
{
   int rank = 0;

   writecpuinfo(&rank);

   return 0;
}*/
