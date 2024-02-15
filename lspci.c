#include <stdio.h>
#include <stdlib.h>


int lspci(int rank)
{

  FILE *fp,*f;
  char path[1035];
  char fname[100];
  //int rank = 13;

  /* Open the command for reading. */
  fp = popen("/bin/lspci", "r");
  if (fp == NULL) {
    printf("Failed to run command\n" );
    exit(1);
  }

  sprintf(fname,"ls_pci%010d.dat",rank);

  if((f = fopen(fname,"wt")) == NULL) return 1;

  /* Read the output a line at a time - output it. */
  while (fgets(path, sizeof(path), fp) != NULL) {
    fprintf(f,"%s", path);
  }

  /* close */
  pclose(fp);

  return 0;
}

int meminfo(int rank)
{

  FILE *fp,*f;
  char path[1035];
  char fname[100];
  //int rank = 13;

  /* Open the command for reading. */
  fp = popen("cat /proc/meminfo", "r");
  if (fp == NULL) {
    printf("Failed to run command\n" );
    exit(1);
  }

  sprintf(fname,"meminfo%010d.dat",rank);

  if((f = fopen(fname,"wt")) == NULL) return 1;

  /* Read the output a line at a time - output it. */
  while (fgets(path, sizeof(path), fp) != NULL) {
    fprintf(f,"%s", path);
  }

  /* close */
  pclose(fp);

  return 0;
}

/*
int main(void)
{
    lspci(13);
    meminfo(67);
    return 0; 
}
*/


