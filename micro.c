#include<sys/time.h>
#include<stdio.h>
#include <unistd.h>
#include <string.h>
//#include<math.h>

double amicro_() 
{
    struct timeval tv;
    
    gettimeofday(&tv,NULL);
    
    return ((double)(tv.tv_sec) + 1e-6*((double)(tv.tv_usec)));
}

void gethostname_(char (*host)[20])
{
     char hostname[1024];
     int i;

//     gethostname(hostname, 1023);
//     printf("hostname %s %d \n",hostname,strlen(hostname));
     hostname[strlen(hostname)] = 0;
     strcpy(*host,hostname);
     for(i = strlen(hostname);i <= 20;i++)
     {
        (*host)[i] = ' ';
     }


}


/*
int main(void)
{
 int i;
 printf("%20.15f \n",amicro_()); 
 for(i = 1;i < 100000;i++) sin(i);
 printf("%20.15f \n",amicro_());  
 char hn[20];

 gethostname_(&hn);
 printf("Hostname: %s\n", hn);

 
 return 0;
}
*/
