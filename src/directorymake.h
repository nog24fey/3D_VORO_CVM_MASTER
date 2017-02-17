#ifndef DIRECTORYMAKE_H_
#define DIRECTORYMAKE_H_ 1

#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>

void directoryMake(char* dirname) {
  /*0->success*/
  struct stat st;
  if(stat(dirname,&st) != 0){
    if (mkdir(dirname, S_IRUSR | S_IWUSR | S_IXUSR |
              S_IRGRP | S_IWGRP | S_IXGRP |
              S_IROTH | S_IXOTH | S_IXOTH) == 0) {
      printf("%s is created\n", dirname);
    } else {
      printf("%s cannot be created\n", dirname);
      perror("");
      exit(1);
    }
  }
}

#endif



