#include <cstdio>
#include <sys/stat.h>

int directoryMake(char* dirname){
  /*0->success*/
  int return_code = 0;
  struct stat st;
  if(stat(dirname,&st) != 0){
    if (mkdir(dirname, S_IRUSR | S_IWUSR | S_IXUSR |
              S_IRGRP | S_IWGRP | S_IXGRP |
              S_IROTH | S_IXOTH | S_IXOTH) == 0) {
      printf("%s is created\n", dirname);
    } else {
      printf("%s cannot be created\n", dirname);
      perror("");
      return_code = 1;
    }
  }
  return return_code;
}



