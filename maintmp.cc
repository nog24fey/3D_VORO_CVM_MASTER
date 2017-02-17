#include "./src/directorymake.h"
#include <string>
#include <iostream>
using namespace std;

int main(int argc, char **argv) {
  directoryMake(argv[1]);
  string av1(argv[1]);
  cout<<av1<<endl;
  return 1;
}
