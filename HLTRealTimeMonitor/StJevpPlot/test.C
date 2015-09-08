#include <stdio.h>

class a
{
public:

  void test() {
    printf("Hello\n");
  }
};

class b : public a
{
public:
  void test() {
    a::test();
    printf("There\n");
  }
};

int main(int argc, char *argv[])
{
  a aa;
  b bb;
  
  aa.test();
  printf("---\n");
  bb.test();
  
  printf("---\n");
  a *p = (a *)&bb;
  p->test();
}
