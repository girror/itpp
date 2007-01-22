#include <itpp/itbase.h>

using namespace itpp;
using namespace std;

#ifndef ITFILE_TEST_FILE

int main() {
  cout << "ITFILE_TEST_FILE not defined. Test skipped." << endl;
}

#else

int main()
{
  int x = 1234567890;
  int y;

  it_file ff;
  ff.open(string(ITFILE_TEST_FILE));
#ifdef SAVE_DATA  
  ff << Name("x") << x;
#endif
  ff >> Name("x") >> y;
  ff.close();

  cout << x << endl;
  cout << y << endl;

  return 0;
}

#endif
