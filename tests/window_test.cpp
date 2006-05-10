#include <itpp/itbase.h>

using std::cout;
using std::endl;
using namespace itpp;

int main(void)
{

 cout << "================================" << endl;
 cout << "    Test of window functions " << endl;
 cout << "================================" << endl;

 cout << "hamming(32) = " << round_to_zero(hamming(32)) << endl;
 cout << "hamming(128) = " << round_to_zero(hamming(128)) << endl;

 cout << "hanning(32) = " << round_to_zero(hanning(32)) << endl;
 cout << "hanning(128) = " << round_to_zero(hanning(128)) << endl;

 cout << "hann(32) = " << round_to_zero(hann(32)) << endl;
 cout << "hann(128) = " << round_to_zero(hann(128)) << endl;

 cout << "blackman(32) = " << round_to_zero(blackman(32)) << endl;
 cout << "blackman(128) = " << round_to_zero(blackman(128)) << endl;

 cout << "triang(32) = " << round_to_zero(triang(32)) << endl;
 cout << "triang(128) = " << round_to_zero(triang(128)) << endl;

 return 0;

}
