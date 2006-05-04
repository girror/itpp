#include <itpp/itcomm.h>

using std::cout;
using std::endl;
using namespace itpp;

int main()
{

 Real_Timer time;
 int N=10000, M=8, k=needed_bits(M - 1);
 
 vec EbN0 = "0:2:20";
 vec EN0ggr = inv_dB(EbN0)*k;
 vec pb(EbN0.size());

 PAM pam(M);
 Array<BERC> berc(EbN0.size());
 
 bvec transmitted_bits, received_bits;
 vec transmitted_symbols, received_symbols, noise;

 transmitted_bits = randb(N*k);

 pam.modulate_bits(transmitted_bits, transmitted_symbols);

 cout << "Energy=" << sum(sqr(transmitted_symbols))/double(transmitted_symbols.size()) << endl;
 randn(transmitted_symbols.size(), noise);

 for (int i=0; i<EbN0.size(); i++) {
   received_symbols = transmitted_symbols + sqrt(1.0/EN0ggr(i)/2.0)*noise;    
   pam.demodulate_bits(received_symbols, received_bits);
   berc(i).count(transmitted_bits, received_bits);
   pb(i) = berc(i).get_errorrate();
 }
 cout << "SNR=" << EbN0 << endl;
 cout << "pb=" << pb << endl;

} 

