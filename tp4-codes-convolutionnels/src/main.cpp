#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include <cstdlib>
#include <ctime>
#include <list>
#include <algorithm>


const int N = 2;
const int K = 1;
const int R = 4;

const int NB_ERRORS = 2; // Nb_errors to add in the transmission channel

#define DEBUG

using namespace std;

////////////////////////////////////////////////////////////
//                  Structure - Path                      //
//                                                        //
//               Represents a path in the tree            //
////////////////////////////////////////////////////////////

class Path
{
  public:
    bitset<R> reg;
    vector< bitset<K> > mess_dec;
    int cumulated_distance;
    Path(){};
    Path( bitset<R> a, int b ) : reg(a), cumulated_distance(b) {}
};

////////////////////////////////////////////////////////////
//             Comparation Weight of Hamming Distance     //
//                                                        //
// A function than return the lightest path               //
////////////////////////////////////////////////////////////

bool weight_comp(Path path_a, Path path_b)
{
  return ( path_a.cumulated_distance < path_b.cumulated_distance);
}

////////////////////////////////////////////////////////////
//                  Hamming Distance                      //
//                                                        //
// A function that calculate the Hamming distance between //
// two bitsets                                            //
////////////////////////////////////////////////////////////

int hamming_distance(bitset<N> bitset_1, bitset<N> bitset_2)
{
  int cpt = 0;
  if (bitset_1[0] != bitset_2[0])
    cpt++;
  if (bitset_1[1] != bitset_2[1])
    cpt++;
  return cpt;
}

bitset<N> code(const bitset<R+1>& reg)
{
  bitset<N> cod_out;
	int g0, g1;
  bitset<R+1> G0(25);
  bitset<R+1> G1(27);
	g0 = (reg & G0).count() % 2;
	g1 = (reg & G1).count() % 2;

	cod_out.set(0, g0);
	cod_out.set(1, g1);
  return cod_out;
}

////////////////////////////////////////////////////////////
//      template<int bits> bitset<bits> randBitset()      //
//                                                        //
//               Generate random bitset                   //
////////////////////////////////////////////////////////////

template <int bits>
bitset<bits> randBitset()
{
  bitset<bits> r(rand());
  for (int i = 0; i < bits / 16 - 1; i++)
  {
    r <<= 16;
    r |= bitset<bits>(rand());
  }
  return r;
}

////////////////////////////////////////////////////////////
// vector< bitset<N> > GSM_code(vector< bitset<K> > mess) //
//                                                        //
//     Convolutional coding of a message (GSM norm)       //
////////////////////////////////////////////////////////////

vector< bitset<N> > GSM_code(vector< bitset<K> > mess)
{
  int i = 0, g0, g1;
  vector< bitset<N> > mess_out;

  bitset<N> cod_out;
  bitset<R + 1> G0(25);
  bitset<R + 1> G1(27);
  bitset<R + 1> reg;
  reg.reset();

#ifdef DEBUG
  cout << "-------------------- Debug Informations (Coding) --------------------" << endl
       << endl;
  cout << "Initial register ( u(i-4)  u(i-3)  u(i-2)  u(i-1)  u(i)  ): " << reg << endl;
  cout << "Polynom G0       ( g0(i-4) g0(i-3) g0(i-2) g0(i-1) g0(i) ): " << G0 << endl;
  cout << "Polynom G1       ( g1(i-4) g1(i-3) g1(i-2) g1(i-1) g1(i) ): " << G1 << endl
       << endl;
#endif

  for (vector< bitset<K> >::iterator it = mess.begin(); it != mess.end(); ++it)
  {
    reg = reg << 1;
    reg.set(0, (*it).count());

    g0 = (reg & G0).count() % 2;
    g1 = (reg & G1).count() % 2;

    cod_out.reset();
    cod_out.set(0, g0);
    cod_out.set(1, g1);

    mess_out.push_back(cod_out);

#ifdef DEBUG
    cout << "Block number: " << ++i << " - In frame: " << *it << endl;
    cout << "\t Current status of registers: " << reg << endl;
    cout << "\t Out : " << cod_out << endl;
#endif
  }
#ifdef DEBUG
  cout << "------------------------------------------------------------------" << endl
       << endl;
#endif

  return mess_out;
}

/////////////////////////////////////////////////////////////////////////
// vector< bitset<N> >  GSM_transmission(vector< bitset<N> > mess_cod) //
//                                                                     //
//         Simulation of a transmission channel => adding errors       //
/////////////////////////////////////////////////////////////////////////

// This function generates randoms errors in the transmitted message in
// order to simulate errors on the transmission channel.

vector< bitset<N> > GSM_transmission(vector< bitset<N> > mess_cod)
{
  vector< bitset<N> > mess_tra = mess_cod;

  for (int i = 0; i < NB_ERRORS; i++)
    mess_tra[rand() % mess_cod.size()] = randBitset<N>();

  return mess_tra;
}

//////////////////////////////////////////////////////////////////
// vector< bitset<K> > GSM_decode(vector< bitset<N> > mess_tra) //
//                                                              //
//     Convolutional decoding of a message (GSM norm)           //
//////////////////////////////////////////////////////////////////

vector< bitset<K> > GSM_decode(vector< bitset<N> > mess_tra)
{
  vector< bitset<K> > mess_dec;

  vector< Path > paths;
  vector< Path > new_paths;

  Path initial_path(0000,0); // Initializing the register
  paths.push_back(initial_path);

  vector< bitset<N> >::iterator it;
  for( it=mess_tra.begin(); it!=mess_tra.end(); it++ )
  {

    new_paths = paths;
    paths.clear();

    vector<Path>::iterator it2;
    for( it2=new_paths.begin(); it2!=new_paths.end(); it2++)
    {
      // Creating the two potentials paths
      Path path1 = *it2;
      Path path2 = *it2;

      // Calcute hamming distance
      bitset<R+1> flow_1(path1.reg.to_string()+"0");
      int hdist1 = hamming_distance( code(flow_1), (*it));
      bitset<R+1> flow_2(path2.reg.to_string()+"1");
      int hdist2 = hamming_distance( code(flow_2), (*it));

      // Adding cumulated distance and updating the register
      path1.cumulated_distance+=hdist1;
      path2.cumulated_distance+=hdist2;
      path1.reg=(path1.reg<<1);
      path2.reg=(path2.reg<<1);
      path1.reg[0]=0;
      path2.reg[0]=1;

      path1.mess_dec.push_back(0);
      path2.mess_dec.push_back(1);

      paths.push_back(path1);
      paths.push_back(path2);
    }
  }
  // Sorting on cumulated hamming distance
  std::sort( paths.begin() , paths.end() , weight_comp );
  mess_dec = paths.at(0).mess_dec;
  
  return mess_dec;
}

//////////////////////////////////////////////////////////////////
//                             MAIN                             //
//////////////////////////////////////////////////////////////////

int main()
{
  int NbMot = 12;

  vector< bitset<K> > mess;
  vector< bitset<N> > mess_cod;
  vector< bitset<N> > mess_tra;
  vector< bitset<K> > mess_dec;

  // Random initialization message
  srand((unsigned)time(NULL));
  for (int i = 0; i < NbMot; ++i)
    mess.push_back(randBitset<K>());
  for (int i = 0; i < R; ++i)
    mess.push_back(bitset<K>(0));

  // Coding of the message => mess_cod
  mess_cod = GSM_code(mess);

  // Simulation of a transmission (errors) => mess_tra
  mess_tra = GSM_transmission(mess_cod);

  // Decoding of the transmitted message => mess_dec
  mess_dec = GSM_decode(mess_tra);

  cout << "Source Message   : ";
  for (vector< bitset<K> >::iterator it = mess.begin(); it != mess.end(); ++it)
    cout << ' ' << *it;
  cout << '\n';

  cout << "Coded Message    : ";
  for (vector< bitset<N> >::iterator it = mess_cod.begin(); it != mess_cod.end(); ++it)
    cout << ' ' << *it;
  cout << '\n';

  cout << "Received Message : ";
  for (vector< bitset<N> >::iterator it = mess_tra.begin(); it != mess_tra.end(); ++it)
    cout << ' ' << *it;
  cout << '\n';

  cout << "Decoded Message  : ";
  for (vector< bitset<K> >::iterator it = mess_dec.begin(); it != mess_dec.end(); ++it)
    cout << ' ' << *it;
  cout << '\n';
}
