#ifndef Set_h
#define Set_h
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>
#include "Hash16.h"

#include <cmath> // for log2f
#include <algorithm> // for max/min
#include <vector>

// didn't want to "typedef kmer_type set_elem" because Set is an independent file
#ifdef _largeint
#include "LargeInt.h"
typedef LargeInt<KMER_PRECISION> set_elem;
#else
#ifdef _ttmath
#include "ttmath/ttmath.h"
typedef ttmath::UInt<KMER_PRECISION> set_elem;
#else
#if (! defined kmer_type) || (! defined _LP64)
typedef uint64_t set_elem;
#else
typedef kmer_type set_elem;
#endif
#endif
#endif

using namespace std;

// abstract class
class Set{
    public:
    virtual void insert(set_elem elem) = 0;
    virtual void finalize() = 0;
    virtual bool contains(set_elem elemn) = 0;
};

class HashSet : public Set{

    Hash16 *hash;

    public:
    void insert(set_elem elem);
    void finalize();
    bool contains(set_elem elemn);
    static int bits_per_element ;

    HashSet(uint64_t taille_approx);
};

class ListSet : public Set{

 protected:
    vector<set_elem> liste;

    public:
    void insert(set_elem elem);
    void finalize();
    bool contains(set_elem elemn);
    //Raluca
    bool containsNotBinary(set_elem elem);

    uint64_t capacity() {return (uint64_t)liste.capacity();}
    
    static int bits_per_element ;

    ListSet(uint64_t taille_approx);
    ListSet();

};


//typedef unsigned char set_value_t;
typedef unsigned short int set_value_t; // need 9 bits for Terminator now


class AssocSet : public ListSet {

    vector<set_value_t> liste_value;
 public:
    int get( set_elem elem, set_value_t * val);
    int set(set_elem graine, set_value_t val);
    void finalize();
    AssocSet();
    void print_total_size();
    void clear();

    void start_iterator();
    bool next_iterator();

    
    vector<set_elem>::iterator iterator;



};

//Raluca
typedef struct nt_kmer{
    char nt;
    kmer_type prev_kmer;
}nt_kmer_t;

typedef struct pair_nt_kmer{
    nt_kmer_t nk1, nk2;
}pair_nt_kmer_t;


void copy_nt_kmer (nt_kmer_t from, nt_kmer_t* to);
void copy_pair_nt_kmer (pair_nt_kmer_t from, pair_nt_kmer_t* to);


class AssocPairedSet : public ListSet {
    
    vector<pair_nt_kmer_t> liste_value;
public:
    int get(set_elem elem, pair_nt_kmer_t * val);
    int set(set_elem graine, pair_nt_kmer_t val);
    void finalize();
    AssocPairedSet();
    void print_total_size();
    
    void start_iterator();
    bool next_iterator();
    
    void load_branching(BinaryBank * branches);
    void load(char * filename);
    void dump(char * filename);

    
    vector<set_elem>::iterator iterator;
    
    void print();
    
};

class FPSetCascading4 : public Set{
 public:
  Bloom *bloom2, *bloom3, *bloom4;
  ListSet *false_positives;
  bool contains(set_elem elemn) 
  { 
    if (bloom2->contains(elemn))
    {
      if (!bloom3->contains(elemn)) 
	return true;
      else if (bloom4->contains(elemn) && !false_positives->contains(elemn)) 
	return true;
    }
    return false;
  };
  void insert(set_elem elem) {fprintf (stderr, "Error can't insert in FPSetCascading!\n"); exit(0); };
  void finalize() {fprintf (stderr, "Error can't finalize in FPSetCascading!\n"); exit(0);};
  bool is_false_positive(set_elem elemn) {return contains(elemn);};
};

#endif
