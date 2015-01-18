#include "Set.h"


int ListSet::bits_per_element = sizeof(set_elem)*8;

ListSet::ListSet(uint64_t taille_approx)
{
  liste.reserve(taille_approx);
}
ListSet::ListSet()
{
}

void ListSet::insert(set_elem elem)
{
    liste.push_back(elem);
}

void ListSet::finalize()
{
    sort(liste.begin(), liste.end());
}

bool ListSet::contains(set_elem elem)
{
    return binary_search(liste.begin(), liste.end(), elem);
}

//Raluca
bool ListSet::containsNotBinary(set_elem elem)
{
    int i;
    for (i=0; i<liste.size(); i++)
        if (liste[i]==elem)
            return true;
    
    return false;
}

//--------------
// plain old expensive hash

int HashSet::bits_per_element = sizeof(cell<hash_elem>)*8;

HashSet::HashSet(uint64_t taille_approx)
{
    int NBITS_HT = max( (int)ceilf(log2f(taille_approx)) , 1);
    hash = new Hash16(NBITS_HT);
}

void HashSet::insert(set_elem elem)
{
    hash->insert(elem,1); // dummy value
}

void HashSet::finalize()
{
    // NOP
}

bool HashSet::contains(set_elem elem)
{
    return hash->has_key(elem);
}



//--------------
// emulates a hash table with two lists
// benefits: lower memory usage
// but logarithmic acess instead of constant


void AssocSet::finalize()
{
    sort(liste.begin(), liste.end());
    liste_value.assign(liste.size(),0); 
}


int AssocSet::get( set_elem elem, set_value_t * val)
{

  vector<set_elem>::iterator it;

  it = lower_bound(liste.begin(), liste.end(),elem); 
  if (it == liste.end() || elem != *it) return 0;

  size_t rank = it - liste.begin();

  *val = liste_value[rank];
  return 1;
}

AssocSet::AssocSet()
{
}

//return -1 if elem is not in the set
int AssocSet::set(set_elem elem, set_value_t val)
{
 vector<set_elem>::iterator it;

  it = lower_bound(liste.begin(), liste.end(),elem); 
  if (it == liste.end() ||elem != *it) return 0;

  size_t rank = it - liste.begin();

  liste_value[rank]=val;

  return 1;
}


void AssocSet::start_iterator()
{
  iterator = liste.begin()-1;

}

bool AssocSet::next_iterator()
{
  iterator++;
  if (iterator==liste.end()) return false;

  return true;

}

void AssocSet::print_total_size()
{
    printf("Assoc set size: %li\n",liste.size());
    printf("Assoc set capacity: listekmer %li  liste val%li\n",liste.capacity(),liste_value.capacity());
    printf("%li *%li  + %li* %li =  %li MB \n",liste.capacity(),sizeof(set_elem),liste_value.capacity(),sizeof(set_value_t),
            (liste.capacity()*sizeof(set_elem)+liste_value.capacity()*sizeof(set_value_t))/1024/1024
          );
}

void AssocSet::clear()
{
    liste_value.assign(liste_value.size(),0);
//    liste_value.clear();
    
}


//Raluca
//--------------
// emulates a hash table with two lists for paired reads (elements=branching kmers in the second read file, values=pairs of corresponding kmers in the first read file + next nucleotide in the second read file)

void copy_nt_kmer (nt_kmer_t from, nt_kmer_t* to){
    to->nt = from.nt;
    to->prev_kmer = from.prev_kmer;
}


void copy_pair_nt_kmer (pair_nt_kmer_t from, pair_nt_kmer_t* to){
    copy_nt_kmer(from.nk1, &(to->nk1));
    copy_nt_kmer(from.nk2, &(to->nk2));
}


void AssocPairedSet::finalize()
{

    sort(liste.begin(), liste.end());
    //liste_value.assign(liste.size(),0);
    liste_value.reserve(liste.size());
  //  printf("finalize %i elems \n",liste_value.size());

    liste_value.resize(liste.size());


    int i;
    for (i=0; i<liste_value.size(); i++)
    {
        liste_value[i].nk1.nt = 0;
        liste_value[i].nk2.nt = 0;
    }
}


int AssocPairedSet::get( set_elem elem, pair_nt_kmer_t * val)
{
    vector<set_elem>::iterator it;
    it = lower_bound(liste.begin(), liste.end(),elem);
    if (it == liste.end() || elem != *it) return 0;
    
    size_t rank = it - liste.begin();
    
    copy_pair_nt_kmer(liste_value[rank], val);
    
    return 1;
}



AssocPairedSet::AssocPairedSet()
{
}


//return -1 if elem is not in the set
int AssocPairedSet::set(set_elem elem, pair_nt_kmer_t val)
{
    vector<set_elem>::iterator it;
    
    it = lower_bound(liste.begin(), liste.end(),elem);
    if (it == liste.end() ||elem != *it) return 0;
    
    size_t rank = it - liste.begin();
    
    copy_pair_nt_kmer(val, &liste_value[rank]);
    
    return 1;
}



void AssocPairedSet::start_iterator()
{
    iterator = liste.begin()-1;
    
}


bool AssocPairedSet::next_iterator()
{
    iterator++;
    if (iterator==liste.end()) return false;
    
    return true;
    
}


void AssocPairedSet::print()
{
    int i;
    char seq[100];
    
    printf("print %lu elems \n",liste.size());

    
    for (i=0; i<liste.size(); i++)
    {
        code2seq(liste[i], seq);
        printf("%s=(", seq);
        
        code2seq(liste_value[i].nk1.prev_kmer, seq);
        printf("%s,%c) (", seq, liste_value[i].nk1.nt);
        
        code2seq(liste_value[i].nk2.prev_kmer, seq);
        printf("%s,%c)\n", seq, liste_value[i].nk2.nt);
    }
}


void AssocPairedSet::load_branching(BinaryBank * branches)
{
    // load branching kmers
    branches->rewind_all();
    kmer_type kmer;
    while (branches->read_element(&kmer))
        liste.push_back(kmer);
    
    printf("done loaded  branching kmers in AssocPairedSet \n");
}


void AssocPairedSet::dump(char * filename)
{
    FILE *file_data;
    file_data = fopen(filename,"wb");
    fwrite(&liste_value[0], sizeof(pair_nt_kmer_t), liste_value.size(), file_data); //1+
    printf("liste_value dumped \n");
    
}


void AssocPairedSet::load(char * filename)
{
    FILE *file_data;
    file_data = fopen(filename,"rb");
    printf("loading assoc paired kmer index from file\n");
    
    liste_value.reserve(liste.size());
    liste_value.resize(liste.size());
    
    fread(&liste_value[0], sizeof(pair_nt_kmer_t), liste_value.size(), file_data);
    printf("assoc paired loaded\n");

}


void AssocPairedSet::print_total_size()
{
    printf("Assoc paired set size: %li\n",liste.size());
    printf("Assoc paired set capacity: listekmer %li  liste val%li\n",liste.capacity(),liste_value.capacity());
    printf("%li *%li  + %li* %li =  %li MB \n",liste.capacity(),sizeof(set_elem),liste_value.capacity(),sizeof(set_value_t),
           (liste.capacity()*sizeof(set_elem)+liste_value.capacity()*sizeof(set_value_t))/1024/1024
           );
    
}

