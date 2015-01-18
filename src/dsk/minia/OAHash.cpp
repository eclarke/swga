// open-addressing hash table with linear probing, follows wikipedia
// to reduce memory, elements with [value == 0] are UNOCCUPIED, deal with it.

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <algorithm> // for max

#include "OAHash.h"

using namespace::std;


OAHash::OAHash(uint64_t max_memory) // in bytes
{
    hash_size = max_memory / sizeof(element_pair);
    if (hash_size == 0)
    {
        printf("empty OAHash allocated\n");
        exit(1);
    }
    nb_inserted_keys = 0;
    data = (element_pair *) calloc( hash_size, sizeof(element_pair));  //create hashtable
}

OAHash::~OAHash()
{
    free(data);
}

int OAHash::size_entry()
{
    return sizeof(element_pair);
}

// hash functions: [ any integer type, e.g. 64 bits, 128 bits or ttmath ] -> [ 64 bits hash ]
#ifdef _largeint
inline uint64_t OAHash::hashcode(LargeInt<KMER_PRECISION> elem)
{
    // hash = XOR_of_series[hash(i-th chunk iof 64 bits)]
    uint64_t result = 0, chunk, mask = ~0;
    LargeInt<KMER_PRECISION> intermediate = elem;
    int i;
    for (i=0;i<KMER_PRECISION;i++)
    {
        chunk = (intermediate & mask).toInt();
        intermediate = intermediate >> 64;
        result ^= hashcode(chunk);
    }
    return result;
}
#endif

#ifdef _ttmath
inline uint64_t OAHash::hashcode(ttmath::UInt<KMER_PRECISION> elem)
{
    // hash = XOR_of_series[hash(i-th chunk iof 64 bits)
    uint64_t result = 0, to_hash;
    ttmath::UInt<KMER_PRECISION> intermediate = elem;
    uint32_t cmask=~0, chunk;
    int i;
    for (i=0;i<KMER_PRECISION/2;i++)
    {
        // retrieve a 64 bits part to hash 
        (intermediate & cmask).ToInt(chunk);
        to_hash = chunk;
        intermediate >>= 32;
        (intermediate & cmask).ToInt(chunk);
        to_hash |= ((uint64_t)chunk) << 32 ;
        intermediate >>= 32;

        result ^= hashcode(to_hash);
    }
    return result;
}
#endif

#ifdef _LP64
inline uint64_t OAHash::hashcode( __uint128_t elem )
{
    // hashcode(uint128) = ( hashcode(upper 64 bits) xor hashcode(lower 64 bits))
    return (hashcode((uint64_t)(elem>>64)) ^ hashcode((uint64_t)(elem&((((__uint128_t)1)<<64)-1))));
}
#endif

inline uint64_t OAHash::hashcode( uint64_t elem )
{
    uint64_t code = elem;
    
    code = code ^ (code >> 14); //supp
    code = (~code) + (code << 18); 
    code = code ^ (code >> 31);
    code = code * 21; 
    code = code ^ (code >> 11);
    code = code + (code << 6);
    code = code ^ (code >> 22);

    return code;
    
}

bool OAHash::is_occupied(element_pair *element)
{
    return (element->value != 0);
}

OAHash::element_pair * OAHash::find_slot(key_type key)
{
    uint64_t ptr = hashcode(key) % hash_size; 
    element_pair *element = data+ptr;
    uint64_t retries = 0;

    // search until we either find the key, or find an empty slot.
    while ( ( is_occupied(element)) && ( element->key != key ) && (retries < hash_size))
    {
        ptr = (ptr + 1) % hash_size;
        element = data+ptr;
        retries++;
    }
    if (retries == hash_size)
    {
        printf("OAHash: max rehashes reached: %lld (notify a developer)\n",(long long)hash_size);
        exit(1);
    }

    return element;
}

//if graine already here, overwrite old value
void OAHash::insert(key_type graine, int value)
{
    element_pair *element = find_slot(graine);
    if (!is_occupied(element))
    {
        element->key = graine;
        nb_inserted_keys++;
    }
    element->value = value;
}

// increment the value of a graine
void OAHash::increment(key_type graine)
{
    element_pair *element = find_slot(graine);
    if (!is_occupied(element))
    {
        element->key = graine;
        nb_inserted_keys++;
    }
    if( element->value == -1)  element->value = 0; //special case, emulate 0 value with -1, (0 is not a valid value, used for empty cell)
    
    element->value = element->value + 1;
}

bool OAHash::get( key_type graine, int * val)
{ 
    element_pair *element = find_slot(graine);
    if (!is_occupied(element))
        return false;
    if ((element->key) == graine && (val != NULL))
        *val = element->value;
    
     if( element->value ==-1) *val = 0; // 0 is emulated with -1
    return true;
}

bool OAHash::has_key(key_type graine)
{
    return get(graine,NULL);
}


// call start_iterator to reinit the iterator, then do a while(next_iterator()) {..} to traverse every cell
void OAHash::start_iterator()
{
    iterator = data-1;
}


// returns true as long as the iterator contains a valid cell
bool OAHash::next_iterator()
{
    while (1)
    {
        iterator++;
        if (iterator == data+hash_size)
            return false;
        if (iterator->value != 0)
            break;
    }
    return true;
}


float OAHash::load_factor()
{
    return (float)nb_inserted_keys/(float)hash_size;
}


uint64_t OAHash::memory_usage()
{
    return hash_size* sizeof(element_pair); // in bits
}

void OAHash::printstat()
{
    fprintf(stderr,"\n----------------------Stat OA Hash Table ---------------------\n");
    
    fprintf(stderr,"max elements: %lld, memory usage: %lld\n",(long long)hash_size,(long long)memory_usage());
    fprintf(stderr,"load factor: %.2f\n",load_factor());
}   
