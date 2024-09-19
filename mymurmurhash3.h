// MurmurHash3 was written by Austin Appleby.
// Imai modified it a little in this code.

#include <stdint.h>

void mymurmurhash3(uint32_t key, uint32_t seed, uint32_t * out){
  uint32_t h = seed;
  uint32_t k = key;

  k *= 0xcc9e2d51;
  k = (k << 15) | (k >> 17);
  k *= 0x1b873593;
  h ^= k;
  h = (h << 13) | (h >> 19);
  h = h*5 + 0xe6546b64;

  h ^= 4;
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  *out = h;
}

