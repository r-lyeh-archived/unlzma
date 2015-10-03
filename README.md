unlzma 
====== 

UNLZMA is a very compact LZMA decoder (C++03).

## Features
- [x] Original lzd.cc code by Antonio Diaz Diaz.
- [x] Simple. Public API is two functions: unlzma() and unlzma_size().
- [x] Portable, header-only and self-contained.
- [x] Tiny (~320 LOC).
- [x] Public domain.

## Cons
- [ ] Four times slower than [official LZMA SDK](https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv%E2%80%93Markov_chain_algorithm).
- [ ] Uncommon LZMA streams with props.lp flag <> 0 are not supported.

## Public API
```c++
   // Unpacks LZMA stream.
   // Destination buffer must be allocated in advance (see function below)
   // Returns: human error string or "\0" if ok.
   const char * const unlzma( *dst, dstlen, const *src, srclen );

   // Retrieves unpacked size of LZMA stream
   // Returns: unpacked size or 0 if invalid
   uint64_t unlzma_size( src, srclen );
```

## Showcase
```c++
#include <iostream>
#include <fstream>
#include <sstream>
#include "unlzma.hpp"

int main( const int argc, const char **argv ) {
    if( argc != 3 ) {
        std::cout << argv[0] << " infile.lzma outfile" << std::endl;
    } else {
        std::ifstream ifs( argv[1], std::ios::binary );
        std::ofstream ofs( argv[2], std::ios::binary );
        if( ifs.good() && ofs.good() ) {
            std::stringstream ss; ss << ifs.rdbuf();
            std::string in = ss.str();
            std::string out;

            const char *errors = unlzma(out, in);
            std::cout << errors << std::endl;
            ofs << out;
        }
    }
}
```

## Changelog
- v1.0.0 (2015/09/20)
  - Initial version
