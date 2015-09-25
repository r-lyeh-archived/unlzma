/* A very compact LZMA decoder (C++03). Original lzd.cc code by Antonio Diaz Diaz.
   - rlyeh, public domain 

   Features:
   [x] Original lzd.cc code by Antonio Diaz Diaz.
   [x] Simple. Public API is two functions: unlzma() and unlzma_size().
   [x] Portable, header-only and self-contained.
   [x] Tiny (~320 LOC).

   Cons:
   [ ] ~4x times slower than official LZMA SDK (https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv%E2%80%93Markov_chain_algorithm)
   [ ] uncommon LZMA streams with props.lp flag <> 0 are not supported 

   Public API:
   // Unpacks LZMA stream.
   // Destination buffer must be allocated in advance (see function below)
   // Returns: human error string or empty string if ok.
   const char *const unlzma( *dst, dstlen, const *src, srclen );

   // Retrieves unpacked size of LZMA stream
   // Returns: unpacked size or 0 if invalid
   uint64_t unlzma_size( src, srclen );
*/

#pragma once
#include <stdint.h>
#include <cstring>
#include <algorithm>
#include <vector>

#define UNLZMA_VERSION "1.0.0" // (2015/09/20) Initial version

namespace {
enum {
    len_states = 4,
    dis_slot_bits = 6,
    start_dis_model = 4,
    end_dis_model = 14,
    dis_modeled = 1 << (end_dis_model / 2), // fb 128
    dis_align_bits = 4,
    dis_align_size = 1 << dis_align_bits,

    len_low_bits = 3,
    len_mid_bits = 3,
    len_high_bits = 8,
    len_low_symbols = 1 << len_low_bits,
    len_mid_symbols = 1 << len_mid_bits,
    len_high_symbols = 1 << len_high_bits,
    max_len_symbols = len_low_symbols + len_mid_symbols + len_high_symbols,

    min_match_len = 2,

    bit_model_move_bits = 5,
    bit_model_total_bits = 11,
    bit_model_total = 1 << bit_model_total_bits,

    num_states = 12
};

struct Bit_model {
    int probability;
    Bit_model() : probability( bit_model_total / 2 ) {}
};

struct io {
    const uint8_t *in0, *in1;
          uint8_t *out0, *out1;
};

struct State {
    int st;
    State() : st( 0 ) {}
    int  operator()() const { return st; }
    bool is_char()    const { return st < 7; }
    void set_match()        { st = ( st < 7 ) ? 7 : 10; }
    void set_rep()          { st = ( st < 7 ) ? 8 : 11; }
    void set_short_rep()    { st = ( st < 7 ) ? 9 : 11; }
    void set_char() {
        const int next[num_states] = { 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 4, 5 };
        st = next[st];
    }
};

struct Len_model {
    Bit_model choice1;
    Bit_model choice2;
    Bit_model bm_high[len_high_symbols];
    std::vector< std::vector<Bit_model> > bm_low, bm_mid;
    explicit Len_model( int pos_states ) : 
        bm_low( std::vector< std::vector<Bit_model> >( pos_states, std::vector<Bit_model>(len_low_symbols) ) ),
        bm_mid( std::vector< std::vector<Bit_model> >( pos_states, std::vector<Bit_model>(len_mid_symbols) ) )
    {}
};

struct Range_decoder {
    io &mem;
    uint32_t code, range;

    explicit Range_decoder(io &mem) : mem(mem), code( 0 ), range( 0xFFFFFFFFU ) {
        for( int i = 0; i < 5; ++i ) code = (code << 8) | get_byte();
    }

    uint8_t get_byte() {
        return *mem.in0++;
    }

    int decode( const int num_bits ) {
        int symbol = 0;
        for( int i = 0; i < num_bits; ++i ) {
            range >>= 1;
            symbol <<= 1;
            if( code >= range ) { code -= range; symbol |= 1; }
            if( range <= 0x00FFFFFFU ) { // normalize
                range <<= 8; code = (code << 8) | get_byte();
            }
        }
        return symbol;
    }

    int decode_bit( Bit_model &bm ) {
        int symbol;
        const uint32_t bound = ( range >> bit_model_total_bits ) * bm.probability;
        if( code < bound ) {
            range = bound;
            bm.probability += (bit_model_total - bm.probability) >> bit_model_move_bits;
            symbol = 0;
        } else {
            range -= bound;
            code -= bound;
            bm.probability -= bm.probability >> bit_model_move_bits;
            symbol = 1;
        }
        if( range <= 0x00FFFFFFU ) { // normalize
            range <<= 8; code = (code << 8) | get_byte();
        }
        return symbol;
    }

    int decode_tree( Bit_model bm[], const int num_bits ) {
        int symbol = 1;
        for( int i = 0; i < num_bits; ++i )
            symbol = ( symbol << 1 ) | decode_bit( bm[symbol] );
        return symbol - (1 << num_bits);
    }

    int decode_tree_reversed( Bit_model bm[], const int num_bits ) {
        int symbol = decode_tree( bm, num_bits );
        int reversed_symbol = 0;
        for( int i = 0; i < num_bits; ++i ) {
            reversed_symbol = ( reversed_symbol << 1 ) | ( symbol & 1 );
            symbol >>= 1;
        }
        return reversed_symbol;
    }

    int decode_matched( Bit_model bm[], const int match_byte ) {
        Bit_model * const bm1 = bm + 0x100;
        int symbol = 1;
        for( int i = 7; i >= 0; --i ) {
            const int match_bit = ( match_byte >> i ) & 1;
            const int bit = decode_bit( bm1[(match_bit<<8)+symbol] );
            symbol = ( symbol << 1 ) | bit;
            if( match_bit != bit ) {
                while( symbol < 0x100 )
                symbol = ( symbol << 1 ) | decode_bit( bm[symbol] );
                break;
            }
        }
        return symbol & 0xFF;
    }

    int decode_len( Len_model &lm, const int pos_state ) {
        if( decode_bit( lm.choice1 ) == 0 )
            return decode_tree( lm.bm_low[pos_state].data(), len_low_bits );
        if( decode_bit( lm.choice2 ) == 0 )
            return len_low_symbols +
        decode_tree( lm.bm_mid[pos_state].data(), len_mid_bits );
        return len_low_symbols + len_mid_symbols + decode_tree( lm.bm_high, len_high_bits );
    }
};

struct LZ_decoder {
    io &mem;
    Range_decoder rdec;
    uint64_t partial_data_pos;
    const unsigned dictionary_size;
    std::vector<uint8_t> buffer;    // output buffer
    unsigned pos, stream_pos;       // current pos in buffer, first unwritten byte, dict size
    const int literal_context_bits, pos_state_bits; // lc [0..8] def(3), pb [0..4] def(2)
    const int pos_states, pos_state_mask;

    void flush_data() {
        if( pos > stream_pos ) {
            unsigned size = (std::min)(pos - stream_pos, (unsigned)(mem.out1 - mem.out0));
            memcpy( mem.out0, buffer.data() + stream_pos, size );
            mem.out0 += size;
            if( pos >= dictionary_size ) {
                partial_data_pos += pos;
                pos = 0;
            }
            stream_pos = pos;
        }
    }

    uint8_t peek( const unsigned distance ) const {
        unsigned i = pos - distance - 1;
        return buffer[ pos > distance ? i : i + dictionary_size ];
    }

    void put_byte( const uint8_t b ) {
        buffer[pos] = b;
        if( ++pos >= dictionary_size ) flush_data();
    }

    uint64_t data_position() const {
        return partial_data_pos + pos;
    }

    LZ_decoder( io &mem, const unsigned dictionary_size, int lc, int pb ) :
        mem(mem),
        rdec(mem),
        partial_data_pos( 0 ),
        dictionary_size( dictionary_size ),
        buffer( dictionary_size ),
        pos( 0 ),
        stream_pos( 0 ),
        literal_context_bits( lc), 
        pos_state_bits( pb), 
        pos_states( 1 << pos_state_bits), 
        pos_state_mask( pos_states - 1) {
            buffer[dictionary_size-1] = 0;  // prev_byte of first byte
        }

    int operator()() {  // Returns line error or 0 if ok
        State state;
        unsigned rep[4] = {}; // rep[0-3] latest four distances used for efficient coding of repeated distances
        Bit_model bm_rep[4][num_states], bm_dis[dis_modeled-end_dis_model], bm_align[dis_align_size];
        Len_model match_len_model(pos_states), rep_len_model(pos_states);
        std::vector< std::vector<Bit_model> > 
            bm_literal( 1<<literal_context_bits, std::vector<Bit_model>(0x300) ),
            bm_match( num_states, std::vector<Bit_model>(pos_states) ),
            bm_len( num_states, std::vector<Bit_model>(pos_states ) ),
            bm_dis_slot( len_states, std::vector<Bit_model>(1<<dis_slot_bits) );

        while( mem.in0 < mem.in1 ) {
            const int pos_state = data_position() & pos_state_mask;
            if( rdec.decode_bit( bm_match[state()][pos_state] ) == 0 ) { // bit1
                const uint8_t prev_byte = peek( 0 );
                const int literal_state = prev_byte >> ( 8 - literal_context_bits );
                Bit_model * const bm = bm_literal[literal_state].data();
                put_byte( state.is_char() ? rdec.decode_tree( bm, 8 ) : rdec.decode_matched( bm, peek( rep[0] ) ) );
                state.set_char();
            } else { // match or repeated match
                int len;
                if( rdec.decode_bit( bm_rep[0][state()] ) != 0 ) { // bit2
                    if( rdec.decode_bit( bm_rep[1][state()] ) != 0 ) { // bit3
                        unsigned distance;
                        if( rdec.decode_bit( bm_rep[2][state()] ) == 0 ) { // bit4
                            distance = rep[1];
                        } else {
                            if( rdec.decode_bit( bm_rep[3][state()] ) == 0 ) { // bit5
                                distance = rep[2];
                            } else {
                                distance = rep[3];
                                rep[3] = rep[2];
                            }
                            rep[2] = rep[1]; 
                        }
                        rep[1] = rep[0];
                        rep[0] = distance;
                    } else {
                        if( rdec.decode_bit( bm_len[state()][pos_state] ) == 0 ) { // bit4
                            state.set_short_rep(); 
                            put_byte( peek( rep[0] ) );
                            continue;
                        }
                    }
                    state.set_rep();
                    len = min_match_len + rdec.decode_len( rep_len_model, pos_state );
                } else { // match
                    rep[3] = rep[2]; rep[2] = rep[1]; rep[1] = rep[0];
                    len = min_match_len + rdec.decode_len( match_len_model, pos_state );
                    const int len_state = (std::min)( len - min_match_len, len_states - 1 );
                    const int dis_slot = rdec.decode_tree( bm_dis_slot[len_state].data(), dis_slot_bits );
                    if( dis_slot >= start_dis_model ) {
                        const int direct_bits = ( dis_slot >> 1 ) - 1;
                        rep[0] = ( 2 | ( dis_slot & 1 ) ) << direct_bits;
                        if( dis_slot >= end_dis_model ) {
                            rep[0] += rdec.decode( direct_bits - dis_align_bits ) << dis_align_bits;
                            rep[0] += rdec.decode_tree_reversed( bm_align, dis_align_bits );
                            if( rep[0] == 0xFFFFFFFFU ) { // End Of Stream marker found
                                return flush_data(), __LINE__ * !(len == min_match_len && mem.out0 == mem.out1); 
                            }
                        } else rep[0] += rdec.decode_tree_reversed( bm_dis + rep[0] - dis_slot - 1, direct_bits );
                    } else rep[0] = dis_slot;
                    state.set_match();
                    if( rep[0] >= dictionary_size || rep[0] >= data_position() ) {
                        return flush_data(), __LINE__;
                    }
                }
                for( int i = 0; i < len; ++i ) {
                    put_byte( peek( rep[0] ) );
                }
            }
        }
        return flush_data(), __LINE__ * !(mem.out0 == mem.out1);
    }
};
}

static inline uint64_t unlzma_size( const void *src, uint64_t ilen ) {
    uint64_t data_size = 0;
    if( src && ilen > 13 ) {
        for( int i = 0; i < 8; i++ ) {
            data_size += *(((uint8_t*)src) + 5 + i) << (i * 8);
        }
    }
    return data_size;
}

static inline const char *const unlzma( void *dst, uint64_t olen, const void *src, uint64_t ilen ) {
    io mem = { (const uint8_t*)src, (const uint8_t*)src+ilen, (uint8_t*)dst, (uint8_t*)dst+olen };
    if( mem.in0 && mem.out0 && mem.in0 < mem.in1 && mem.out0 < mem.out1 ) {
        uint8_t lc, lp, pb, prop = *mem.in0++;
        if( prop <= (4 * 5 + 4) * 9 + 8 ) {
            pb = (prop) / (9 * 5);
            lp = (prop - pb * 9 * 5) / 9;
            lc = (prop - pb * 9 * 5) - lp * 9;
            if( lp == 0 ) {
                uint64_t dict = mem.in0[0] | (mem.in0[1] << 8) | (mem.in0[2] << 16) | (mem.in0[3] << 24);
                mem.in0 += 12;
                int errline = LZ_decoder( mem, dict, lc, pb )();
                return errline ? "Data error" : "";
            } else return "Unsupported props.lp flag";
        } else return "Unsupported props byte";
    } else return "Bad pointer";
}

template<typename T, typename U>
static inline const char *const unlzma( T &t, const U &u ) {
    t.resize( unlzma_size( u.data(), u.size() ) );
    return unlzma( &t[0], t.size(), u.data(), u.size() );
}

#ifdef UNLZMA_BUILD_SAMPLE
#include <iostream>
#include <fstream>
#include <sstream>

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
#endif
